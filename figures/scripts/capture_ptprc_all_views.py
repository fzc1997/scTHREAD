#!/usr/bin/env python3
"""Capture consistent audit crops for all six PTPRC browser evidence views."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from playwright.sync_api import Page, sync_playwright


VIEWS = (
    ("junctions", "01_junctions.png"),
    ("cellmap", "02_cell_map.png"),
    ("expression", "03_expression.png"),
    ("isoforms", "04_isoforms.png"),
    ("pas", "05_polya.png"),
    ("ase", "06_ase.png"),
)


def wait_for_gene(page: Page, gene: str) -> None:
    page.wait_for_selector("#browserSearch", state="visible")
    page.wait_for_function(
        """gene => {
            const location = document.querySelector("#browserLocation")?.textContent || "";
            const evidence = document.querySelector("#geneEvidence")?.textContent || "";
            return location.includes(gene) && !location.toLowerCase().includes("loading")
                && evidence.includes("junctions") && !evidence.includes("Loading");
        }""",
        arg=gene,
        timeout=60_000,
    )


def wait_for_view(page: Page, view: str) -> None:
    page.wait_for_function(
        """view => {
            const panel = document.querySelector(`#${view}View`);
            if (!panel?.classList.contains("active")) return false;
            if (view === "junctions") {
                return document.querySelectorAll("#browserTableBody tr").length > 0;
            }
            if (view === "cellmap") {
                const legend = document.querySelector("#umapSignalLegend");
                const canvas = document.querySelector("#umapSignalCanvas");
                return Boolean(
                    legend?.dataset?.loadState?.startsWith("rendered:")
                    &&
                    legend?.textContent?.trim()
                    && !legend.textContent.includes("Loading")
                    && canvas
                    && canvas.width > 0
                    && canvas.height > 0
                );
            }
            const output = document.querySelector(`#${view}Output`);
            if (!output || output.textContent.includes("Loading")) return false;
            return Boolean(
                output.querySelector(".track-row")
                || output.querySelector(".empty-state")
                || output.querySelector(".error-state")
            );
        }""",
        arg=view,
        timeout=60_000,
    )
    page.wait_for_timeout(700)


def _select_with_option(
    page: Page,
    *,
    option_text: str,
    startswith: bool = False,
) -> dict[str, object]:
    """Select the first Cell map control containing the requested option."""
    controls = page.locator("#cellmapView select")
    for control_index in range(controls.count()):
        control = controls.nth(control_index)
        options = control.locator("option")
        for option_index, text in enumerate(options.all_inner_texts()):
            normalized = text.strip()
            matched = (
                normalized.startswith(option_text)
                if startswith
                else normalized == option_text
            )
            if not matched:
                continue
            option = options.nth(option_index)
            value = option.get_attribute("value")
            if value is None:
                raise RuntimeError(f"Option {normalized!r} has no value")
            control.select_option(value=value)
            return {
                "control_index": control_index,
                "control_id": control.get_attribute("id"),
                "selected_value": value,
                "selected_text": normalized,
            }
    raise RuntimeError(f"Could not find Cell map option {option_text!r}")


def configure_cellmap(
    page: Page,
    *,
    right_view: str,
    isoform: str,
    source: str,
) -> dict[str, object]:
    """Configure the right UMAP before capture and wait for its rerender."""
    view_labels = {
        "gene": "Gene expression",
        "isoform": "Isoform expression",
        "ase": "Allelic balance (ASE)",
    }
    mode_selection = _select_with_option(
        page,
        option_text=view_labels[right_view],
    )
    isoform_selection = None
    if right_view == "isoform":
        isoform_selection = _select_with_option(
            page,
            option_text=isoform,
            startswith=True,
        )
    source_selection = None
    if source:
        source_selection = _select_with_option(
            page,
            option_text=source,
            startswith=True,
        )
    page.wait_for_function(
        """expected => {
            const legend = document.querySelector("#umapSignalLegend");
            const canvas = document.querySelector("#umapSignalCanvas");
            const state = legend?.dataset?.loadState || "";
            return state.startsWith("rendered:")
                && (!expected.isoform || document.querySelector("#umapIsoform")?.value === expected.isoform)
                && (!expected.source || document.querySelector("#umapStudy")?.value === expected.source)
                && Boolean(legend?.textContent?.trim())
                && canvas
                && canvas.width > 0
                && canvas.height > 0;
        }""",
        arg={"isoform": isoform if right_view == "isoform" else "", "source": source},
        timeout=60_000,
    )
    page.wait_for_timeout(1_200)
    return {
        "right_view": right_view,
        "mode_selection": mode_selection,
        "isoform_selection": isoform_selection,
        "source_selection": source_selection,
        "load_state": page.locator("#umapSignalLegend").get_attribute(
            "data-load-state"
        ),
        "legend_text": page.locator("#umapSignalLegend").inner_text().strip(),
    }


def audit_clip(page: Page, view: str, max_height: float) -> dict[str, float]:
    workspace = page.locator("#transcriptBrowser .browser-workspace").bounding_box()
    tabs = page.locator(".browser-view-tabs").bounding_box()
    panel = page.locator(f"#{view}View").bounding_box()
    if not workspace or not tabs or not panel:
        raise RuntimeError(f"Could not resolve crop geometry for {view}")
    pad = 2.0
    top = max(0.0, tabs["y"] - pad)
    natural_bottom = panel["y"] + panel["height"] + pad
    bottom = min(natural_bottom, top + max_height)
    return {
        "x": max(0.0, workspace["x"] - pad),
        "y": top,
        "width": workspace["width"] + 2 * pad,
        "height": bottom - top,
    }


def view_summary(page: Page, view: str) -> dict[str, object]:
    panel = page.locator(f"#{view}View")
    text = panel.inner_text().strip()
    return {
        "view": view,
        "active_button": page.locator(
            f'button[data-browser-view="{view}"]'
        ).inner_text(),
        "junction_rows": page.locator("#browserTableBody tr").count()
        if view == "junctions"
        else None,
        "track_rows": panel.locator(".track-row").count(),
        "text_excerpt": text[:4_000],
    }


def export_junction_table(page: Page, output_dir: Path) -> tuple[Path, int]:
    """Export every currently rendered junction row as a traceable TSV."""
    headers = ("junction", "span", "molecules", "reads", "runs", "studies")
    rows: list[list[str]] = []
    table_rows = page.locator("#browserTableBody tr")
    for row_index in range(table_rows.count()):
        cells = [
            text.strip()
            for text in table_rows.nth(row_index).locator("td").all_inner_texts()
        ]
        if len(cells) != len(headers):
            raise RuntimeError(
                f"Unexpected junction row {row_index + 1}: "
                f"{len(cells)} cells, expected {len(headers)}"
            )
        rows.append(cells)
    if not rows:
        raise RuntimeError("No rendered junction rows were available for export")

    output_path = output_dir / "junction_first100_live.tsv"
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(rows)
    return output_path, len(rows)


def capture(args: argparse.Namespace) -> None:
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    page_url = (
        f"{args.base_url.rstrip('/')}/browse"
        f"?species={args.species}&query={args.gene}"
    )

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=args.chrome,
            headless=True,
            args=["--disable-gpu", "--force-device-scale-factor=1"],
        )
        context = browser.new_context(
            viewport={
                "width": args.viewport_width,
                "height": args.viewport_height,
            },
            device_scale_factor=args.device_scale_factor,
            color_scheme="light",
            reduced_motion="reduce",
        )
        page = context.new_page()
        page.goto(page_url, wait_until="domcontentloaded", timeout=60_000)
        wait_for_gene(page, args.gene)
        page.add_style_tag(
            content="""
                .browser-shell { grid-template-columns: minmax(0, 1fr) !important; }
                .browser-sidebar { display: none !important; }
                .browser-workspace { min-width: 0 !important; }
                * { animation: none !important; transition: none !important; }
            """
        )

        selected_views = (
            tuple(item for item in VIEWS if item[0] == args.only_view)
            if args.only_view
            else VIEWS
        )
        summaries: list[dict[str, object]] = []
        for view, filename in selected_views:
            page.click(f'button[data-browser-view="{view}"]')
            wait_for_view(page, view)
            cellmap_config = None
            if view == "cellmap":
                cellmap_config = configure_cellmap(
                    page,
                    right_view=args.cellmap_right_view,
                    isoform=args.cellmap_isoform,
                    source=args.cellmap_source,
                )
            output_path = output_dir / filename
            page.screenshot(
                path=str(output_path),
                clip=audit_clip(page, view, args.max_crop_height),
                animations="disabled",
            )
            summary = view_summary(page, view)
            summary["file"] = filename
            summary["png_bytes"] = output_path.stat().st_size
            if view == "junctions":
                table_path, table_rows = export_junction_table(page, output_dir)
                summary["junction_table_file"] = table_path.name
                summary["junction_table_rows"] = table_rows
            if cellmap_config is not None:
                summary["cellmap_config"] = cellmap_config
            summaries.append(summary)

        metadata = {
            "captured_at_utc": datetime.now(timezone.utc).isoformat(),
            "source_url": page.url,
            "gene": args.gene,
            "species": args.species,
            "viewport_css_px": [args.viewport_width, args.viewport_height],
            "device_scale_factor": args.device_scale_factor,
            "max_crop_height_css_px": args.max_crop_height,
            "only_view": args.only_view,
            "cellmap_right_view": args.cellmap_right_view,
            "cellmap_isoform": args.cellmap_isoform,
            "cellmap_source": args.cellmap_source,
            "browser_location": page.locator("#browserLocation").inner_text(),
            "evidence_text": page.locator("#geneEvidence").inner_text(),
            "crop_note": (
                "Each crop begins at the six-tab strip and includes the active "
                "view up to a common maximum height. The gene sidebar was hidden "
                "only to use the full publication width; no data were altered."
            ),
            "views": summaries,
        }
        metadata_path = output_dir / args.metadata_filename
        metadata_path.write_text(
            json.dumps(metadata, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        context.close()
        browser.close()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-url", default="http://10.168.3.4:4173")
    parser.add_argument("--gene", default="PTPRC")
    parser.add_argument("--species", default="human")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--chrome",
        default="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    )
    parser.add_argument("--viewport-width", type=int, default=1500)
    parser.add_argument("--viewport-height", type=int, default=1200)
    parser.add_argument("--device-scale-factor", type=float, default=2.0)
    parser.add_argument("--max-crop-height", type=float, default=900.0)
    parser.add_argument(
        "--only-view",
        choices=[view for view, _ in VIEWS],
        default=None,
        help="Capture only one view; default captures all six.",
    )
    parser.add_argument(
        "--cellmap-right-view",
        choices=("gene", "isoform", "ase"),
        default="gene",
        help="Signal shown in the Cell map right panel.",
    )
    parser.add_argument(
        "--cellmap-isoform",
        default="ENST00000367364",
        help="Isoform selected when --cellmap-right-view=isoform.",
    )
    parser.add_argument(
        "--cellmap-source",
        default="",
        help="Optional Cell map source/stage option, matched by visible prefix.",
    )
    parser.add_argument(
        "--metadata-filename",
        default="ptprc_all_views_metadata.json",
        help="Output metadata filename within --output-dir.",
    )
    return parser


if __name__ == "__main__":
    capture(build_parser().parse_args())
