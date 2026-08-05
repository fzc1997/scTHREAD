#!/usr/bin/env python3
"""Capture publication-ready PTPRC portal crops from a live scTHREAD instance."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from playwright.sync_api import Page, sync_playwright


def _wait_for_ptprc(page: Page, gene: str) -> None:
    page.wait_for_selector("#browserSearch", state="visible")
    try:
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
    except Exception:
        page.fill("#browserSearch", gene)
        page.click("#browserGo")
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


def _clip_between(page: Page, top_selector: str, bottom_selector: str) -> dict[str, float]:
    workspace = page.locator("#transcriptBrowser .browser-workspace").bounding_box()
    top = page.locator(top_selector).bounding_box()
    bottom = page.locator(bottom_selector).bounding_box()
    if not workspace or not top or not bottom:
        raise RuntimeError("Could not resolve the requested publication crop.")
    pad = 2.0
    return {
        "x": max(0.0, workspace["x"] - pad),
        "y": max(0.0, top["y"] - pad),
        "width": workspace["width"] + 2 * pad,
        "height": bottom["y"] + bottom["height"] - top["y"] + 2 * pad,
    }


def capture(args: argparse.Namespace) -> None:
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    base_url = args.base_url.rstrip("/")
    page_url = f"{base_url}/browse?species=human&query={args.gene}"

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=args.chrome or None,
            headless=True,
            args=["--disable-gpu", "--no-sandbox", "--force-device-scale-factor=1"],
        )
        context = browser.new_context(
            viewport={"width": args.viewport_width, "height": args.viewport_height},
            device_scale_factor=args.device_scale_factor,
            color_scheme="light",
            reduced_motion="reduce",
        )
        page = context.new_page()
        page.goto(page_url, wait_until="domcontentloaded", timeout=60_000)
        _wait_for_ptprc(page, args.gene)

        # The sidebar and detailed data views are outside the evidence-card crop.
        # Hiding them changes no data; it only removes unused browser chrome.
        page.add_style_tag(
            content="""
                .browser-shell { grid-template-columns: minmax(0, 1fr) !important; }
                .browser-sidebar { display: none !important; }
                .browser-workspace { min-width: 0 !important; }
                * { animation: none !important; transition: none !important; }
            """
        )

        card_path = output_dir / "ptprc_gene_card_live_v3.png"
        evidence_path = output_dir / "ptprc_evidence_strip_live_v3.png"
        isoform_path = output_dir / "ptprc_isoforms_live_v3.png"

        page.screenshot(
            path=str(card_path),
            clip=_clip_between(page, ".browser-toolbar", ".browser-view-tabs"),
            animations="disabled",
        )
        page.locator("#geneEvidence").screenshot(
            path=str(evidence_path),
            animations="disabled",
        )

        page.click('button[data-browser-view="isoforms"]')
        page.wait_for_selector("#isoformsOutput .track-row", state="visible", timeout=60_000)
        page.locator("#isoformsView").screenshot(
            path=str(isoform_path),
            animations="disabled",
        )

        metadata = {
            "captured_at_utc": datetime.now(timezone.utc).isoformat(),
            "source_url": page.url,
            "gene": args.gene,
            "viewport_css_px": [args.viewport_width, args.viewport_height],
            "device_scale_factor": args.device_scale_factor,
            "browser_location": page.locator("#browserLocation").inner_text(),
            "evidence_text": page.locator("#geneEvidence").inner_text(),
            "crop_note": (
                "The gene-list sidebar and detailed data views were hidden only to crop "
                "the live toolbar, location, evidence strip and tabs without altering data."
            ),
            "files": [card_path.name, evidence_path.name, isoform_path.name],
        }
        (output_dir / "ptprc_live_capture_v3.json").write_text(
            json.dumps(metadata, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        context.close()
        browser.close()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-url", default="https://scthread.ai4sc.ac.cn")
    parser.add_argument("--gene", default="PTPRC")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--chrome",
        default=None,
    )
    parser.add_argument("--viewport-width", type=int, default=1500)
    parser.add_argument("--viewport-height", type=int, default=1200)
    parser.add_argument("--device-scale-factor", type=float, default=2.0)
    return parser


if __name__ == "__main__":
    capture(build_parser().parse_args())
