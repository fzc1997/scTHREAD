#!/usr/bin/env python3
"""Capture a current, internally consistent scTHREAD portal walkthrough."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from playwright.sync_api import sync_playwright


PAGES = (
    ("01_home_current.png", "/", "scTHREAD"),
    ("02_search_PTPRC_current.png", "/search?q=PTPRC&species=human", "PTPRC"),
    ("05_download_current.png", "/download", "Frozen release TSV"),
    ("06_analyze_current.png", "/analyze", "Compare cell types."),
    ("07_about_current.png", "/about", "About scTHREAD"),
    ("08_docs_current.png", "/docs", "scTHREAD API"),
)


def capture(args: argparse.Namespace) -> None:
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=args.chrome or None,
            headless=True,
            args=["--disable-gpu", "--no-sandbox", "--force-device-scale-factor=1"],
        )
        context = browser.new_context(
            viewport={"width": 1440, "height": 900},
            device_scale_factor=2.0,
            color_scheme="light",
            reduced_motion="reduce",
        )
        for filename, route, expected in PAGES:
            page = context.new_page()
            url = f"{args.base_url.rstrip('/')}{route}"
            page.goto(url, wait_until="domcontentloaded", timeout=60_000)
            page.wait_for_function(
                "expected => document.body.innerText.includes(expected)",
                arg=expected,
                timeout=60_000,
            )
            page.wait_for_timeout(1_000)
            output = output_dir / filename
            page.screenshot(path=str(output), full_page=False, animations="disabled")
            body_text = page.locator("body").inner_text()
            records.append(
                {
                    "file": filename,
                    "url": page.url,
                    "expected_text": expected,
                    "expected_text_present": expected in body_text,
                    "body_text_excerpt": body_text[:2_000],
                    "png_bytes": output.stat().st_size,
                }
            )
            page.close()
        context.close()
        browser.close()
    metadata = {
        "captured_at_utc": datetime.now(timezone.utc).isoformat(),
        "base_url": args.base_url,
        "viewport_css_px": [1440, 900],
        "device_scale_factor": 2.0,
        "pages": records,
    }
    (output_dir / "portal_walkthrough_current_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--base-url", default="https://scthread.ai4sc.ac.cn")
    result.add_argument("--output-dir", required=True)
    result.add_argument(
        "--chrome",
        default=None,
        help="browser executable; omit to use Playwright's bundled Chromium",
    )
    return result


if __name__ == "__main__":
    capture(parser().parse_args())
