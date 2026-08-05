#!/usr/bin/env python3
"""Run phone/tablet acceptance tasks against a deployed scTHREAD portal."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from playwright.sync_api import sync_playwright


VIEWPORTS = {
    "phone": {"width": 390, "height": 844},
    "tablet": {"width": 820, "height": 1180},
}


def require_no_document_overflow(page, label: str) -> None:
    overflow = page.evaluate(
        """() => ({
          documentWidth: document.documentElement.scrollWidth,
          viewportWidth: document.documentElement.clientWidth,
          offenders: [...document.querySelectorAll('body *')]
            .map(element => {
              const rect = element.getBoundingClientRect();
              return {
                tag: element.tagName,
                id: element.id,
                className: String(element.className || '').slice(0, 100),
                left: Math.round(rect.left),
                right: Math.round(rect.right),
                width: Math.round(rect.width)
              };
            })
            .filter(item => item.right > document.documentElement.clientWidth + 2 || item.left < -2)
            .slice(0, 20)
        })"""
    )
    if overflow["documentWidth"] > overflow["viewportWidth"] + 2:
        raise AssertionError(f"{label}: document overflow {overflow}")


def timed_goto(page, url: str) -> float:
    started = time.monotonic()
    page.goto(url, wait_until="networkidle", timeout=60_000)
    return round(time.monotonic() - started, 3)


def run(base_url: str, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    report: dict[str, object] = {
        "base_url": base_url,
        "viewports": {},
        "status": "pass",
    }
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            headless=True,
            executable_path="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
        )
        try:
            for name, viewport in VIEWPORTS.items():
                context = browser.new_context(
                    viewport=viewport,
                    device_scale_factor=1,
                    locale="en-US",
                )
                page = context.new_page()
                console_errors: list[str] = []
                page.on(
                    "console",
                    lambda message: (
                        console_errors.append(message.text)
                        if message.type == "error"
                        else None
                    ),
                )
                timings: dict[str, float] = {}

                timings["home_seconds"] = timed_goto(page, f"{base_url}/")
                page.locator("#releaseContextBanner").wait_for(state="visible")
                page.get_by_text("Manuscript snapshot", exact=False).wait_for()
                page.locator(".navbar-toggle").click()
                page.locator("#navbarContent").wait_for(state="visible")
                require_no_document_overflow(page, f"{name} home")
                page.screenshot(path=output / f"{name}_01_home.png")

                timings["search_seconds"] = timed_goto(
                    page, f"{base_url}/search?q=PTPRC"
                )
                page.locator("#globalSearchResults .result-list").wait_for(
                    state="visible", timeout=30_000
                )
                page.get_by_text("PTPRC", exact=True).first.wait_for()
                require_no_document_overflow(page, f"{name} search")
                page.screenshot(path=output / f"{name}_02_search_ptprc.png")

                timings["browse_seconds"] = timed_goto(
                    page, f"{base_url}/browse?query=PTPRC&species=human"
                )
                page.get_by_text("Database-wide evidence", exact=True).wait_for(
                    timeout=30_000
                )
                page.locator("#browserTableBody tr").first.wait_for(
                    state="visible", timeout=30_000
                )
                require_no_document_overflow(page, f"{name} browse")
                page.screenshot(path=output / f"{name}_03_browse_ptprc.png")

                timings["download_seconds"] = timed_goto(
                    page, f"{base_url}/download"
                )
                page.locator("#sampleManifestBody tr").first.wait_for(
                    state="visible", timeout=30_000
                )
                page.get_by_text("Frozen 469-run release TSV", exact=True).wait_for()
                require_no_document_overflow(page, f"{name} download")
                page.screenshot(path=output / f"{name}_04_download.png")

                report["viewports"][name] = {
                    "viewport": viewport,
                    "tasks": [
                        "home_and_mobile_navigation",
                        "PTPRC_search",
                        "PTPRC_gene_browser",
                        "release_download",
                    ],
                    "timings": timings,
                    "console_errors": console_errors,
                }
                context.close()
        finally:
            browser.close()

    (output / "responsive_acceptance.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base_url")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(args.base_url.rstrip("/"), args.output)


if __name__ == "__main__":
    main()
