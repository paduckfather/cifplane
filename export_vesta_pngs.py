#!/usr/bin/env python3
"""Batch-export VESTA PNGs from *_vesta_view.cif files via macOS GUI automation."""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path


VESTA_APP = "/Applications/VESTA.app"


def run_command(args: list[str], check: bool = True) -> subprocess.CompletedProcess[str]:
    """Run one subprocess and optionally raise with combined stderr/stdout on failure."""
    completed = subprocess.run(args, capture_output=True, text=True)
    if check and completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip() or f"Command failed: {args!r}"
        raise RuntimeError(message)
    return completed


def run_osascript(lines: list[str], check: bool = True) -> str:
    """Run AppleScript lines through osascript."""
    cmd = ["osascript"]
    for line in lines:
        cmd.extend(["-e", line])
    completed = run_command(cmd, check=check)
    return completed.stdout.strip()


def get_window_names() -> list[str]:
    """Return current VESTA window names."""
    output = run_osascript(
        [
            'tell application "System Events"',
            'tell process "VESTA"',
            "set windowNames to name of every window",
            "end tell",
            "end tell",
            "set AppleScript's text item delimiters to linefeed",
            "return windowNames as text",
        ],
        check=False,
    )
    if not output:
        return []
    return [line.strip() for line in output.splitlines() if line.strip()]


def vesta_is_running() -> bool:
    """Return whether a VESTA process is currently running."""
    output = run_osascript(
        [
            'tell application "System Events" to return exists process "VESTA"',
        ],
        check=False,
    )
    return output.lower() == "true"


def wait_for(predicate, timeout: float, poll_interval: float, description: str) -> None:
    """Wait until predicate() becomes truthy or raise TimeoutError."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        if predicate():
            return
        time.sleep(poll_interval)
    raise TimeoutError(f"Timed out while waiting for {description}.")


def quit_vesta(grace_timeout: float = 1.2) -> None:
    """Quit VESTA quickly, falling back to a forced kill when needed."""
    if not vesta_is_running():
        return
    run_osascript(['tell application "VESTA" to quit'], check=False)
    try:
        wait_for(lambda: not vesta_is_running(), timeout=grace_timeout, poll_interval=0.2, description="VESTA to quit")
        return
    except TimeoutError:
        run_command(["pkill", "-x", "VESTA"], check=False)
        wait_for(lambda: not vesta_is_running(), timeout=2.0, poll_interval=0.2, description="VESTA to force-quit")


def open_vesta_file(cif_path: Path) -> None:
    """Open one CIF file in a fresh VESTA instance."""
    run_command(["open", "-a", VESTA_APP, "--args", "-open", str(cif_path)], check=True)


def click_export_raster_image() -> None:
    """Trigger File -> Export Raster Image..."""
    run_osascript(
        [
            'tell application "VESTA" to activate',
            'tell application "System Events" to tell process "VESTA" to click menu item "Export Raster Image..." of menu "File" of menu bar item "File" of menu bar 1',
        ],
        check=True,
    )


def get_save_as_filename() -> str:
    """Read the default filename from the Save dialog."""
    return run_osascript(
        [
            'tell application "System Events" to tell process "VESTA" to get value of text field "Save As:" of splitter group 1 of window "Save"',
        ],
        check=True,
    )


def click_save_button() -> None:
    """Click the Save button in the Save dialog."""
    run_osascript(
        [
            'tell application "VESTA" to activate',
            'tell application "System Events" to tell process "VESTA" to click button "Save" of splitter group 1 of window "Save"',
        ],
        check=True,
    )


def set_transparent_background_if_needed() -> None:
    """Enable transparent background in the Export image dialog when present."""
    if "Export image" not in get_window_names():
        return
    current_value = run_osascript(
        [
            'tell application "System Events" to tell process "VESTA" to get value of checkbox "Let the background transparent" of window 1',
        ],
        check=False,
    ).strip()
    if current_value == "1":
        return
    run_osascript(
        [
            'tell application "VESTA" to activate',
            'tell application "System Events" to tell process "VESTA" to click checkbox "Let the background transparent" of window 1',
        ],
        check=False,
    )
    time.sleep(0.1)


def dismiss_export_dialog_if_needed() -> None:
    """Dismiss the transient Export image dialog after Save."""
    names = get_window_names()
    if "Export image" not in names:
        return
    time.sleep(0.15)
    set_transparent_background_if_needed()
    for key_code in (36, 76):
        run_osascript(
            [
                'tell application "VESTA" to activate',
                f'tell application "System Events" to key code {key_code}',
            ],
            check=False,
        )
        time.sleep(0.15)
        if "Export image" not in get_window_names():
            return

    run_osascript(
        [
            'tell application "VESTA" to activate',
            'tell application "System Events" to tell process "VESTA" to click button "OK" of window 1',
        ],
        check=False,
    )
    time.sleep(0.15)


def export_one(cif_path: Path) -> tuple[Path, str]:
    """Export one VESTA PNG next to one *_vesta_view.cif file."""
    png_path = cif_path.with_suffix(".png")
    quit_vesta()
    if png_path.exists():
        return png_path, "skipped-existing"

    time.sleep(0.2)
    open_vesta_file(cif_path)

    expected_window = f"{cif_path.name} - VESTA"
    wait_for(
        lambda: expected_window in get_window_names() or "VESTA" in get_window_names(),
        timeout=30.0,
        poll_interval=0.2,
        description=expected_window,
    )

    click_export_raster_image()
    wait_for(lambda: "Save" in get_window_names(), timeout=10.0, poll_interval=0.2, description="Save dialog")

    default_filename = get_save_as_filename()
    if default_filename != png_path.name:
        raise RuntimeError(f"Unexpected default export filename: {default_filename!r} != {png_path.name!r}")

    click_save_button()

    def export_completed() -> bool:
        dismiss_export_dialog_if_needed()
        names = get_window_names()
        dialog_closed = "Export image" not in names and "Save" not in names
        return png_path.exists() and png_path.stat().st_size > 0 and dialog_closed

    wait_for(export_completed, timeout=15.0, poll_interval=0.2, description=str(png_path))
    quit_vesta()
    return png_path, "exported"


def collect_targets(paths: list[str]) -> list[Path]:
    """Expand files/directories into a sorted list of *_vesta_view.cif targets."""
    targets: list[Path] = []
    for raw_path in paths:
        path = Path(raw_path).expanduser().resolve()
        if path.is_dir():
            targets.extend(sorted(path.rglob("*_vesta_view.cif")))
        elif path.is_file():
            targets.append(path)
        else:
            raise FileNotFoundError(f"Path not found: {path}")
    targets = sorted(dict.fromkeys(targets))
    if not targets:
        raise FileNotFoundError("No *_vesta_view.cif targets were found.")
    return targets


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description="Export VESTA-native PNGs for *_vesta_view.cif files.")
    parser.add_argument("paths", nargs="+", help="One or more *_vesta_view.cif files or directories containing them.")
    return parser.parse_args()


def main() -> int:
    """CLI entry point."""
    args = parse_args()
    targets = collect_targets(args.paths)

    exported = 0
    skipped = 0
    for index, cif_path in enumerate(targets, start=1):
        png_path, status = export_one(cif_path)
        print(f"[{index}/{len(targets)}] {status}: {png_path}")
        if status == "exported":
            exported += 1
        else:
            skipped += 1

    print(f"exported={exported} skipped={skipped} total={len(targets)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
