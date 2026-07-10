"""Reserved browser automation interface for future external tool adapters."""

from __future__ import annotations

from pathlib import Path


class BrowserAutomationAdapter:
    """Non-operational browser automation adapter.

    Browser automation is intentionally disabled by default in AbDev-Lite v0.9.
    This class does not install browser drivers, open websites, submit sequences,
    bypass CAPTCHA, or save credentials.
    """

    def __init__(self, enabled: bool = False):
        self.enabled = bool(enabled)

    def run(self, tool_id: str, input_file: Path):
        if not self.enabled:
            return {
                "tool_id": tool_id,
                "input_file": str(input_file),
                "run_status": "BROWSER_AUTOMATION_DISABLED",
                "message": "Browser automation is disabled by default in AbDev-Lite v0.9.",
            }
        return {
            "tool_id": tool_id,
            "input_file": str(input_file),
            "run_status": "BROWSER_AUTOMATION_NOT_IMPLEMENTED",
            "message": "No real browser automation is implemented in v0.9.",
        }
