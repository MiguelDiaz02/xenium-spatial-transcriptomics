"""Rich-based logger for pipeline scripts."""

from __future__ import annotations

import logging
import sys

from rich.console import Console
from rich.logging import RichHandler


def get_logger(name: str, log_file: str | None = None) -> logging.Logger:
    """
    Return a logger that writes to stderr (rich) and optionally to a file.

    Parameters
    ----------
    name
        Logger name (use __name__ in each script).
    log_file
        If provided, also write plain-text logs to this path.
    """
    handlers: list[logging.Handler] = [
        RichHandler(console=Console(stderr=True), show_path=False)
    ]

    if log_file:
        from pathlib import Path
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        datefmt="[%X]",
        handlers=handlers,
    )
    return logging.getLogger(name)
