"""Logging configuration."""

from __future__ import annotations

import logging


def get_logger(name: str) -> logging.Logger:
    """Get a logger for the given module name.

    All eukan loggers are children of the 'eukan' root logger,
    which is configured once by setup_logging().
    """
    return logging.getLogger(name)


def setup_logging(verbosity: int = 0) -> None:
    """Configure the eukan logging hierarchy.

    Args:
        verbosity: 0=INFO, 1=DEBUG, -1=WARNING (quiet).
    """
    level = {-1: logging.WARNING, 0: logging.INFO, 1: logging.DEBUG}.get(
        verbosity, logging.INFO
    )

    root = logging.getLogger("eukan")
    if root.handlers:
        return  # already configured

    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s", datefmt="%H:%M:%S")
    )
    root.addHandler(handler)
    root.setLevel(level)
