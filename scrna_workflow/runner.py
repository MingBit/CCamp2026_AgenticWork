"""Compatibility imports for the orchestrator and legacy CLI."""
from .orchestrator import DEPS, REVIEW_TASKS, digest, write_json, valid_checkpoint, execute


def main():
    from .cli import main as cli_main
    return cli_main()
