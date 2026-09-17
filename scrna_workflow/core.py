"""Compatibility imports for core analysis tools.

Scientific implementations live in :mod:`scrna_workflow.tools`.
"""
import json
from pathlib import Path

from .tools.inspection import inspect_data
from .tools.qc import qc
from .tools.preprocess import preprocess
from .tools.neighbors import build_neighbors
from .tools.clustering import clustering
from .tools._core_common import _countlike

representation = preprocess
graph = build_neighbors


def _dir(ctx, name):
    """Ensure directory exists within context root and return Path."""
    p = Path(ctx["root"]) / name
    p.mkdir(parents=True, exist_ok=True)
    return p


def _json(path, data):
    """Write structured JSON payload to target path."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2, default=str) + "\n")
    return str(p)


def _previous(ctx, task):
    """Extract dataset output path from preceding artifact context."""
    item = ctx.get("artifacts", {}).get(task, {})
    outputs = item.get("outputs", item)
    return outputs.get("dataset")


def _result(outputs=None, metrics=None, warnings=None, status="completed", inputs=None):
    """Construct standardized task result dictionary."""
    return {
        "status": status,
        "input_references": inputs or [],
        "outputs": outputs or {},
        "metrics": metrics or {},
        "warnings": warnings or [],
        "recommended_next_actions": []
    }