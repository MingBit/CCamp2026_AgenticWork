"""Executable scientific tools used by the workflow DAG.

Each public stage accepts the runner context and returns a structured result.
Large AnnData objects travel through file paths, not through result dictionaries.
"""
from .inspection import inspect_data
from .qc import qc
from .preprocess import preprocess
from .neighbors import build_neighbors
from .clustering import clustering
from .regulon import regulon
from .discovery import discovery
from .validation import validation
from .report import report

TASK_TOOLS = {
    "inspection": inspect_data,
    "qc": qc,
    "representation": preprocess,
    "graph": build_neighbors,
    "clustering": clustering,
    "regulon": regulon,
    "discovery": discovery,
    "validation": validation,
    "report": report,
}
