"""Compatibility imports for core analysis tools.

Scientific implementations live in :mod:`scrna_workflow.tools`.
"""
from .tools.inspection import inspect_data
from .tools.qc import qc
from .tools.preprocess import preprocess
from .tools.neighbors import build_neighbors
from .tools.clustering import clustering
from .tools._core_common import _countlike

representation = preprocess
graph = build_neighbors
