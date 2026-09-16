"""Owns inspection, QC, expression preprocessing, and cell identity."""
from .base import Specialist

AGENT = Specialist("cell_state", ("inspection", "qc", "representation", "clustering"))
