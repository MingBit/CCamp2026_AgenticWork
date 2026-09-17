"""Owns independent integrity checks and critique."""
from .base import Specialist

AGENT = Specialist("evaluator", ("validation",))
