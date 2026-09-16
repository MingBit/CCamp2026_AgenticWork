"""Deterministic specialist ownership for the workflow DAG."""
from .cell_state import AGENT as CELL_STATE
from .network import AGENT as NETWORK
from .regulon import AGENT as REGULON
from .discovery import AGENT as DISCOVERY
from .evaluator import AGENT as EVALUATOR
from .reporter import AGENT as REPORTER

SPECIALISTS = (CELL_STATE, NETWORK, REGULON, DISCOVERY, EVALUATOR, REPORTER)
TASK_AGENTS = {task: agent for agent in SPECIALISTS for task in agent.tasks}
if len(TASK_AGENTS) != sum(len(agent.tasks) for agent in SPECIALISTS):
    raise RuntimeError("A workflow task has more than one specialist owner")
