"""Small, explicit boundary between a specialist and its scientific tools."""
from dataclasses import dataclass

from ..tools import TASK_TOOLS


@dataclass(frozen=True)
class Specialist:
    name: str
    tasks: tuple[str, ...]

    def execute(self, task: str, context: dict) -> dict:
        if task not in self.tasks:
            raise ValueError(f"{self.name} does not own task {task!r}")
        return TASK_TOOLS[task](context)
