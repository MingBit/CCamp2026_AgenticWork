"""Shared result and artifact helpers for downstream tools."""
from pathlib import Path
import json

def _get(ctx, key, default=None):
    return ctx.get(key, default) if isinstance(ctx, dict) else getattr(ctx, key, default)


def _base(ctx, name):
    directory = Path(_get(ctx, "root")) / name
    directory.mkdir(parents=True, exist_ok=True)
    return directory, _get(ctx, "config", {}), _get(ctx, "artifacts", {})


def _result(status="completed", inputs=None, outputs=None, metrics=None, warnings=None, actions=None):
    return dict(status=status, input_references=inputs or [], outputs=outputs or {},
                metrics=metrics or {}, warnings=warnings or [], recommended_next_actions=actions or [])


def _dataset(artifacts):
    for task in ("clustering", "qc", "inspection"):
        outputs = artifacts.get(task, {}).get("outputs", {})
        for key in ("dataset", "annotated_dataset", "filtered_dataset", "adata"):
            candidate = outputs.get(key)
            if candidate and Path(candidate).suffix == ".h5ad" and Path(candidate).is_file():
                return str(candidate)
    return None


def _json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False, default=str))
    return str(path)
