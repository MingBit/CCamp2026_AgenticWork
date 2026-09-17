"""Shared file-backed interface for core analysis tools."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy import sparse

def _result(outputs=None, metrics=None, warnings=None, status="completed", inputs=None):
    return {
        "status": status,
        "input_references": inputs or [],
        "outputs": outputs or {},
        "metrics": metrics or {},
        "warnings": warnings or [],
        "recommended_next_actions": [],
    }


def _dir(ctx, name):
    directory = Path(ctx["root"]) / name
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def _previous(ctx, task, key="dataset"):
    result = ctx["artifacts"][task]
    return result.get("outputs", result)[key]


def _load_step(ctx, parent, name):
    import anndata as ad

    source = _previous(ctx, parent)
    return ad.read_h5ad(source), _dir(ctx, name), source


def _json(path, value):
    path.write_text(json.dumps(value, indent=2, default=str))
    return str(path)


def _save(adata, path):
    adata.write_h5ad(path)
    return str(path)


def _digest(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _countlike(matrix):
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    return bool(
        np.isfinite(values).all()
        and (values >= 0).all()
        and np.allclose(values, np.round(values), atol=1e-6, rtol=0)
    )
