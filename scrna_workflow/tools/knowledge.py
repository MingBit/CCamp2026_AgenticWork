"""Versioned local marker knowledge bases for grounded annotation."""
from pathlib import Path
import yaml


def load_knowledge_base(path):
    path = Path(path).expanduser().resolve()
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or not isinstance(data.get("markers"), dict):
        raise ValueError("Knowledge base must contain a marker mapping")
    markers = {}
    for label, genes in data["markers"].items():
        if not isinstance(label, str) or not label.strip() or not isinstance(genes, list):
            raise ValueError("Knowledge base labels need nonempty gene lists")
        cleaned = list(dict.fromkeys(str(gene).strip().upper() for gene in genes))
        if len(cleaned) < 2 or any(not gene for gene in cleaned):
            raise ValueError(f"Knowledge base label {label!r} needs two distinct marker genes")
        markers[label] = cleaned
    return markers, {"path": str(path), "name": data.get("name"),
                     "organism": data.get("organism"), "tissue": data.get("tissue"),
                     "source": data.get("source")}


def built_in_pbmc_knowledge(config, question):
    """Offer the bundled human PBMC panel only when species and tissue are explicit."""
    organism = str(config.get("organism") or "").lower()
    tissue = str(config.get("tissue") or "").lower()
    text = question.lower()
    human = organism in {"human", "homo sapiens"} or "human" in text
    pbmc = any(term in tissue or term in text for term in
               ("pbmc", "peripheral blood mononuclear"))
    if human and pbmc:
        return str((Path(__file__).resolve().parents[1] / "knowledge" /
                    "human_pbmc_markers.yaml").resolve())
    return None
