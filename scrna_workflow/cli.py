"""Command-line interface for questions, planning, and running analyses."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import yaml

from .orchestrator import DEPS, execute

PATH_KEYS = ("input_path", "metadata_path", "output_dir", "tf_list", "motif_evidence_path",
             "marker_path", "knowledge_base_path")


def _marker_sets(path):
    """Load a local YAML/JSON marker dictionary; never download references."""
    data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if isinstance(data, dict) and "markers" in data:
        data = data["markers"]
    if not isinstance(data, dict) or not data:
        raise ValueError("Marker file must contain a nonempty mapping of cell type to gene list")
    markers = {}
    for cell_type, genes in data.items():
        if not isinstance(cell_type, str) or not cell_type.strip():
            raise ValueError("Each cell type needs a nonempty name")
        if not isinstance(genes, list) or any(not isinstance(g, str) or not g.strip() for g in genes):
            raise ValueError(f"Markers for {cell_type!r} must be a list of gene symbols")
        cleaned = list(dict.fromkeys(g.strip() for g in genes))
        if len(cleaned) < 2:
            raise ValueError(f"Markers for {cell_type!r} need at least two distinct genes")
        markers[cell_type.strip()] = cleaned
    return markers


def _config(args):
    config = {}
    if args.config:
        config_path = Path(args.config).expanduser().resolve()
        config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        if not isinstance(config, dict):
            raise ValueError("Configuration must be a YAML mapping")
        for key in PATH_KEYS:
            if config.get(key):
                config[key] = str((config_path.parent / Path(config[key]).expanduser()).resolve())
    for arg, key in (("input", "input_path"), ("metadata", "metadata_path"),
                     ("output", "output_dir"), ("markers", "marker_path"),
                     ("knowledge_base", "knowledge_base_path"),
                     ("seed", "seed"), ("workers", "workers")):
        value = getattr(args, arg, None)
        if value is not None:
            config[key] = str(Path(value).expanduser().resolve()) if key in PATH_KEYS else value
    if config.get("marker_path"):
        config["markers"] = _marker_sets(config["marker_path"])
    elif config.get("markers"):
        # Apply the same input checks to inline YAML markers.
        config["markers"] = _validated_inline_markers(config["markers"])
    return config


def _validated_inline_markers(data):
    if not isinstance(data, dict) or not data:
        raise ValueError("markers must map cell types to gene lists")
    result = {}
    for name, genes in data.items():
        if (not isinstance(name, str) or not name.strip() or not isinstance(genes, list)
                or any(not isinstance(g, str) or not g.strip() for g in genes)):
            raise ValueError("Each marker set needs a cell type name and gene-symbol list")
        cleaned = list(dict.fromkeys(g.strip() for g in genes))
        if len(cleaned) < 2:
            raise ValueError(f"Markers for {name!r} need at least two distinct genes")
        result[name.strip()] = cleaned
    return result


def _add_run_options(parser):
    parser.add_argument("--config", help="YAML workflow configuration")
    parser.add_argument("--input", help="AnnData or 10x count input")
    parser.add_argument("--metadata", help="Barcode-indexed CSV/TSV metadata")
    parser.add_argument("--output", help="New output directory")
    parser.add_argument("--markers", help="YAML/JSON cell-type marker sets")
    parser.add_argument("--knowledge-base", help="Local marker knowledge-base YAML")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--workers", type=int)
    parser.add_argument("--resume", action="store_true")


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0].startswith("-") and argv[0] not in {"-h", "--help"}:
        argv.insert(0, "run")  # Preserve the pre-subcommand CLI.
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    _add_run_options(sub.add_parser("run", help="Run or resume analysis"))
    ask = sub.add_parser("ask", help="Analyze a dataset from a plain-language question")
    ask.add_argument("question", nargs="+", help="Question, including a dataset path unless --input is used")
    _add_run_options(ask)
    ask.add_argument("--planner", choices=("ollama", "rules"), default="ollama",
                     help="Local LLM planner (default) or explicit offline rules")
    ask.add_argument("--llm-model", default="qwen2.5:7b", help="Installed local Ollama model")
    ask.add_argument("--annotation-backend", choices=("ollama", "markers"),
                     help="Use local Ollama or direct marker overlap for cell labels")
    chat = sub.add_parser("chat", help="Cluster a dataset once, then ask follow-up questions")
    _add_run_options(chat)
    chat.add_argument("--run-dir", help="Reopen a completed clustering run without rerunning it")
    chat.add_argument("--llm-model", default="qwen2.5:7b", help="Installed local Ollama model")
    chat.add_argument("--annotation-backend", choices=("ollama", "markers"),
                      help="Initial cluster annotation method")
    plan = sub.add_parser("plan", help="Show agent ownership and task dependencies")
    plan.add_argument("--config", help="Optional YAML configuration")
    sub.add_parser("tools", help="List executable scientific tools")
    args = parser.parse_args(argv)
    try:
        if args.command in {"plan", "tools"}:
            from .agents import TASK_AGENTS
            from .tools import TASK_TOOLS
            rows = [
                {"task": name, "agent": TASK_AGENTS[name].name,
                 "tool": f"{TASK_TOOLS[name].__module__}.{TASK_TOOLS[name].__name__}",
                 "depends_on": DEPS[name]}
                for name in DEPS
            ]
            if args.command == "plan" and args.config:
                config = _config(args)
                print(json.dumps({"tasks": rows, "marker_sets": list(config.get("markers", {})),
                                  "input_path": config.get("input_path")}, indent=2))
            else:
                print(json.dumps(rows, indent=2))
            return 0
        if args.command == "ask":
            from .question import prepare_question, summarize_run
            base_config = _config(args)
            if args.resume and not base_config.get("output_dir"):
                raise ValueError("Use --output (or a config with output_dir) to resume a question run.")
            question = " ".join(args.question)
            if args.planner == "ollama":
                from .llm_planner import plan_question
                config = plan_question(question, base_config, model=args.llm_model)
                targets = config["requested_tasks"]
                planner_name = "local_ollama"
            else:
                from .planner import plan_question
                config = prepare_question(question, base_config)
                stage_plan = plan_question(question, "rules")
                targets = stage_plan["targets"]
                config["requested_tasks"] = targets
                config["question_plan"] = stage_plan
                planner_name = "rules"
            if config.get("marker_path") and not config.get("markers"):
                config["markers"] = _marker_sets(config["marker_path"])
            if "clustering" in targets or "report" in targets:
                from .tools.knowledge import built_in_pbmc_knowledge
                if (not config.get("knowledge_base_path") and not config.get("marker_path")
                        and not config.get("markers")):
                    config["knowledge_base_path"] = built_in_pbmc_knowledge(config, question)
                config["annotation_backend"] = (
                    args.annotation_backend or ("ollama" if args.planner == "ollama" else "markers")
                )
                config["annotation_model"] = args.llm_model
            print("Selected stages: " + ", ".join(targets) +
                  f" (planner: {planner_name})", flush=True)
            status = execute(config, resume=args.resume)
            print("\n" + summarize_run(config["output_dir"], question))
            return status
        if args.command == "chat":
            from .interactive import run_chat
            if args.run_dir:
                return run_chat(args.run_dir, model=args.llm_model)
            config = _config(args)
            if not config.get("input_path"):
                raise ValueError("Start chat with --input or a config containing input_path")
            if not config.get("output_dir"):
                stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
                config["output_dir"] = str((Path.cwd() / "results" / f"chat_{stamp}").resolve())
            config["requested_tasks"] = ["clustering"]
            config.setdefault("biological_question", "Interactive clustering and annotation")
            if not config.get("knowledge_base_path") and not config.get("markers"):
                from .tools.knowledge import built_in_pbmc_knowledge
                config["knowledge_base_path"] = built_in_pbmc_knowledge(config, "")
            config["annotation_backend"] = args.annotation_backend or "ollama"
            config["annotation_model"] = args.llm_model
            status = execute(config, resume=args.resume)
            return status if status else run_chat(config["output_dir"], model=args.llm_model)
        return execute(_config(args), resume=args.resume)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
