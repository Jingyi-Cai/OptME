# coding:utf-8
"""
Author: Jingyi Cai (2024-2026)
Function: Run the C. glutamicum iCW773 single-task pipeline. Asks for the solver, optional solver parameters, and optional LLM settings, then calls Error.py, summary.py, outputjson_status.py, and sum.py.
Input: Environment variables OPTME_SOLVER_PATH and OPTME_SOLVER. Optional OPTME_SOLVER_OPTIONS holds extra solver parameters such as a time limit. Thread count is chosen from the CPU count. Optional OPTME_LLM_MODEL, OPTME_LLM_BASE, and OPTME_LLM_KEY. Optional OPTME_METHODS limits the run to a comma-separated subset, for example iBridge,llm. Task file INPUT/task/iCW773R_task.json, plus INPUT/model and INPUT/map.
Output: Method results under output/iCW773R_task/ and a runtime log at output/iCW773R_task.log.
"""
import os
import sys

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_EXAMPLE_ROOT = os.path.dirname(_SCRIPT_DIR)
_LAUNCH_DIR = os.getcwd()
os.chdir(_EXAMPLE_ROOT)
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

from task_config import (  # noqa: E402
    OPTME_PYTHON,
    OPTME_DIR,
    OUTPUT_DIR,
    MODEL_DIR,
    TASK_DIR,
    MAP_DIR,
)

import json
import subprocess
import time

from solver_budget import use_sequential_threads

SEP = "=" * 60
TASK_FILE = "iCW773R_task.json"
JOB_ID = "iCW773R_task"
_SOLVER_WALK_DEPTH = 6


def prompt(message, default=None):
    """Print a prompt and return stripped input; return default on blank entry."""
    if default:
        suffix = f" [{default}]: "
    else:
        suffix = ": "
    value = input(message + suffix).strip()
    return value if value else default


def confirm(message):
    """Ask a yes/no question; return True for yes."""
    while True:
        ans = input(f"{message} (y/n): ").strip().lower()
        if ans in ("y", "yes"):
            return True
        if ans in ("n", "no"):
            return False
        print("  Please enter y or n.")


def _format_duration(seconds):
    total = int(max(0.0, float(seconds)))
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f"{hours}:{minutes:02d}:{secs:02d}"
    return f"{minutes:02d}:{secs:02d}"


def _append_runtime_summary(log_file, step_times, pipeline_elapsed):
    """Append method and pipeline-step timings to the end of the log."""
    method_rows = []
    if os.path.isfile(log_file):
        with open(log_file, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith("METHOD_RUNTIME\t"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                try:
                    elapsed = float(parts[2])
                except ValueError:
                    continue
                method_rows.append((parts[1], elapsed, parts[3]))

    lines = ["", "=== Runtime summary ===", "Methods:"]
    if method_rows:
        method_total = 0.0
        for name, elapsed, status in method_rows:
            method_total += elapsed
            lines.append(
                f"  {name:<36} {_format_duration(elapsed):>10}  {elapsed:10.1f} s  {status}"
            )
        lines.append(
            f"  {'methods total':<36} {_format_duration(method_total):>10}  {method_total:10.1f} s"
        )
    else:
        lines.append("  (no method timings recorded)")
    lines.append("Pipeline steps:")
    for label, elapsed, status in step_times:
        lines.append(
            f"  {label:<36} {_format_duration(elapsed):>10}  {elapsed:10.1f} s  {status}"
        )
    lines.append(
        f"  {'pipeline total':<36} {_format_duration(pipeline_elapsed):>10}  {pipeline_elapsed:10.1f} s"
    )
    lines.append("")
    text = "\n".join(lines) + "\n"
    with open(log_file, "a", encoding="utf-8") as handle:
        handle.write(text)
    print(text, end="")


def run_command(command, log_file):
    """Run a shell command, stream output to log_file, return (ok, elapsed_seconds)."""
    print(f"  Executing: {' '.join(command)}")
    started = time.perf_counter()
    with open(log_file, "a") as log:
        log.write(f"--- Running command: {' '.join(command)} ---\n")
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=os.environ.copy(),
        )
        stdout, _ = process.communicate()
        elapsed = time.perf_counter() - started
        log.write(stdout or "")
        log.write(f"--- Exit Code: {process.returncode} ---\n")
        log.write(f"--- Elapsed: {elapsed:.3f} s ---\n\n")
        if stdout:
            print(stdout, end="" if stdout.endswith("\n") else "\n")
        if process.returncode != 0:
            print(
                f"  ERROR: Command failed with exit code {process.returncode}. See {log_file} for details."
            )
            print(f"  Elapsed: {elapsed:.1f} s")
            return False, elapsed
    print(f"  SUCCESS: Command finished successfully. Elapsed: {elapsed:.1f} s")
    return True, elapsed


def _split_paths(raw):
    """Keep the path the user typed. Resolve a relative solver path from the launch directory."""
    parts = []
    for chunk in raw.replace(";", os.pathsep).split(os.pathsep):
        chunk = os.path.expanduser(chunk.strip().strip('"').strip("'"))
        if not chunk:
            continue
        if not os.path.isabs(chunk):
            chunk = os.path.normpath(os.path.join(_LAUNCH_DIR, chunk))
        parts.append(chunk)
    return parts


def _walk_named(root, names):
    """Yield files under root whose basename is in names, with a depth cap."""
    if os.path.isfile(root):
        if os.path.basename(root) in names:
            yield root
        return
    root = os.path.normpath(root)
    for dirpath, dirnames, filenames in os.walk(root):
        depth = dirpath[len(root):].count(os.sep)
        if depth > _SOLVER_WALK_DEPTH:
            dirnames[:] = []
            continue
        dirnames[:] = [d for d in dirnames if d not in ("__pycache__", ".git", "node_modules")]
        for filename in filenames:
            if filename in names:
                yield os.path.join(dirpath, filename)


def configure_solver_path(solver_path):
    """Put CPLEX and/or Gurobi from solver_path onto PATH and PYTHONPATH.

    Returns (ok, message).
    """
    entries = _split_paths(solver_path)
    if not entries:
        return False, "Solver path is empty."
    missing = [p for p in entries if not os.path.exists(p)]
    if missing:
        return False, "Path does not exist: " + ", ".join(missing)

    bin_dirs = []
    python_dirs = []
    found = []
    license_file = None

    for entry in entries:
        gurobi_hit = False
        for gurobi_sh in _walk_named(entry, {"gurobi.sh"}):
            bindir = os.path.dirname(gurobi_sh)
            bin_dirs.append(bindir)
            os.environ["GUROBI_HOME"] = os.path.dirname(bindir)
            found.append("Gurobi " + gurobi_sh)
            gurobi_hit = True
            break
        if not gurobi_hit:
            for gurobi_cl in _walk_named(entry, {"gurobi_cl"}):
                bindir = os.path.dirname(gurobi_cl)
                bin_dirs.append(bindir)
                os.environ["GUROBI_HOME"] = os.path.dirname(bindir)
                found.append("Gurobi " + gurobi_cl)
                break

        for cplex_bin in _walk_named(entry, {"cplex"}):
            if os.path.isfile(cplex_bin) and os.access(cplex_bin, os.X_OK):
                bin_dirs.append(os.path.dirname(cplex_bin))
                found.append("CPLEX " + cplex_bin)
                break

        cplex_python = []
        for init_py in _walk_named(entry, {"__init__.py"}):
            package_dir = os.path.dirname(init_py)
            if os.path.basename(package_dir) == "cplex":
                cplex_python.append(os.path.dirname(package_dir))
        if cplex_python:
            want = f"{sys.version_info.major}.{sys.version_info.minor}"
            matched = [p for p in cplex_python if os.sep + want + os.sep in p]
            chosen = matched[0] if matched else cplex_python[0]
            python_dirs.append(chosen)
            found.append("CPLEX Python " + os.path.join(chosen, "cplex"))
            if not matched:
                found.append(
                    f"warning: no CPLEX bindings for Python {want}; using {chosen}"
                )

        for lic in _walk_named(entry, {"gurobi.lic"}):
            license_file = lic
            break

    if not found:
        return False, (
            "No CPLEX or Gurobi installation was found under that path. "
            "Point to a CPLEX Studio directory and/or a Gurobi linux64 directory."
        )

    if bin_dirs:
        os.environ["PATH"] = os.pathsep.join(bin_dirs + [os.environ.get("PATH", "")])
    if python_dirs:
        os.environ["PYTHONPATH"] = os.pathsep.join(
            python_dirs + [os.environ.get("PYTHONPATH", "")]
        )
        for python_dir in python_dirs:
            if python_dir not in sys.path:
                sys.path.insert(0, python_dir)
    if license_file and not os.environ.get("GRB_LICENSE_FILE"):
        os.environ["GRB_LICENSE_FILE"] = license_file

    os.environ["OPTME_SOLVER_PATH"] = solver_path
    has_cplex = any(item.startswith("CPLEX") for item in found)
    has_gurobi = any(item.startswith("Gurobi") for item in found)
    available = []
    if has_cplex:
        available.append("cplex")
    if has_gurobi:
        available.append("gurobi")
    os.environ["OPTME_SOLVERS_FOUND"] = ",".join(available)
    message = "Using " + "; ".join(found)
    if available:
        message += ". Available solvers: " + ", ".join(available)
    return True, message


def _coerce_option(value):
    text = value.strip()
    lowered = text.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    try:
        if any(ch in text for ch in ".eE"):
            return float(text)
        return int(text)
    except ValueError:
        return text


def parse_solver_options(raw):
    """Parse solver parameters from a JSON object or comma-separated key=value pairs."""
    text = (raw or "").strip()
    if not text:
        return {}
    if text.startswith("{"):
        data = json.loads(text)
        if not isinstance(data, dict):
            raise ValueError("Solver parameters JSON must be an object.")
        return data
    options = {}
    for part in text.split(","):
        item = part.strip()
        if not item:
            continue
        if "=" not in item:
            raise ValueError(
                "Solver parameters must look like timelimit=3600 or a JSON object."
            )
        key, value = item.split("=", 1)
        key = key.strip()
        if not key:
            raise ValueError("A solver parameter is missing its name.")
        options[key] = _coerce_option(value)
    return options


def _drop_thread_options(options):
    """Thread count is chosen from the CPU count, not from the user."""
    removed = False
    for key in ("threads", "Threads"):
        if key in options:
            options.pop(key)
            removed = True
    return removed


def ask_solver_selection():
    """Choose which installed solver to use, then optional solver parameters."""
    found = [name for name in os.environ.get("OPTME_SOLVERS_FOUND", "").split(",") if name]
    print()
    print("  The same solver is used for FSEOF, OptForce, and the enzyme LPs.")
    print("  Installed solvers: " + ", ".join(found))
    print("  Thread count is chosen from the number of CPUs.")
    print("  Extra parameters are optional. Example: timelimit=3600")
    print()

    preset = (
        os.environ.get("OPTME_SOLVER", "").strip().lower()
        or os.environ.get("OPTME_COBRA_SOLVER", "").strip().lower()
    )
    if not sys.stdin.isatty():
        if not preset:
            if len(found) == 1:
                preset = found[0]
            else:
                print("  ERROR: set OPTME_SOLVER to one of: " + ", ".join(found))
                return False
        if preset not in found:
            print(f"  ERROR: solver '{preset}' was not found under the solver path.")
            return False
        options_raw = os.environ.get("OPTME_SOLVER_OPTIONS", "")
        try:
            options = parse_solver_options(options_raw)
        except (ValueError, json.JSONDecodeError) as exc:
            print(f"  ERROR: {exc}")
            return False
    else:
        default_choice = preset if preset in found else (found[0] if len(found) == 1 else None)
        while True:
            choice = (prompt("  Solver", default_choice) or "").strip().lower()
            if choice in found:
                preset = choice
                break
            print("  Choose one of: " + ", ".join(found))
        while True:
            options_raw = prompt("  Extra solver parameters", "") or ""
            try:
                options = parse_solver_options(options_raw)
                break
            except (ValueError, json.JSONDecodeError) as exc:
                print(f"  ERROR: {exc}")

    if _drop_thread_options(options):
        print("  Thread count is set automatically and was not taken from the parameters.")
    if options:
        os.environ["OPTME_SOLVER_OPTIONS"] = json.dumps(options)
    else:
        os.environ.pop("OPTME_SOLVER_OPTIONS", None)
    n = use_sequential_threads()
    os.environ["OPTME_SOLVER"] = preset
    os.environ["OPTME_COBRA_SOLVER"] = preset
    os.environ["OPTME_PYOMO_SOLVER"] = preset
    print(f"  Using solver {preset} with parameters {os.environ['OPTME_SOLVER_OPTIONS']}")
    print(
        f"  CPUs: {n}. A single solve uses {n} threads. "
        f"Parallel enzyme LPs use {n} processes with 1 thread each."
    )
    return True


def ask_solver_path():
    """Ask for the solver path before any LLM settings. Required."""
    print(f"\n{SEP}")
    print("Step 1 / 4  —  Solver")
    print(SEP)
    print("  Give a CPLEX Studio directory, a Gurobi linux64 directory, or both")
    print("  separated by ':' . You then choose which solver to use and may set")
    print("  its parameters. This is requested before any LLM model, URL, or key.")
    print()
    preset = os.environ.get("OPTME_SOLVER_PATH", "").strip()
    if not sys.stdin.isatty():
        if not preset:
            print("  ERROR: set OPTME_SOLVER_PATH before a non-interactive run.")
            return None
        ok, message = configure_solver_path(preset)
        print("  " + message)
        return preset if ok else None

    while True:
        raw = prompt("  Solver path", preset)
        if not raw:
            print("  Solver path is required. Enter it before LLM settings.")
            continue
        ok, message = configure_solver_path(raw)
        print("  " + message)
        if ok:
            return raw


def ask_llm(task):
    """Ask for LLM model, URL, and key. Skip llm if any value is missing."""
    methods = task.get("taskname", [])
    if isinstance(methods, str):
        methods = [methods]
    if "llm" not in methods:
        return task

    print(f"\n{SEP}")
    print("Step 2 / 4  —  LLM settings")
    print(SEP)
    print("  The llm method needs a model name, API base URL, and API key.")
    print("  Leave any field blank to skip llm and run the other methods.")
    print()

    if not sys.stdin.isatty():
        model = os.environ.get("OPTME_LLM_MODEL", "").strip()
        url = os.environ.get("OPTME_LLM_BASE", "").strip()
        key = os.environ.get("OPTME_LLM_KEY", "").strip()
    else:
        model = prompt("  LLM model name", os.environ.get("OPTME_LLM_MODEL"))
        url = prompt("  LLM API base URL", os.environ.get("OPTME_LLM_BASE"))
        key = prompt("  LLM API key") or ""

    if model and url and key:
        os.environ["LLM_MODEL_NAME"] = model
        os.environ["LLM_API_BASE"] = url.rstrip("/")
        os.environ["LLM_API_KEY"] = key
        os.environ["OPTME_LLM_MODEL"] = model
        os.environ["OPTME_LLM_BASE"] = url.rstrip("/")
        print("  LLM method will run with the supplied model and URL.")
        return task

    print("  LLM model, URL, or key was not given. Skipping llm.")
    for name in ("LLM_MODEL_NAME", "LLM_API_BASE", "LLM_API_KEY"):
        os.environ.pop(name, None)
    task = dict(task)
    task["taskname"] = [m for m in methods if m != "llm"]
    return task


def configure_dirs(default_task_dir, default_model_dir, default_map_dir, default_output_dir):
    """Let the user confirm or change the four key directories."""
    print(f"\n{SEP}")
    print("Step 3 / 4  —  Confirm Directories")
    print(SEP)

    dirs = {
        "Task   dir (INPUT/task) ": default_task_dir,
        "Model  dir (INPUT/model)": default_model_dir,
        "Map    dir (INPUT/map)  ": default_map_dir,
        "Output dir              ": default_output_dir,
    }
    print("  Current paths:")
    for label, path in dirs.items():
        status = "OK" if os.path.exists(path) else "MISSING"
        print(f"    {label}  {path}  [{status}]")
    print()
    if not sys.stdin.isatty() or confirm("  Use these paths?"):
        return default_task_dir, default_model_dir, default_map_dir, default_output_dir

    task_dir = prompt("  Task   dir", default_task_dir)
    model_dir = prompt("  Model  dir", default_model_dir)
    map_dir = prompt("  Map    dir", default_map_dir)
    output_dir = prompt("  Output dir", default_output_dir)
    return task_dir, model_dir, map_dir, output_dir


def _task_value(task, key):
    value = task.get(key, "")
    if value is None or value == "" or value == []:
        return "(none)"
    return value


def review_and_confirm(job_id, task, task_dir, model_dir, map_dir, output_dir, methods):
    """Show the task inputs and ask for confirmation before running."""
    task_path = os.path.join(task_dir, TASK_FILE)
    excluded = task.get("excluded_rxns") or []
    if isinstance(excluded, (list, tuple)):
        excluded_text = "(none)" if not excluded else ", ".join(str(item) for item in excluded)
    else:
        excluded_text = str(excluded)
    print(f"\n{SEP}")
    print("Step 4 / 4  —  Review & Run")
    print(SEP)
    print(f"  Task file  : {task_path}")
    print("  Edit this file and start again if any of the following is not the case you want.")
    print(f"  Substrate  : {_task_value(task, 'substrate')} ({_task_value(task, 'substrate_name')})")
    print(f"  Uptake     : {_task_value(task, 'substrate_uptake_rate')} mmol/gDW/h")
    print(f"  Product    : {_task_value(task, 'product')} ({_task_value(task, 'product_name')})")
    print(f"  Biomass    : {_task_value(task, 'biomass')}")
    print(f"  Min growth : {_task_value(task, 'min_growth')}")
    print(f"  ATPM       : {_task_value(task, 'ATPM')}")
    print(f"  Oxygen     : {_task_value(task, 'oxygenstate')} ({_task_value(task, 'O2')})")
    print(f"  Model      : {_task_value(task, 'model')}")
    print(f"  Species    : {_task_value(task, 'species')} ({_task_value(task, 'ID')})")
    print(f"  Excluded   : {excluded_text}")
    print(f"  Methods    : {', '.join(methods)}")
    print(f"  Job ID     : {job_id}")
    print(f"  Solver path: {os.environ.get('OPTME_SOLVER_PATH', '')}")
    print(f"  Solver     : {os.environ.get('OPTME_SOLVER', '')}")
    print(f"  Parameters : {os.environ.get('OPTME_SOLVER_OPTIONS', '(defaults)')}")
    print(f"  Model dir  : {model_dir}")
    print(f"  Map dir    : {map_dir}")
    print(f"  Output dir : {output_dir}")
    print()
    if not sys.stdin.isatty():
        return True
    return confirm("  Start pipeline?")


def _read_bytes(path):
    with open(path, "rb") as handle:
        return handle.read()


def _write_bytes(path, data):
    with open(path, "wb") as handle:
        handle.write(data)


def main():
    print(f"\n{SEP}")
    print("  OptME  —  C. glutamicum iCW example")
    print(SEP)
    print("  Task file : " + TASK_FILE)

    task_path = os.path.join(TASK_DIR, TASK_FILE)
    if not os.path.isfile(task_path):
        print(f"\n  ERROR: {task_path} is missing.")
        return 1

    original_task = _read_bytes(task_path)
    model_name = json.loads(original_task.decode("utf-8")).get("model", "ala_iCW.json")
    model_path = os.path.join(MODEL_DIR, model_name)
    original_model = _read_bytes(model_path) if os.path.isfile(model_path) else None

    solver_path = ask_solver_path()
    if not solver_path or not ask_solver_selection():
        print("  Aborted: solver path or solver selection was not accepted.")
        return 1

    with open(task_path, encoding="utf-8") as handle:
        task = json.load(handle)
    task = ask_llm(task)
    methods = task.get("taskname", [])
    if isinstance(methods, str):
        methods = [methods]
    only_raw = os.environ.get("OPTME_METHODS", "").strip()
    if only_raw:
        wanted = [part.strip() for part in only_raw.split(",") if part.strip()]
        methods = [name for name in methods if name in wanted]
        task["taskname"] = methods
        print("  Methods limited to: " + ", ".join(methods))
    if not methods:
        print("  No methods left to run.")
        _write_bytes(task_path, original_task)
        return 1

    task_dir, model_dir, map_dir, output_dir = configure_dirs(
        TASK_DIR, MODEL_DIR, MAP_DIR, OUTPUT_DIR
    )
    active_task = os.path.join(task_dir, TASK_FILE)
    if not os.path.isfile(active_task) and active_task != task_path:
        print(f"\n  ERROR: {TASK_FILE} not found in {task_dir}. Aborting.")
        _write_bytes(task_path, original_task)
        return 1

    if not review_and_confirm(JOB_ID, task, task_dir, model_dir, map_dir, output_dir, methods):
        print("  Aborted by user.")
        _write_bytes(task_path, original_task)
        return 1

    with open(active_task, "w", encoding="utf-8") as handle:
        json.dump(task, handle)

    print(f"\n{SEP}")
    print("  Running Pipeline")
    print(SEP)

    os.makedirs(output_dir, exist_ok=True)
    os.environ.setdefault("MPLBACKEND", "Agg")

    log_file = os.path.join(output_dir, f"{JOB_ID}.log")
    if os.path.exists(log_file):
        os.remove(log_file)

    steps = [
        ("Step 1: Pre-check (Error.py)", [OPTME_PYTHON, os.path.join(OPTME_DIR, "Error.py"), model_dir, task_dir, JOB_ID]),
        ("Step 2: summary.py", [OPTME_PYTHON, os.path.join(OPTME_DIR, "summary.py"), model_dir, task_dir, map_dir, output_dir, JOB_ID]),
        ("Step 3: outputjson_status.py", [OPTME_PYTHON, os.path.join(OPTME_DIR, "outputjson_status.py"), task_dir, output_dir, JOB_ID]),
        ("Step 4: sum.py", [OPTME_PYTHON, os.path.join(OPTME_DIR, "sum.py"), model_dir, task_dir, map_dir, output_dir, JOB_ID]),
    ]
    step_times = []
    pipeline_started = time.perf_counter()
    ran_pipeline = False
    try:
        for label, command in steps:
            ran_pipeline = True
            print(f"\n  {label}...")
            ok, elapsed = run_command(command, log_file)
            step_times.append((label, elapsed, "ok" if ok else "failed"))
            if not ok:
                print(f"  Stopped. See {log_file}")
                return 1
    finally:
        if ran_pipeline:
            _append_runtime_summary(
                log_file, step_times, time.perf_counter() - pipeline_started
            )
        _write_bytes(task_path, original_task)
        if original_model is not None and os.path.isfile(model_path):
            _write_bytes(model_path, original_model)

    print(f"\n{SEP}")
    print(f"  Completed: {JOB_ID}")
    print(f"  Log    : {log_file}")
    print(f"  Output : {output_dir}")
    print(SEP)
    return 0


if __name__ == "__main__":
    sys.exit(main())
