# coding:utf-8
"""Choose how many CPUs a solve may use.

N is the logical CPU count. If the machine does not report one, N is 10.
A step that solves one linear program at a time gives that program N threads.
A step that solves many independent linear programs at once uses N processes
and 1 thread in each, so the two layers do not both claim every core.
"""
import json
import os

_THREADS_ENV = "OPTME_CPLEX_THREADS"
_OPTIONS_ENV = "OPTME_SOLVER_OPTIONS"
_INSTALLED = False


def cpu_budget():
    """Return the number of CPUs this job should occupy."""
    raw = os.cpu_count()
    try:
        n = int(raw)
    except (TypeError, ValueError):
        n = 0
    if n < 1:
        return 10
    return n


def _load_options():
    raw = os.environ.get(_OPTIONS_ENV, "").strip()
    if not raw:
        return {}
    try:
        data = json.loads(raw)
    except json.JSONDecodeError:
        return {}
    if not isinstance(data, dict):
        return {}
    return data


def _store_threads(n):
    n = int(n)
    options = _load_options()
    options.pop("Threads", None)
    options["threads"] = n
    os.environ[_OPTIONS_ENV] = json.dumps(options)
    os.environ[_THREADS_ENV] = str(n)
    return n


def current_threads():
    """Thread count for the solve that is about to run."""
    raw = os.environ.get(_THREADS_ENV, "").strip()
    try:
        n = int(raw)
    except ValueError:
        n = 0
    if n < 1:
        return cpu_budget()
    return n


def use_sequential_threads():
    """One solver at a time: it may use every CPU."""
    return _store_threads(cpu_budget())


def use_parallel_lp_threads():
    """Many independent solvers at once: each one uses a single CPU."""
    return _store_threads(1)


def options_for_solver(solver):
    """Solver options for this process, with the thread count for this step."""
    options = _load_options()
    options.pop("threads", None)
    options.pop("Threads", None)
    n = current_threads()
    if solver and "gurobi" in str(solver).lower():
        options["Threads"] = n
    else:
        options["threads"] = n
    return options


def apply_cobra_threads(model):
    """Set the live solver's thread count. Ignore interfaces that have none."""
    n = current_threads()
    problem = model.solver.problem
    try:
        problem.parameters.threads.set(n)
        return
    except Exception:
        pass
    try:
        problem.Params.Threads = n
    except Exception:
        pass


def install_cobra_threads():
    """Make each COBRA solve pick up the thread count for this process."""
    global _INSTALLED
    if _INSTALLED:
        return
    import cobra

    def _wrap(name):
        original = getattr(cobra.Model, name)

        def wrapped(self, *args, **kwargs):
            apply_cobra_threads(self)
            return original(self, *args, **kwargs)

        wrapped.__name__ = getattr(original, "__name__", name)
        wrapped.__doc__ = getattr(original, "__doc__", None)
        setattr(cobra.Model, name, wrapped)

    _wrap("optimize")
    _wrap("slim_optimize")
    _INSTALLED = True


def parallel_fva(model, **kwargs):
    """Run flux variability with N processes and 1 thread in each."""
    from cobra.flux_analysis.variability import flux_variability_analysis

    kwargs.setdefault("processes", cpu_budget())
    use_parallel_lp_threads()
    try:
        return flux_variability_analysis(model, **kwargs)
    finally:
        use_sequential_threads()
