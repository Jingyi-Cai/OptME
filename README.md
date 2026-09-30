# OptME example

Standalone example of the OptME single-task pipeline for the *C. glutamicum* genome-scale model iCW773. The case produces fumarate from aerobic D-glucose and runs FSEOF, loopless OptForce (MUST), enzyme-constrained FSEOF, enzyme-constrained OptForce, and an optional LLM target method.

## Layout


| Path           | Contents                                                                                       |
| -------------- | ---------------------------------------------------------------------------------------------- |
| `script/`      | Pipeline scripts and `run_icw_example.py`                                                      |
| `INPUT/task/`  | `iCW773R_task.json` and the carbon-metabolite table                                            |
| `INPUT/model/` | `ala_iCW.json`, kcat/MW, thermodynamic data, gene-reaction map, and cached irreversible models |
| `INPUT/map/`   | Escher central-metabolism map                                                                  |
| `output/`      | Logs and method results                                                                        |


## Task setup

The example task is `INPUT/task/iCW773R_task.json`. It asks for fumarate from aerobic D-glucose on the *C. glutamicum* model `ala_iCW.json`. The job id is the file name without `.json` (`iCW773R_task`), and results are written to `output/iCW773R_task/`.

`taskname` is one method or a list of methods. This example runs FSEOF, loopless OptForce, the two enzyme-constrained methods, and `llm`.

```json
{
  "substrate": "EX_glc_e",
  "substrate_name": "D-Glucose",
  "product": "EX_fum_e",
  "product_name": "Fumarate",
  "min_growth": 10,
  "substrate_uptake_rate": 10,
  "ATPM": 6.86,
  "oxygenstate": "aerobic",
  "O2": "EX_o2_e",
  "CO2_flag": "False",
  "model": "ala_iCW.json",
  "biomass": "CG_biomass_cgl_ATCC13032",
  "species": "C.glutamicum",
  "ID": "iCW773R",
  "email": "user_without_email",
  "excluded_rxns": [],
  "yield_na": false,
  "yield_na_reason": "",
  "taskname": ["FSEOF", "loopless_optforce_MUST", "E_FSEOF", "E_OptForce", "llm"]
}
```


| Field                         | Example                    | Role                                              |
| ----------------------------- | -------------------------- | ------------------------------------------------- |
| `substrate`, `substrate_name` | `EX_glc_e`, D-Glucose      | Uptake reaction and its name                      |
| `product`, `product_name`     | `EX_fum_e`, Fumarate       | Target exchange reaction and its name             |
| `substrate_uptake_rate`       | `10`                       | Substrate uptake bound (mmol/gDW/h)               |
| `min_growth`                  | `10`                       | Minimum growth setting stored with the task       |
| `ATPM`                        | `6.86`                     | ATP maintenance bound                             |
| `oxygenstate`, `O2`           | `aerobic`, `EX_o2_e`       | Oxygen condition and the oxygen exchange reaction |
| `model`                       | `ala_iCW.json`             | Model file in `INPUT/model/`                      |
| `biomass`                     | `CG_biomass_cgl_ATCC13032` | Biomass reaction                                  |
| `species`, `ID`               | `C.glutamicum`, `iCW773R`  | Organism label and model id                       |
| `excluded_rxns`               | `[]`                       | Reactions left out of target prediction           |
| `taskname`                    | list above                 | Methods to run                                    |


## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Install [IBM CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio), [Gurobi](https://www.gurobi.com/), or both. They may need a license (e.g. `GRB_LICENSE_FILE`) that matches the machine. The runner asks which of the installed solvers to use. That choice is used for FSEOF, OptForce, and the enzyme LPs. The CPLEX Python package must match the interpreter you use to launch the example.

## Run

From the example directory, with CPLEX and four threads. Clearing the LLM variables skips the `llm` method. The other methods still run.

```bash
cd "/path/to/OptME" && \
unset OPTME_LLM_MODEL OPTME_LLM_BASE OPTME_LLM_KEY LLM_API_KEY LLM_API_BASE LLM_MODEL_NAME && \
OPTME_SOLVER_PATH="/opt/ibm/ILOG/CPLEX_Studio221" \
OPTME_SOLVER="cplex" \
OPTME_SOLVER_OPTIONS="threads=4" \
MPLBACKEND=Agg \
/path/to/python_interpreter -u script/run_icw_example.py
```

`OPTME_SOLVER_PATH` is a CPLEX Studio directory, a Gurobi `linux64` directory, or both separated by `:`. `OPTME_SOLVER` is `cplex` or `gurobi`. `OPTME_SOLVER_OPTIONS` is optional; `threads=4,timelimit=3600` is the same form. To run `llm` as well, export `OPTME_LLM_MODEL`, `OPTME_LLM_BASE`, and `OPTME_LLM_KEY` instead of unsetting them. Those values stay in the process environment and are not written into the repository.

In a terminal the script still asks you to confirm the solver, the parameters, the directories, and the start of the run, using the variables above as the defaults. Results and the runtime summary are written to `output/iCW773R_task/` and `output/iCW773R_task.log`.

## Methods


| Task name                | Script               | Solver                                     |
| ------------------------ | -------------------- | ------------------------------------------ |
| `FSEOF`                  | `script/target.py`   | selected solver                            |
| `loopless_optforce_MUST` | `script/target.py`   | selected solver                            |
| `E_FSEOF`                | `script/EToptme.py`  | selected solver                            |
| `E_OptForce`             | `script/EToptme.py`  | selected solver                            |
| `llm`                    | `script/llmoptme.py` | none (skipped without model, URL, and key) |


`script/Error.py` checks the model, then `summary.py` runs each method, `outputjson_status.py` records status, and `sum.py` builds the combined tables.

`script/modelbuilder.py` holds the model-building helpers imported by `target.py`, `EToptme.py`, and `summary.py`.

## Model files used by this case

- `INPUT/model/ala_iCW.json` — iCW773 reconstruction used by the task
- `INPUT/model/reaction_change_by_enzuse_PDH_n.csv` — kcat and molecular weight
- `INPUT/model/iCW773_uniprot_modification_del_reaction_g0.csv` — reaction ΔG°
- `INPUT/model/metabolites_lnC_cg1.txt` — metabolite concentration bounds
- `INPUT/model/generxnmapping4.json` — gene and enzyme names for the LLM method
- `INPUT/model/iCW773R_task_model_irrev.json` and `iCW773R_task_model_irr_splited.json` — cached irreversible and isoenzyme-split models
- `INPUT/map/iJO1366.Central_metabolism.json` — map used for pathway figures
- `INPUT/task/met_contain_C_df.tsv` — metabolites that contain carbon
- `INPUT/task/reaction_groups_mapping.csv` — reaction-to-pathway map used by `sum.py`

## License

MIT. See [LICENSE](LICENSE). The iCW773 reconstruction and the Escher map remain subject to the terms of their original sources.