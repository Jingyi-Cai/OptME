# coding:utf-8
"""
Author: Jingyi Cai (2024-2026)
Function: Predict metabolic engineering targets with an LLM. Skips the run when the model name, API base URL, or API key is missing.
Input: python llmoptme.py <model_dir> <task_dir> <map_dir> <results_dir> <task_id>. Reads the task JSON, the model, script/llmsetup.md, and environment variables LLM_MODEL_NAME, LLM_API_BASE, and LLM_API_KEY.
Output: <results_dir>/<task_id>/llm/output.json. Credentials are read from the environment and are not written to disk.
"""
import sys
import os
import json
import re

import cobra
import requests
import json_repair

from task_config import LLM_API_BASE, LLM_API_KEY, LLM_MODEL_NAME, LLM_SETUP_FILE


def extract_assistant_text(resp):
    """Extract text content from an OpenAI-compatible chat response."""
    try:
        obj = resp.json()
        msg = obj["choices"][0]["message"]
        content = msg.get("content") or ""
        reasoning = msg.get("reasoning") or ""
        # Strip <think>...</think> blocks emitted by reasoning models (e.g. MiniMax-M2.7)
        import re as _re
        content = _re.sub(r"<think>.*?</think>", "", content, flags=_re.DOTALL).strip()
        return (content or reasoning).strip(), obj
    except Exception as exc:
        try:
            body = resp.json()
        except Exception:
            body = resp.text[:500] if hasattr(resp, "text") else str(resp)
        print(f"[llmoptme] extract_assistant_text error: {exc}; status={getattr(resp, 'status_code', '?')}; body={body}", flush=True)
        return "", None


def normalize_modification(value):
    if not isinstance(value, str):
        return value
    v = value.strip().upper()
    if v in ("UP", "OVER", "OVEREXPRESSION"):
        return "UP"
    if v in ("DOWN", "DOWNREGULATE", "DOWN-REGULATED", "DOWN_REGULATED"):
        return "Down"
    if v in ("DEL", "DELETE", "DELETION", "KO", "KNOCKOUT"):
        return "DEL"
    return value


def clean_text(value):
    if not isinstance(value, str):
        return ""
    text = value.strip()
    if text.lower() in {"", "none", "null", "n/a", "na", "gene", "reaction"}:
        return ""
    return text


def choose_candidate_key(item, index):
    raw_key = clean_text(str(item.get("gene_or_reaction") or item.get("id") or ""))
    gene_name = clean_text(item.get("gene_name") or item.get("Gene_Name") or "")
    enzyme_name = clean_text(item.get("enzyme_name") or item.get("Enzyme_Name") or "")
    if raw_key:
        return raw_key
    if gene_name:
        return gene_name
    if enzyme_name:
        return enzyme_name
    return f"candidate_{index + 1}"


def to_genes_dict(parsed_obj):
    """Convert LLM JSON to a normalised {gene_key: info} dict."""
    if not isinstance(parsed_obj, dict):
        return {}

    if isinstance(parsed_obj.get("candidates"), list):
        genes = {}
        for i, item in enumerate(parsed_obj["candidates"]):
            if not isinstance(item, dict):
                continue
            key = choose_candidate_key(item, i)
            enzyme_name = clean_text(item.get("enzyme_name") or item.get("Enzyme_Name", ""))
            gene_name = clean_text(item.get("gene_name") or item.get("Gene_Name") or key) or key
            genes[key] = {
                "Reaction Formula": item.get("Reaction Formula", ""),
                "Type": item.get("Type", "Enzymatic"),
                "Modification": normalize_modification(item.get("modification") or item.get("Modification", "UP")),
                "Fold_change": item.get("Fold_change", 0),
                "Reaction_Name": item.get("Reaction_Name") or key,
                "Rationale": item.get("rationale") or item.get("Rationale", ""),
                "Gene_Name": gene_name,
                "Enzyme_Name": enzyme_name,
            }
        return genes

    genes = {}
    for k, v in parsed_obj.items():
        if not isinstance(v, dict):
            continue
        if "Modification" not in v:
            continue
        gene_name = clean_text(v.get("Gene_Name") or k) or k
        genes[k] = {
            "Reaction Formula": v.get("Reaction Formula", ""),
            "Type": v.get("Type", "Enzymatic"),
            "Modification": normalize_modification(v.get("Modification", "UP")),
            "Fold_change": v.get("Fold_change", 0),
            "Reaction_Name": v.get("Reaction_Name") or k,
            "Rationale": v.get("Rationale", ""),
            "Gene_Name": gene_name,
            "Enzyme_Name": clean_text(v.get("Enzyme_Name", "")),
        }
    return genes


def summarize_genes(genes):
    summary = {"UP": 0, "Down": 0, "KO": 0}
    for v in genes.values():
        m = v.get("Modification")
        if m == "UP":
            summary["UP"] += 1
        elif m == "Down":
            summary["Down"] += 1
        elif m == "DEL":
            summary["KO"] += 1
    return summary


def _chat_completions_url(base):
    """Accept either https://host/v1 or a host without the version prefix."""
    base = (base or "").rstrip("/")
    if base.endswith("/chat/completions"):
        return base
    if base.endswith("/v1"):
        return base + "/chat/completions"
    return base + "/v1/chat/completions"


def call_llm(setup_text, inputdic):
    """Send the metabolic engineering prompt to the LLM and return parsed result."""
    url = _chat_completions_url(LLM_API_BASE)
    headers = {
        "Authorization": f"Bearer {LLM_API_KEY}",
        "Content-Type": "application/json",
    }

    task_user_prompt = (
        f"Overproduce {inputdic['product_name']} in {inputdic['species']} "
        f"under {inputdic['oxygenstate']} using {inputdic['substrate_name']}. "
        "Output ONLY one valid JSON object with key 'candidates' (array). "
        "Return 3 to 8 realistic metabolic engineering candidates. "
        "Each candidate must include: gene_or_reaction, gene_name, enzyme_name, "
        "modification(UP/Down/DEL), rationale. "
        "Use actual gene symbols for gene_or_reaction and gene_name. "
        "Never use placeholders like Gene, Reaction, None, null, or candidate_1. "
        "Prefer real enzyme names such as glutamate dehydrogenase instead of "
        "repeating the gene symbol. "
        "Do not include markdown, comments, or any extra text outside the JSON."
    )

    payload = {
        "model": LLM_MODEL_NAME,
        "messages": [
            {
                "role": "system",
                "content": setup_text + "\nOutput only final JSON. No thinking text.",
            },
            {"role": "user", "content": task_user_prompt},
        ],
        "temperature": 0,
        "max_tokens": 2000,
        "response_format": {"type": "text"},
    }

    import time as _time
    _max_retries = int(os.environ.get("LLM_MAX_RETRIES", "12"))
    _retry_delay = int(os.environ.get("LLM_RETRY_DELAY", "20"))
    _timeout = int(os.environ.get("LLM_REQUEST_TIMEOUT", "180"))

    for _attempt in range(1, _max_retries + 1):
        resp = requests.post(url, json=payload, headers=headers, timeout=_timeout)
        if resp.status_code == 200:
            break
        try:
            err_body = resp.json()
        except Exception:
            err_body = resp.text[:500]

        # New API frequently returns 503 when shared CPU is overloaded.
        # Retry with backoff so batch tasks can eventually get a response.
        if resp.status_code in (400, 422) and "response_format" in payload:
            print("[llmoptme] API rejected response_format; retrying without it.", flush=True)
            payload.pop("response_format", None)
            continue
        if resp.status_code in (429, 503) and _attempt < _max_retries:
            _sleep_s = _retry_delay * _attempt
            print(
                f"[llmoptme] HTTP {resp.status_code} (attempt {_attempt}/{_max_retries}), retrying in {_sleep_s}s: {err_body}",
                flush=True,
            )
            _time.sleep(_sleep_s)
            continue

        raise RuntimeError(
            f"LLM API error: HTTP {resp.status_code} — {err_body}"
        )
    txt, obj = extract_assistant_text(resp)

    parsed_obj = {}
    try:
        parsed_obj = json_repair.loads(txt)
    except Exception:
        pass

    genes = to_genes_dict(parsed_obj)
    return {
        "summary": summarize_genes(genes),
        "genes": genes,
        "modelname": LLM_MODEL_NAME,
    }


def normalize_lookup_text(value):
    if not isinstance(value, str):
        return ""
    text = value.strip().lower()
    text = re.sub(r"\s+", " ", text)
    return text


def first_formula_from_reactions(reactions):
    if not isinstance(reactions, dict):
        return ""
    for reaction_info in reactions.values():
        if isinstance(reaction_info, dict):
            formula = (reaction_info.get("Reaction Formula") or "").strip()
            if formula:
                return formula
    return ""


def extract_from_rationale(text):
    patterns = [
        r"activity of ([A-Za-z0-9\-\(\) ,]+?) \(",
        r"Upregulating the ([A-Za-z0-9\-\(\) ,]+?) \(",
        r"Increasing the expression of ([A-Za-z0-9\-\(\) ,]+?) \(",
        r"Enhancing the activity of ([A-Za-z0-9\-\(\) ,]+?) \(",
    ]
    for pat in patterns:
        match = re.search(pat, text, re.IGNORECASE)
        if match:
            return match.group(1).strip()
    return ""


def enrich_with_mappings(llm_result, mapping4, model):
    """Fill in Reaction Formula, Gene_Name, and Enzyme_Name from mapping data and model."""
    reaction_name_by_id = {rxn.id: rxn.name for rxn in model.reactions}
    formula_to_name = {
        rxn.build_reaction_string(use_metabolite_names=False): rxn.name
        for rxn in model.reactions
    }

    gene_to_formula_4 = {}
    gene_to_formula_4_lower = {}
    enzyme_to_formula_4 = {}
    for gene_key, gene_entry in mapping4.items():
        if not isinstance(gene_entry, dict):
            continue
        formula = first_formula_from_reactions(gene_entry.get("reactions"))
        if not formula:
            continue
        gene_name = (gene_entry.get("gene_name") or gene_key or "").strip()
        if gene_key:
            gene_to_formula_4.setdefault(gene_key, formula)
            gene_to_formula_4_lower.setdefault(gene_key.lower(), formula)
        if gene_name:
            gene_to_formula_4.setdefault(gene_name, formula)
            gene_to_formula_4_lower.setdefault(gene_name.lower(), formula)
        enzyme_name = normalize_lookup_text(gene_entry.get("enzyme_name") or "")
        if enzyme_name:
            enzyme_to_formula_4.setdefault(enzyme_name, formula)

    for gene_key, gene_info in llm_result.get("genes", {}).items():
        if not isinstance(gene_info, dict):
            continue
        map_entry = mapping4.get(gene_key, {}) if isinstance(mapping4.get(gene_key), dict) else {}

        gene_name = (gene_info.get("Gene_Name") or "").strip()
        if not gene_name:
            gene_name = map_entry.get("gene_name") or gene_key
            gene_info["Gene_Name"] = gene_name

        enzyme_name = (gene_info.get("Enzyme_Name") or "").strip()
        if not enzyme_name:
            enzyme_name = map_entry.get("enzyme_name") or ""

        current_formula = str(gene_info.get("Reaction Formula", "") or "").strip()
        if not current_formula:
            formula = ""
            gene_key_l = gene_key.lower()
            gene_name_l = gene_name.lower() if isinstance(gene_name, str) else ""
            enzyme_name_l = normalize_lookup_text(enzyme_name)

            if gene_key in gene_to_formula_4:
                formula = gene_to_formula_4[gene_key]
            if not formula and gene_key_l in gene_to_formula_4_lower:
                formula = gene_to_formula_4_lower[gene_key_l]
            if not formula and gene_name in gene_to_formula_4:
                formula = gene_to_formula_4[gene_name]
            if not formula and gene_name_l in gene_to_formula_4_lower:
                formula = gene_to_formula_4_lower[gene_name_l]
            if not formula and enzyme_name_l in enzyme_to_formula_4:
                formula = enzyme_to_formula_4[enzyme_name_l]

            if formula:
                gene_info["Reaction Formula"] = formula

        if not enzyme_name:
            reaction_formula = (gene_info.get("Reaction Formula") or "").strip()
            if ":" in reaction_formula:
                prefix, rest = reaction_formula.split(":", 1)
                prefix = prefix.strip()
                if prefix in reaction_name_by_id and reaction_name_by_id[prefix]:
                    enzyme_name = reaction_name_by_id[prefix]
                    gene_info["Reaction Formula"] = rest.strip()

        if not enzyme_name:
            reaction_formula = (gene_info.get("Reaction Formula") or "").strip()
            if reaction_formula in formula_to_name:
                enzyme_name = formula_to_name[reaction_formula]

        if not enzyme_name:
            enzyme_name = extract_from_rationale(gene_info.get("Rationale", ""))

        if enzyme_name:
            gene_info["Enzyme_Name"] = enzyme_name
        else:
            gene_info.setdefault("Enzyme_Name", "")


def main():
    path_model = sys.argv[1]
    path_task = sys.argv[2]
    path_map = sys.argv[3]
    path_results = sys.argv[4]
    taskname = sys.argv[5]

    if not (LLM_API_BASE and LLM_API_KEY and LLM_MODEL_NAME):
        print("[llmoptme] LLM model, URL, or key was not provided. Skipping llm.")
        return

    task_json_path = os.path.join(path_task, taskname) + ".json"
    with open(task_json_path, encoding="utf-8") as fp:
        inputdic = json.load(fp)

    setup_text = ""
    if os.path.isfile(LLM_SETUP_FILE):
        with open(LLM_SETUP_FILE, encoding="utf-8") as fp:
            setup_text = fp.read()
    else:
        print(f"Warning: LLM setup file not found: {LLM_SETUP_FILE}, using empty system prompt.")

    print(f"[llmoptme] Calling LLM ({LLM_MODEL_NAME}) for {inputdic.get('product_name', '')} ...")
    llm_result = call_llm(setup_text, inputdic)
    print(f"[llmoptme] LLM returned {len(llm_result.get('genes', {}))} candidates.")

    mapping4_path = os.path.join(path_model, "generxnmapping4.json")
    if os.path.isfile(mapping4_path):
        with open(mapping4_path, encoding="utf-8") as fp:
            mapping4 = json.load(fp)
    else:
        print(f"Warning: Mapping file not found: {mapping4_path}, skipping enrichment.")
        mapping4 = {}

    model_file = os.path.join(path_model, inputdic["model"])
    if ".mat" in inputdic["model"]:
        model = cobra.io.load_matlab_model(model_file)
    elif ".xml" in inputdic["model"] or ".sbml" in inputdic["model"]:
        model = cobra.io.read_sbml_model(model_file)
    elif ".json" in inputdic["model"]:
        model = cobra.io.load_json_model(model_file)
    else:
        model = cobra.io.load_model(inputdic["model"])

    if mapping4:
        enrich_with_mappings(llm_result, mapping4, model)
        llm_result["summary"] = summarize_genes(llm_result.get("genes", {}))

    output_dir = os.path.join(path_results, taskname, "llm")
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "output.json")
    with open(output_path, "w", encoding="utf-8") as fp:
        json.dump(llm_result, fp, indent=4, ensure_ascii=False)

    print(f"[llmoptme] Results saved to {output_path}")


if __name__ == "__main__":
    main()
