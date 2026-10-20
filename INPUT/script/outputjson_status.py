# coding:utf-8
"""
Author: Jingyi Cai (2024-2026)
Function: Record whether each method produced a target prediction, and clear visualization paths for methods that failed.
Input: python outputjson_status.py <task_dir> <results_dir> <task_id>. Reads <task_dir>/<task_id>.json and <results_dir>/<task_id>/output.json.
Output: Updates <results_dir>/<task_id>/output.json with per-method status and an overall status.
"""
import os
import json
import sys
import pandas as pd

# Methods that have escher/d3flux visualization paths in output.json
METHODS_WITH_VIZ = {'FSEOF', 'loopless_optforce_MUST', 'iBridge', 'E_FSEOF', 'E_OptForce'}

def remove_path(output, method, type):
    """Blank visualization paths for a failed method. Only keys that exist are cleared."""
    if type == "d3flux":
        for suffix in ['wild', 'over', 'fc']:
            key = f"{method}_d3flux_{suffix}"
            if key in output.get("data", {}).get("Visualization", {}):
                output["data"]["Visualization"][key] = ''
    if type == "escher":
        for suffix in ['wild', 'over', 'fc']:
            key = f"{method}_escher_{suffix}"
            if key in output.get("data", {}).get("Visualization", {}):
                output["data"]["Visualization"][key] = ''
    return output

if __name__ == "__main__":
    path_task = sys.argv[1]
    path_results = sys.argv[2]
    taskname = sys.argv[3]
    path_results2 = os.path.join(path_results, taskname)

    if not os.path.exists(path_results2):
        os.makedirs(path_results2)

    with open(os.path.join(path_task, taskname) + '.json', encoding='utf-8') as fp:
        inputdic = json.load(fp)
    method = inputdic['taskname']
    if isinstance(method, str):
        method = [method]

    taskidmap = {
        'FSEOF': 'FSEOF',
        'loopless_optforce_MUST': 'OptForce',
        'iBridge': 'iBridge',
        'E_FSEOF': 'E_FSEOF',
        'E_OptForce': 'E_OptForce',
        'llm': 'llm'
    }

    checkdic = pd.DataFrame()
    checkdic['target_prediction'] = [0] * len(method)
    checkdic['errorstring'] = '' * len(method)
    checkdic.index = method

    for i in range(len(method)):
        m = method[i]
        m_key = taskidmap.get(m, m)
        output_path = os.path.join(path_results2, m_key, 'output.json')
        if os.path.exists(output_path):
            if os.path.getsize(output_path) > 0:
                checkdic.loc[m, 'target_prediction'] = 1
            else:
                checkdic.loc[m, 'errorstring'] += f' {m_key} target prediction failed;'
        else:
            checkdic.loc[m, 'errorstring'] += f' {m_key} target prediction failed;'

    if all(checkdic['target_prediction'] == 1):
        overall_status = 'success'
        errorstring = ''
    elif all(checkdic['target_prediction'] == 0):
        overall_status = 'failed'
        errorstring = checkdic['errorstring'].sum()
    else:
        overall_status = 'partial success'
        errorstring = checkdic['errorstring'].sum()

    checkdic2 = checkdic.to_dict(orient='index')

    combined_output_path = os.path.join(path_results2, 'output.json')
    if not os.path.isfile(combined_output_path):
        print(f"Warning: combined output.json not found at {combined_output_path}")
        sys.exit(1)

    with open(combined_output_path, 'r') as f:
        output = json.load(f)

    for i in checkdic2:
        if checkdic2[i]['target_prediction'] == 0:
            m_key = taskidmap.get(i, i)
            if i in METHODS_WITH_VIZ:
                output = remove_path(output, m_key, "escher")

    output['overall_status'] = overall_status
    output['errorstring'] = errorstring
    output['error_detail'] = checkdic2

    with open(combined_output_path, 'w') as f:
        json.dump(output, f)
