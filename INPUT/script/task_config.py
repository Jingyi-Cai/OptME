# coding:utf-8
"""
Author: Jingyi Cai (2024-2026)
Function: Path and environment settings for the iCW example. Directories are relative to the example root. run_icw_example.py changes to that root before the pipeline starts.
Input: Optional environment variables OPTME_DATA_DIR, OPTME_PYTHON, OPTME_TEMP_DIR, LLM_API_BASE, LLM_API_KEY, and LLM_MODEL_NAME.
Output: Path constants imported by the other scripts. This module does not write result files.
"""
import os
import sys

# Relative to the example root.
OPTME_DIR = "script"
DATA_DIR = os.environ.get("OPTME_DATA_DIR", ".")
OPTME_PYTHON = os.environ.get("OPTME_PYTHON", sys.executable)
IBRIDGE_PYTHON = OPTME_PYTHON

TARGET_SCRIPT = os.path.join(OPTME_DIR, "target.py")
IBRIDGE_SCRIPT = os.path.join(OPTME_DIR, "iBridge.py")
ETOPTME_SCRIPT = os.path.join(OPTME_DIR, "EToptme.py")
LLMOPTME_SCRIPT = os.path.join(OPTME_DIR, "llmoptme.py")

LLM_PYTHON = OPTME_PYTHON
LLM_API_BASE = os.environ.get("LLM_API_BASE", "").rstrip("/")
LLM_API_KEY = os.environ.get("LLM_API_KEY", "")
LLM_MODEL_NAME = os.environ.get("LLM_MODEL_NAME", "")
LLM_SETUP_FILE = os.path.join(OPTME_DIR, "llmsetup.md")

TEMP_DIR = os.environ.get("OPTME_TEMP_DIR", "output/tmp")
RESULTS_DIR = os.path.join(OPTME_DIR, "results")
LPFILES_DIR = os.path.join(OPTME_DIR, "lpfiles")
ECOECM_PROTAINMODEL_LP = os.path.join(RESULTS_DIR, "EcoECM_protainmodel.lp")
MDF_LP_FILE = os.path.join(LPFILES_DIR, "mdf.lp")
GEFBA_LP_FILE = os.path.join(LPFILES_DIR, "GEFBA.lp")

INPUT_DIR = "INPUT" if DATA_DIR in (".", "") else os.path.join(DATA_DIR, "INPUT")
OUTPUT_DIR = "output" if DATA_DIR in (".", "") else os.path.join(DATA_DIR, "output")
MODEL_DIR = os.path.join(INPUT_DIR, "model")
TASK_DIR = os.path.join(INPUT_DIR, "task")
MAP_DIR = os.path.join(INPUT_DIR, "map")
