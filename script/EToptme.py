# coding:utf-8
"""
Author: Jingyi Cai, Wenqi Xu (2024-2026)
Function: Predict targets with enzyme-constrained FSEOF (E_FSEOF) or enzyme-constrained OptForce (E_OptForce). Linear programs use the solver named by OPTME_PYOMO_SOLVER and the options in OPTME_SOLVER_OPTIONS.
Input: python EToptme.py <model_dir> <task_dir> <map_dir> <results_dir> <task_id>. Reads the task JSON (including substrate_uptake_rate), ala_iCW.json, reaction ΔG, metabolite lnC bounds, and kcat/MW.
Output: <results_dir>/<task_id>/E_FSEOF/ or E_OptForce/, including output.json, results.xlsx, enzyme range files, and target figures.
"""
import sys
# sys.path.append('/hpcfs/fhome/xuwenqi/project/ET-OptME_websit')
from modelbuilder import *
import pandas as pd
import cobra
import json
import multiprocessing
import os
from multiprocessing import Pool
from sympy import subsets
import pandas as pd
import re
import matplotlib
matplotlib.use(os.environ.get("MPLBACKEND", "Agg"))
import matplotlib.pyplot as plt


def get_constraint_params(model_name):
    """Return model-specific defaults for E_total, B_value and K_value."""
    b_value = 0
    k_value = 1249
    model_name = str(model_name)

    if 'iCW773' in model_name or 'iCW' in model_name:
        e_total = 0.56 * np.mean([0.45311986236929197, 0.4622348377433211, 0.4600801040374112])
    elif 'iML1515' in model_name or 'iML' in model_name:
        e_total = 0.228
    elif 'ecBSU1' in model_name or 'BSU' in model_name:
        e_total = 0.28
    elif 'eciZM547' in model_name or 'ZM' in model_name:
        #  https://github.com/AaroncrowAries/eciZM547
        e_total = 0.454 * 0.605 * 0.5
    elif 'ecMTM' in model_name or 'MTM' in model_name:
        # https://link.springer.com/article/10.1186/s12934-024-02415-z
        e_total = 0.55 * 0.4653 * 0.5
    else:
        raise ValueError(f"Unsupported model for constraint defaults: {model_name}")

    return e_total, b_value, k_value  

def load_model(path_model, model_filename):
    # Determine the model type and load accordingly
    if '.mat' in model_filename:
        return cobra.io.load_matlab_model(os.path.join(path_model, model_filename))
    elif '.xml' in model_filename or '.sbml' in model_filename:
        return cobra.io.read_sbml_model(os.path.join(path_model, model_filename))
    elif '.json' in model_filename:
        return cobra.io.load_json_model(os.path.join(path_model, model_filename))
    else:
        return cobra.io.load_model(os.path.join(path_model, model_filename))
    
def load_input_dictionary(path_task, taskname):
    with open(os.path.join(path_task, taskname) + '.json', encoding='utf-8') as fp:
        return json.load(fp)


def disable_thermo_data(Concretemodel_Need_Data):
    """Clear thermo-related tables when g0/lnC data is unavailable or unmapped."""
    rg0_cols = []
    lnc_cols = ['lnClb', 'lnCub']

    if 'reaction_g0' in Concretemodel_Need_Data and isinstance(Concretemodel_Need_Data['reaction_g0'], pd.DataFrame):
        rg0_cols = list(Concretemodel_Need_Data['reaction_g0'].columns)
    if 'metabolites_lnC' in Concretemodel_Need_Data and isinstance(Concretemodel_Need_Data['metabolites_lnC'], pd.DataFrame):
        existing = list(Concretemodel_Need_Data['metabolites_lnC'].columns)
        if existing:
            lnc_cols = existing

    if not rg0_cols:
        rg0_cols = ['g0']

    Concretemodel_Need_Data['reaction_g0'] = pd.DataFrame(columns=rg0_cols)
    Concretemodel_Need_Data['metabolites_lnC'] = pd.DataFrame(columns=lnc_cols)
    Concretemodel_Need_Data['use_thermo_data'] = False
    return Concretemodel_Need_Data


def get_fba_mode(Concretemodel_Need_Data):
    """Use SET only when thermo data is usable; otherwise fall back to SE."""
    return 'SET' if Concretemodel_Need_Data.get('use_thermo_data', True) else 'SE'


def get_substrate_uptake_id(inputdic):
    """Resolved uptake reaction on the irreversible/split model (set in prepare_model)."""
    if inputdic.get('substrate_uptake'):
        return inputdic['substrate_uptake']
    substrate = inputdic['substrate']
    if str(substrate).endswith('_reverse'):
        return substrate
    return substrate + '_reverse'


def _normalize_substrate_exchange_ids(substrate_id):
    """Return (base_exchange_id, reverse_exchange_id) without double _reverse suffix."""
    substrate_id = str(substrate_id).strip()
    if substrate_id.endswith('_reverse'):
        base_id = substrate_id[:-len('_reverse')]
        return base_id, substrate_id
    return substrate_id, substrate_id + '_reverse'


def _exchange_external_metabolite(reaction):
    """Return the extracellular metabolite for a transport/exchange reaction."""
    for met in reaction.metabolites:
        comp = getattr(met, 'compartment', '') or ''
        if comp == 'e' or comp.endswith('_e') or 'extracellular' in comp.lower():
            return met
    mets = list(reaction.metabolites)
    return mets[0] if mets else None


def _reaction_uptake_flux_sign(reaction):
    """
    Infer uptake flux sign on an irreversible (or effective one-way) exchange.

    Returns 'positive' if positive flux imports substrate into the model,
    'negative' if negative flux imports substrate, else None.
    """
    for met, coeff in reaction.metabolites.items():
        comp = getattr(met, 'compartment', '') or ''
        if comp == 'e' or comp.endswith('_e') or 'extracellular' in comp.lower():
            if coeff > 0:
                return 'positive'
            if coeff < 0:
                return 'negative'
    for met, coeff in reaction.metabolites.items():
        if coeff > 0:
            return 'positive'
        if coeff < 0:
            return 'negative'
    return None


def _reaction_allows_uptake(reaction, direction, min_uptake=1e-9):
    min_uptake = float(min_uptake)
    if direction == 'positive':
        return float(reaction.upper_bound) >= min_uptake
    if direction == 'negative':
        return float(reaction.lower_bound) <= -min_uptake
    return False


def _configure_uptake_bounds(reaction, substrate_uptake_rate, direction):
    """Set exchange bounds so uptake up to substrate_uptake_rate is allowed."""
    conc = float(substrate_uptake_rate)
    if direction == 'positive':
        reaction.lower_bound = 0
        reaction.upper_bound = conc
    elif direction == 'negative':
        reaction.lower_bound = -conc
        reaction.upper_bound = 0
    return reaction


def _guess_extracellular_metabolite(model, exchange_base_id):
    """Infer extracellular metabolite from an exchange ID or existing reaction."""
    if exchange_base_id in model.reactions:
        met = _exchange_external_metabolite(model.reactions.get_by_id(exchange_base_id))
        if met is not None:
            return met
    base = exchange_base_id
    if base.endswith('_reverse'):
        base = base[:-len('_reverse')]
    stem = base[3:] if base.startswith('EX_') else base
    met_candidates = [
        stem,
        stem.replace('_e', '__D_e'),
        stem.replace('glc_e', 'glc__D_e'),
        stem.replace('__L_e', '__L_e'),
    ]
    seen = set()
    for met_id in met_candidates:
        if met_id in seen:
            continue
        seen.add(met_id)
        if met_id in model.metabolites:
            return model.metabolites.get_by_id(met_id)
    raise KeyError(
        f"Cannot infer extracellular metabolite for exchange '{exchange_base_id}'. "
        f"Tried metabolite IDs: {sorted(seen)}"
    )


def _add_substrate_uptake_reaction(model, exchange_base_id, extracellular_met, substrate_uptake_rate):
    """
    Add an irreversible import reaction (positive flux = uptake) for split models.
    Stoichiometry: --> met_e  (same as EX_*_reverse on iCW models).
    """
    uptake_id = f"{exchange_base_id}_uptake_import"
    if uptake_id in model.reactions:
        reaction = model.reactions.get_by_id(uptake_id)
    else:
        reaction = cobra.Reaction(uptake_id)
        reaction.name = f"Substrate uptake import for {exchange_base_id}"
        reaction.lower_bound = 0
        reaction.upper_bound = float(substrate_uptake_rate)
        reaction.add_metabolites({extracellular_met: 1.0})
        model.add_reactions([reaction])
        print(f"[substrate_uptake] Added uptake reaction {uptake_id}: {reaction.reaction}")
    _configure_uptake_bounds(reaction, substrate_uptake_rate, 'positive')
    return uptake_id, 'positive'


def resolve_substrate_uptake_reaction(model, inputdic):
    """
    Validate/configure substrate uptake on the split irreversible model before FBA.

    Order: task substrate ID -> base ID -> base+'_reverse'; if none allow uptake,
    add a dedicated import reaction and store its ID in inputdic['substrate_uptake'].
    """
    substrate_uptake_rate = float(inputdic.get('substrate_uptake_rate', 10))
    base_id, reverse_id = _normalize_substrate_exchange_ids(inputdic['substrate'])
    candidate_ids = []
    for rid in (inputdic['substrate'], base_id, reverse_id):
        if rid and rid not in candidate_ids:
            candidate_ids.append(rid)

    print(f"[substrate_uptake] Resolving uptake for task substrate={inputdic['substrate']!r}, "
          f"candidates={candidate_ids}, substrate_uptake_rate={substrate_uptake_rate}")

    for rxn_id in candidate_ids:
        if rxn_id not in model.reactions:
            print(f"[substrate_uptake]   {rxn_id}: not in model")
            continue
        reaction = model.reactions.get_by_id(rxn_id)
        direction = _reaction_uptake_flux_sign(reaction)
        allows = direction and _reaction_allows_uptake(reaction, direction)
        print(
            f"[substrate_uptake]   {rxn_id}: direction={direction}, "
            f"lb={reaction.lower_bound}, ub={reaction.upper_bound}, allows_uptake={allows}"
        )
        if allows:
            _configure_uptake_bounds(reaction, substrate_uptake_rate, direction)
            inputdic['substrate_uptake'] = rxn_id
            inputdic['substrate_uptake_direction'] = direction
            print(
                f"[substrate_uptake] Using {rxn_id} (direction={direction}, "
                f"lb={reaction.lower_bound}, ub={reaction.upper_bound})"
            )
            return rxn_id

    print("[substrate_uptake] No existing reaction allows uptake; adding import reaction")
    e_met = _guess_extracellular_metabolite(model, base_id)
    uptake_id, direction = _add_substrate_uptake_reaction(
        model, base_id, e_met, substrate_uptake_rate
    )
    inputdic['substrate_uptake'] = uptake_id
    inputdic['substrate_uptake_direction'] = direction
    reaction = model.reactions.get_by_id(uptake_id)
    print(
        f"[substrate_uptake] Using new {uptake_id} (direction={direction}, "
        f"lb={reaction.lower_bound}, ub={reaction.upper_bound})"
    )
    return uptake_id


def prepare_model(path_model,path_task,taskname):
    # model_file0="/home/sun/ETGEMS-10.20/data/iML1515_new.json"
    # model_file="/home/sun/ETGEMS-10.20/data/iML1515_irr_enz_constraint_adj_round2.json"
    try:
        # Load the input dictionary
        inputdic = load_input_dictionary(path_task, taskname)

        # Load the model
        model0 = load_model(path_model, inputdic['model'])

    except Exception as e:
        print(f"prepare model failed: {e}")
        
    # model_file0 = os.path.join(path_model, 'iCW773_uniprot_modification_del.json')
    # name_mapping_flie='./data/name_mapping.json'
    # model_file0 = os.path.join(path_model, 'iCW773_uniprot_modification_del.json')
    # name_mapping_flie='./data/name_mapping.json'


    E_total, Bvalue0, Kvalue0 = get_constraint_params(inputdic['model'])
    # Default: no thermo (SE mode) unless a model-specific g0 file is configured.
    g0file = ''
    metconcfile = 'metabolites_lnC.txt'
    kcatmwfile = ''

    if 'iCW773' in inputdic['model'] or 'iCW' in inputdic['model'] or 'ala_iCW' in inputdic['model']:
        g0file='iCW773_uniprot_modification_del_reaction_g0.csv'
        metconcfile='metabolites_lnC_cg1.txt'
        kcatmwfile='reaction_change_by_enzuse_PDH_n.csv'
    elif 'iML1515' in inputdic['model'] or 'iML' in inputdic['model']:
        g0file='iML1515_github_reaction_g0.csv'
        metconcfile='metabolites_lnC.txt'
        kcatmwfile='kcatmw_eciML1515.csv'
    elif 'ecBSU1' in inputdic['model'] or 'BSU' in inputdic['model']:
        g0file=''
        metconcfile='metabolites_lnC.txt'
        kcatmwfile='kcatmw_ecBSU1.csv'
    elif 'eciZM547' in inputdic['model'] or 'ZM' in inputdic['model']:
        g0file=''
        metconcfile='metabolites_lnC.txt'
        kcatmwfile='kcatmw_eciZM547.csv'
    elif 'ecMTM' in inputdic['model'] or 'MTM' in inputdic['model']:
        g0file=''
        metconcfile='metabolites_lnC.txt'
        kcatmwfile='kcatmw_ecMTM.csv'
    model_file = os.path.join(path_task, inputdic['model'])
    model=remove_unused_metabolites(model0)

    model_irrev_json = os.path.join(path_model, taskname+'_'+'model_irrev.json')
    model_irrev_xml = os.path.join(path_model, taskname+'_'+'model_irrev.xml')
    model_split_json = os.path.join(path_model, taskname+'_'+'model_irr_splited.json')
    model_split_xml = os.path.join(path_model, taskname+'_'+'model_irr_splited.xml')
    if os.path.exists(model_split_json):
        model = cobra.io.load_json_model(model_split_json)
        print(f"Loading cached isoenzyme-split model from: {model_split_json}")
    elif os.path.exists(model_split_xml):
        model = cobra.io.read_sbml_model(model_split_xml)
        print(f"Loading cached isoenzyme-split model from: {model_split_xml}")
    else:
        print("No cached isoenzyme-split model found; generating a new one")       
        if os.path.exists(model_irrev_json):
            print(f"Loading cached irreversible model from: {model_irrev_json}")
            model = cobra.io.load_json_model(model_irrev_json)
        elif os.path.exists(model_irrev_xml):
            print(f"Loading cached irreversible model from: {model_irrev_xml}")
            model = cobra.io.read_sbml_model(model_irrev_xml)
        else:
            print("No cached irreversible model found; generating a new one")
            model = convert_to_irreversible(model)
            cobra.io.save_json_model(model, model_irrev_json)
            print(f"Saved irreversible model cache to: {model_irrev_json}")
        model = isoenzyme_split(model)
        cobra.io.save_json_model(model, model_split_json)
        print(f"Saved isoenzyme-split model cache to: {model_split_json}")
    print(f"Model after isoenzyme split: {len(model.reactions)} reactions")
    resolve_substrate_uptake_reaction(model, inputdic)
    dictionary_model, model = trans_model2standard_json_etgem2(model,model_file)
    reaction_g0_file=os.path.join(path_model, g0file) if g0file else ''
    metabolites_lnC_file = os.path.join(path_model, metconcfile)
    reaction_kcat_MW_file=os.path.join(path_model, kcatmwfile)
    if inputdic['oxygenstate']=='aerobic':
        model.reactions.get_by_id(inputdic['O2']).upper_bound = 1000
    if inputdic['oxygenstate']=='micro_aerobic': 
        model.reactions.get_by_id(inputdic['O2']).upper_bound = 2 
    if inputdic['oxygenstate']=='anaerobic': 
        model.reactions.get_by_id(inputdic['O2']).upper_bound = 0

    # add data to dictionary_model
    #get kcat_dict and mw_dict
    kcat_mw=pd.read_csv(reaction_kcat_MW_file,index_col=0)
    kcat_dict={}
    mw_dict={}
    for i in kcat_mw.index:
        if i in model.reactions:
            kcat_dict[i]=kcat_mw.loc[i,'kcat']
            enz=str(model.reactions.get_by_id(i).gpr)
            flag1=enz.split(' and ')
            flag1.sort(key=None, reverse=False)
            flag2=' and '.join(flag1)
            mw_dict[flag2]=kcat_mw.loc[i,'MW']

    def add_parameter(dictionarymodel,kcatdata,kmdata,mwdata):
        for enz in dictionarymodel['enzyme']:
            if enz in mwdata.keys():
                dictionarymodel['enzyme'][enz]['MW']=mwdata[enz]
            for rea in dictionarymodel['enzyme'][enz]['reactions']:
                if  rea in kcatdata.keys():
                    dictionarymodel['enzyme'][enz]['reactions'][rea]['kcat']=kcatdata[rea]
                if  rea in kmdata.keys():
                    dictionarymodel['enzyme'][enz]['reactions'][rea]['km']=kmdata[rea]
        for enz in mwdata:
            if enz not in dictionarymodel['enzyme']:
                print(enz)
    add_parameter(dictionary_model,kcat_dict,{},mw_dict)

    # convert dictionary model to concretemodel_need_data

    Concretemodel_Need_Data=Get_Concretemodel_Need_Data2(model)
    get_dictionarymodel_data2(dictionary_model,Concretemodel_Need_Data,[])

    use_thermo_data = True
    if not g0file:
        print('No g0 file configured for this model, forcing SE mode.')
        use_thermo_data = False
    else:
        try:
            Get_Concretemodel_Need_Data_g0(Concretemodel_Need_Data, reaction_g0_file, metabolites_lnC_file, reaction_kcat_MW_file)
            rg0_df = Concretemodel_Need_Data.get('reaction_g0', pd.DataFrame()).copy()
            Inc = Concretemodel_Need_Data.get('metabolites_lnC', pd.DataFrame()).copy()

            if 'g0' in rg0_df.columns:
                rg0_df['g0'] = rg0_df['g0'].replace(0, np.nan)
                rg0_df.dropna(subset=['g0'], inplace=True)

                if 'AIRC3_reverse' in rg0_df.index:
                    rg0_df.at['AIRC3_reverse', 'g0'] = 0
                if 'ATPS4rpp_reverse_num2' in rg0_df.index:
                    rg0_df.at['ATPS4rpp_reverse_num2', 'g0'] = 0

            if not Inc.empty:
                for i in Concretemodel_Need_Data['metabolite_list']:
                    if i not in Inc.index:
                        Inc.loc[i, 'lnClb'] = -14.508658
                        Inc.loc[i, 'lnCub'] = -3.912023

            # If mapping fails or files are empty, disable thermo constraints.
            if rg0_df.empty or Inc.empty:
                use_thermo_data = False
            else:
                Concretemodel_Need_Data['reaction_g0'] = rg0_df
                Concretemodel_Need_Data['metabolites_lnC'] = Inc
        except Exception as e:
            print(f"Warning: thermodynamic data unavailable or unmapped, fallback to SE mode: {e}")
            use_thermo_data = False

    if not use_thermo_data:
        Concretemodel_Need_Data = disable_thermo_data(Concretemodel_Need_Data)
        Inc = Concretemodel_Need_Data['metabolites_lnC']
        print('Thermodynamic constraints disabled: reaction_g0 and metabolites_lnC are set empty, using SE mode.')
    else:
        Concretemodel_Need_Data['use_thermo_data'] = True
        print('Thermodynamic constraints enabled: using SET mode.')

    Concretemodel_Need_Data['E_total']=E_total
    Concretemodel_Need_Data['B_value']=Bvalue0
    Concretemodel_Need_Data['K_value']=Kvalue0
    # add oxygen state constraint
    if inputdic['oxygenstate']=='aerobic':
        model.reactions.get_by_id(inputdic['O2']).upper_bound = 1000
    if inputdic['oxygenstate']=='micro_aerobic': 
        model.reactions.get_by_id(inputdic['O2']).upper_bound = 2 
    if inputdic['oxygenstate']=='anaerobic': 
        model.reactions.get_by_id(inputdic['O2']).upper_bound = 0

    return Concretemodel_Need_Data,get_dictionarymodel_data2,Inc,model0,model,dictionary_model,inputdic 


# if inputdic['taskname']  ==  "optforce" :
def _log_calculate_biomass_step(step_label, eco_model, inputdic, Concretemodel_Need_Data=None, constr_coeff=None, extra=None):
    """Diagnostic logging only; does not change model state."""
    biomass_id = inputdic['biomass']
    substrate_rev = get_substrate_uptake_id(inputdic)
    lines = [f"[calculate_biomass] {step_label}"]
    try:
        lines.append(f"  objective_value={eco_model.obj()}")
    except Exception as exc:
        lines.append(f"  objective_value=ERROR ({exc})")
    try:
        lines.append(f"  biomass_flux[{biomass_id}]={value(eco_model.reaction[biomass_id])}")
    except Exception as exc:
        lines.append(f"  biomass_flux[{biomass_id}]=ERROR ({exc})")
    try:
        lines.append(f"  substrate_flux[{substrate_rev}]={value(eco_model.reaction[substrate_rev])}")
    except Exception as exc:
        lines.append(f"  substrate_flux[{substrate_rev}]=ERROR ({exc})")
    if hasattr(eco_model, 'e1') and Concretemodel_Need_Data is not None:
        try:
            mw_dict = Concretemodel_Need_Data.get('mw_dict', {})
            total_enzyme = sum(
                value(eco_model.e1[gene]) * mw_dict[gene]
                for gene in eco_model.e1
                if gene in mw_dict
            )
            lines.append(f"  total_enzyme_mass(sum e1*MW)={total_enzyme}")
        except Exception as exc:
            lines.append(f"  total_enzyme_mass=ERROR ({exc})")
    if constr_coeff is not None:
        lines.append(f"  constr_coeff={constr_coeff}")
    if extra:
        for key, val in extra.items():
            lines.append(f"  {key}={val}")
    print("\n".join(lines))


def calculate_biomass(Concretemodel_Need_Data,inputdic,model):
    print("[calculate_biomass] ===== start =====")
    print(
        f"[calculate_biomass] model={inputdic.get('model')}, "
        f"substrate={inputdic.get('substrate')}, biomass={inputdic.get('biomass')}, "
        f"substrate_uptake_rate={inputdic.get('substrate_uptake_rate')}, oxygen={inputdic.get('oxygenstate')}"
    )
    fba_mode = get_fba_mode(Concretemodel_Need_Data)
    print(
        f"[calculate_biomass] FBA mode={fba_mode}, "
        f"use_thermo_data={Concretemodel_Need_Data.get('use_thermo_data')}"
    )
    uptake_id = get_substrate_uptake_id(inputdic)
    try:
        uptake_rxn = model.reactions.get_by_id(uptake_id)
        print(
            f"[calculate_biomass] cobra uptake {uptake_id} bounds: "
            f"lb={uptake_rxn.lower_bound}, ub={uptake_rxn.upper_bound}, "
            f"reaction={uptake_rxn.reaction}"
        )
    except Exception as exc:
        print(f"[calculate_biomass] cobra uptake bounds: ERROR ({exc})")
    base_sub = inputdic['substrate']
    if base_sub in model.reactions and base_sub != uptake_id:
        base_rxn = model.reactions.get_by_id(base_sub)
        print(
            f"[calculate_biomass] cobra base {base_sub} bounds: "
            f"lb={base_rxn.lower_bound}, ub={base_rxn.upper_bound}"
        )

    print("[calculate_biomass] step 0: cobra reference (no enzyme/thermo constraints)")
    model.objective = inputdic['biomass']
    cobra_bio_sol = model.optimize()
    print(
        f"[calculate_biomass]   cobra max biomass: "
        f"obj={cobra_bio_sol.objective_value}, status={cobra_bio_sol.status}"
    )

    model.objective = uptake_id
    cobra_sub_sol = model.optimize()
    objvalue2 = cobra_sub_sol.objective_value
    print(
        f"[calculate_biomass]   cobra max uptake ({uptake_id}): "
        f"objvalue2={objvalue2}, status={cobra_sub_sol.status}"
    )

    constr_coeff={}
    constr_coeff['fix_reactions']={}
    # constr_coeff['fix_reactions']['EX_glc__D_e_reverse']=10
    constr_coeff['substrate_constrain'] =(get_substrate_uptake_id(inputdic),objvalue2)
    obj_name=inputdic['biomass']
    obj_target='maximize'
    e_total, b_value, k_value = get_constraint_params(inputdic['model'])
    Concretemodel_Need_Data['E_total']=e_total
    Concretemodel_Need_Data['B_value']=b_value
    Concretemodel_Need_Data['K_value']=k_value
    print(
        f"[calculate_biomass] constraint defaults: E_total={e_total}, "
        f"B_value={b_value}, K_value={k_value}"
    )

    print("[calculate_biomass] step 1: ET maximize biomass (substrate uptake from cobra)")
    EcoECM_FBA_protainmodel=FBA_template2(set_obj_value=True,obj_name=obj_name,obj_target=obj_target,mode=fba_mode,constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
    Model_Solve(EcoECM_FBA_protainmodel)
    print(EcoECM_FBA_protainmodel.obj())
    _log_calculate_biomass_step("step 1 result", EcoECM_FBA_protainmodel, inputdic, Concretemodel_Need_Data, constr_coeff=constr_coeff)

    biomass_flux_step1 = value(EcoECM_FBA_protainmodel.reaction[inputdic['biomass']])
    constr_coeff['fix_reactions'][inputdic['biomass']]=[biomass_flux_step1,np.inf]
    if 'iML' in inputdic['model']:
        constr_coeff['fix_reactions'][inputdic['biomass']]=[biomass_flux_step1*0.9,np.inf]
    print(
        f"[calculate_biomass] step 2 prep: fix biomass flux for B optimization: "
        f"{constr_coeff['fix_reactions'][inputdic['biomass']]}"
    )

    mode_now = get_fba_mode(Concretemodel_Need_Data)
    if mode_now == 'SET':
        print("[calculate_biomass] step 2: ET maximize thermodynamic B (biomass fixed)")
        EcoECM_FBA_protainmodel_B1=FBA_template2(set_obj_B_value=True,obj_name=obj_name,obj_target=obj_target,mode=mode_now,constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
        Model_Solve(EcoECM_FBA_protainmodel_B1)
        B_value1=EcoECM_FBA_protainmodel_B1.obj()
        print(EcoECM_FBA_protainmodel_B1.obj())
        _log_calculate_biomass_step("step 2 result", EcoECM_FBA_protainmodel_B1, inputdic, Concretemodel_Need_Data, constr_coeff=constr_coeff)
    else:
        B_value1 = Concretemodel_Need_Data.get('B_value', 0)
        print(f"SE mode: skip B optimization in calculate_biomass, use B_value={B_value1}")

    constr_coeff={}
    constr_coeff['fix_reactions']={}
    Concretemodel_Need_Data['B_value']=B_value1
    constr_coeff['substrate_constrain']=(get_substrate_uptake_id(inputdic),inputdic["substrate_uptake_rate"])
    print(
        f"[calculate_biomass] step 3 prep: B_value={B_value1}, "
        f"substrate_constrain={constr_coeff['substrate_constrain']}"
    )

    print("[calculate_biomass] step 3: ET maximize biomass -> v0_biomass")
    EcoECM_FBA_protainmodel=FBA_template2(set_obj_value=True,obj_name=obj_name,obj_target=obj_target,mode=get_fba_mode(Concretemodel_Need_Data),constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
    Model_Solve(EcoECM_FBA_protainmodel)
    v0_biomass=EcoECM_FBA_protainmodel.obj()
    print(EcoECM_FBA_protainmodel.obj())
    _log_calculate_biomass_step(
        "step 3 result (v0_biomass)",
        EcoECM_FBA_protainmodel,
        inputdic,
        Concretemodel_Need_Data,
        constr_coeff=constr_coeff,
        extra={'v0_biomass': v0_biomass},
    )

    constr_coeff['fix_E_total']=True
    print("[calculate_biomass] step 4: ET pFBA minimize sum of fluxes (fix_E_total=True)")
    EcoECM_PFBA_protainmodel_wild=FBA_template2(set_obj_V_value=True,mode=get_fba_mode(Concretemodel_Need_Data),constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
    Model_Solve(EcoECM_PFBA_protainmodel_wild)
    EcoECM_PFBA_protainmodel_wild.obj()
    _log_calculate_biomass_step("step 4 result (pFBA)", EcoECM_PFBA_protainmodel_wild, inputdic, Concretemodel_Need_Data, constr_coeff=constr_coeff)
    bio=showflux(EcoECM_PFBA_protainmodel_wild)
# mini enzyme
# mini enzyme
    constr_coeff={}
    constr_coeff['fix_reactions']={}
    constr_coeff['substrate_constrain']=(get_substrate_uptake_id(inputdic),inputdic["substrate_uptake_rate"])
    constr_coeff['biomass_constrain']=(inputdic['biomass'],v0_biomass*0.95)
    obj_name=inputdic['biomass']
    obj_target='minimize'
    print(
        f"[calculate_biomass] step 5: ET minimize total enzyme -> totalE "
        f"(biomass_constrain={constr_coeff['biomass_constrain']})"
    )
    EcoECM_FBA_protainmodel=FBA_template2(set_obj_sum_e=True,obj_name=obj_name,obj_target=obj_target,mode=get_fba_mode(Concretemodel_Need_Data),constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
    Model_Solve(EcoECM_FBA_protainmodel)
    totalE=EcoECM_FBA_protainmodel.obj()
    print(EcoECM_FBA_protainmodel.obj())
    _log_calculate_biomass_step(
        "step 5 result (totalE)",
        EcoECM_FBA_protainmodel,
        inputdic,
        Concretemodel_Need_Data,
        constr_coeff=constr_coeff,
        extra={'totalE': totalE},
    )
    print(
        f"[calculate_biomass] ===== done: B_value1={B_value1}, v0_biomass={v0_biomass}, "
        f"totalE={totalE}, objvalue2={objvalue2} ====="
    )
    return B_value1,v0_biomass,bio,totalE,objvalue2


def calculate_wildrange(Concretemodel_Need_Data,obj_name,inputdic): 
    results = {}
    constr_coeff = {}
    constr_coeff['fix_reactions'] = {} 
    Concretemodel_Need_Data['E_total'] = totalE * 1.01
    Concretemodel_Need_Data['B_value'] = B_value1 * 0.98
    constr_coeff['substrate_constrain'] = (get_substrate_uptake_id(inputdic),objvalue2)
    if 'iML' in inputdic['model']:
        Concretemodel_Need_Data['E_total']=totalE*1.01
        Concretemodel_Need_Data['B_value']=B_value1*0.99 
        constr_coeff['substrate_constrain'] = (get_substrate_uptake_id(inputdic), objvalue2)
    constr_coeff['biomass_constrain'] = (inputdic['biomass'], v0_biomass * 0.95) 
    obj_target = 'minimize'  
    EcoECM_FBA_protainmodel = FBA_template2(set_obj_value_e=True, obj_name=obj_name, obj_target=obj_target, mode=get_fba_mode(Concretemodel_Need_Data), constr_coeff=constr_coeff, Concretemodel_Need_Data=Concretemodel_Need_Data)
    try:
        Model_Solve(EcoECM_FBA_protainmodel) 
        min_value = EcoECM_FBA_protainmodel.obj()
    except Exception as e:
        print(f"Error occurred while solving: {e}")
        min_value = 0

    print(f"Objective: {obj_name}, Minimize: {min_value}")
    obj_target = 'maximize'  
    EcoECM_FBA_protainmodel = FBA_template2(set_obj_value_e=True, obj_name=obj_name, obj_target=obj_target,mode=get_fba_mode(Concretemodel_Need_Data), constr_coeff=constr_coeff, Concretemodel_Need_Data=Concretemodel_Need_Data)
    try:
        Model_Solve(EcoECM_FBA_protainmodel)  
        max_value = EcoECM_FBA_protainmodel.obj()
    except Exception as e:
        print(f"Error occurred while solving: {e}")
        max_value = 100

    print(f"Objective: {obj_name}, Maximize: {max_value}")

    results[obj_name] = {'range': [min_value, max_value]}


    return results


def calculate_product(Concretemodel_Need_Data,inputdic):
    constr_coeff={}
    constr_coeff['fix_reactions'] = {}
    # constr_coeff['fix_reactions']['EX_glc__D_e_reverse']=10
    constr_coeff['substrate_constrain'] = (get_substrate_uptake_id(inputdic),objvalue2)
    constr_coeff['biomass_constrain'] = (inputdic['biomass'],v0_biomass*0.1)
    obj_name = inputdic['product']
    obj_target = 'maximize'
    e_total, b_value, k_value = get_constraint_params(inputdic['model'])
    Concretemodel_Need_Data['E_total']=e_total
    Concretemodel_Need_Data['B_value']=b_value
    Concretemodel_Need_Data['K_value']=k_value
    EcoECM_FBA_protainmodel=FBA_template2(set_obj_value=True,obj_name=obj_name,obj_target=obj_target,mode=get_fba_mode(Concretemodel_Need_Data),constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
    Model_Solve(EcoECM_FBA_protainmodel)
    print(EcoECM_FBA_protainmodel.obj())

    constr_coeff['fix_reactions'][inputdic['product']]=[value(EcoECM_FBA_protainmodel.reaction[inputdic['product']])*0.9,np.inf]
    mode_now = get_fba_mode(Concretemodel_Need_Data)
    if mode_now == 'SET':
        EcoECM_FBA_protainmodel_B2=FBA_template2(set_obj_B_value=True,obj_name=obj_name,obj_target=obj_target,mode=mode_now,constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
        Model_Solve(EcoECM_FBA_protainmodel_B2)
        B_value2=EcoECM_FBA_protainmodel_B2.obj()
    else:
        B_value2 = Concretemodel_Need_Data.get('B_value', 0)
        EcoECM_FBA_protainmodel_B2 = EcoECM_FBA_protainmodel
        print(f"SE mode: skip B optimization in calculate_product, use B_value={B_value2}")
    print(EcoECM_FBA_protainmodel.obj())

    constr_coeff={}
    constr_coeff['fix_reactions']={}
    Concretemodel_Need_Data['B_value']=B_value2
    constr_coeff['biomass_constrain']=(inputdic['biomass'],v0_biomass*0.1)
    constr_coeff['substrate_constrain']=(get_substrate_uptake_id(inputdic),inputdic["substrate_uptake_rate"])
    EcoECM_FBA_protainmodel=FBA_template2(set_obj_value=True,obj_name=obj_name,obj_target=obj_target,mode=get_fba_mode(Concretemodel_Need_Data),constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
    Model_Solve(EcoECM_FBA_protainmodel)
    v1_product_max=EcoECM_FBA_protainmodel.obj()
    print(EcoECM_FBA_protainmodel.obj())
    pro=showflux(EcoECM_FBA_protainmodel)
# mini enzyme
    constr_coeff={}
    constr_coeff['fix_reactions']={}
    constr_coeff['substrate_constrain']=(get_substrate_uptake_id(inputdic),inputdic["substrate_uptake_rate"])
    constr_coeff['biomass_constrain']=(inputdic['biomass'],v0_biomass*0.1)
    constr_coeff['product_constrain']=(inputdic['product'],v1_product_max*0.95)
    obj_name=inputdic['product']
    obj_target='minimize'
    EcoECM_FBA_protainmodel=FBA_template2(set_obj_sum_e=True,obj_name=obj_name,obj_target=obj_target,mode=get_fba_mode(Concretemodel_Need_Data),constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
    Model_Solve(EcoECM_FBA_protainmodel)
    totalE2=EcoECM_FBA_protainmodel.obj()
    print(EcoECM_FBA_protainmodel.obj())
    return B_value2,v1_product_max,pro,totalE2,EcoECM_FBA_protainmodel_B2



def calculate_over(Concretemodel_Need_Data,obj_name,inputdic):
    results = {}
    constr_coeff={}
    constr_coeff['fix_reactions']={} 
    Concretemodel_Need_Data['E_total']=totalE2*1.01
    Concretemodel_Need_Data['B_value']=B_value2*0.9     
    constr_coeff['substrate_constrain']=(get_substrate_uptake_id(inputdic),objvalue2)
    constr_coeff['biomass_constrain']=(inputdic['biomass'],v0_biomass*0.1) 
    constr_coeff['product_constrain']=(inputdic['product'],v1_product_max*0.95) 
    obj_target = 'minimize'   

    EcoECM_FBA_protainmodel = FBA_template2(set_obj_value_e=True, obj_name=obj_name, obj_target=obj_target, mode=get_fba_mode(Concretemodel_Need_Data), constr_coeff=constr_coeff, Concretemodel_Need_Data=Concretemodel_Need_Data)
    try:
        Model_Solve(EcoECM_FBA_protainmodel) 
        min_value = EcoECM_FBA_protainmodel.obj()
    except Exception as e:
        print(f"Error occurred while solving: {e}")
        min_value = 0

    print(f"Objective: {obj_name}, Minimize: {min_value}")

    obj_target = 'maximize'  

    EcoECM_FBA_protainmodel = FBA_template2(set_obj_value_e=True, obj_name=obj_name, obj_target=obj_target,mode=get_fba_mode(Concretemodel_Need_Data), constr_coeff=constr_coeff, Concretemodel_Need_Data=Concretemodel_Need_Data)
    try:
        Model_Solve(EcoECM_FBA_protainmodel)  
        max_value = EcoECM_FBA_protainmodel.obj()
    except Exception as e:
        print(f"Error occurred while solving: {e}")
        max_value = 100

    print(f"Objective: {obj_name}, Maximize: {max_value}")

    results[obj_name] = {'range': [min_value, max_value]}
    return results

def read_file(path_results):
    wild_file_path = os.path.join(path_results, 'wild-enzyme.json')
    over_file_path = os.path.join(path_results, 'over-enzyme.json') 
    with open(wild_file_path, 'r') as enzyme_results_file:
        enzyme_results_data = json.load(enzyme_results_file)
    with open(over_file_path, 'r') as enzyme_overresults2_file:
        enzyme_overresults2_data = json.load(enzyme_overresults2_file)
    return enzyme_results_data,enzyme_overresults2_data


def compare_results(enzyme_results_data, enzyme_overresults_data): 
    # ko_flux
    mean_values = []
    for values in enzyme_results_data.values():
        if values.get('range') and len(values['range']) >= 2:
            range_values = values['range']
            mean_value = sum(range_values) / len(range_values)
            mean_values.append(mean_value)
    mean = np.mean(mean_values)
    std_value = np.std(mean_values)
    R = mean - 3 * std_value  

    ko_data = []

    for enzyme, values in enzyme_results_data.items():
        if values.get('range') and len(values['range']) >= 2:
            range_values = values['range']
            first_value = range_values[0]


            if enzyme in enzyme_overresults_data:
                overresults_values = enzyme_overresults_data[enzyme]
                if overresults_values.get('range') and len(overresults_values['range']) >= 2:
                    over_second_value = overresults_values['range'][1]


                    if first_value >= R and over_second_value <= R and range_values != [0, 0]:
                        ko_data.append((enzyme, range_values))

    # up flux
    up_data = []

    for enzyme, values in enzyme_results_data.items():
        if values.get('range') and len(values['range']) >= 2:
            range_values = values['range']
            second_value = range_values[1]
            if enzyme in enzyme_overresults_data:
                overresults_values = enzyme_overresults_data[enzyme]
                if overresults_values.get('range') and len(overresults_values['range']) >= 2:
                    over_first_value = overresults_values['range'][0]
                    over_second_value = overresults_values['range'][1]


                    if second_value <= over_first_value and range_values != [0, 0]:
                        up_data.append((enzyme, range_values))

    # down_flux
    down_data = []

    for enzyme, values in enzyme_results_data.items():
        if values.get('range') and len(values['range']) >= 2:
            range_values = values['range']
            first_value = range_values[0]

            if enzyme in enzyme_overresults_data:
                overresults_values = enzyme_overresults_data[enzyme]
                if overresults_values.get('range') and len(overresults_values['range']) >= 2:
                    over_first_value = overresults_values['range'][0]
                    over_second_value = overresults_values['range'][1]

                    if first_value >= over_second_value and first_value >= over_first_value and range_values != [0, 0]:
                        down_data.append((enzyme, range_values))


    range_change = []
    mean_change = []
    for enzyme, values in enzyme_results_data.items():
        if values.get('range') and len(values['range']) >= 2:
            range_values = values['range']
            first_value = range_values[0]
            second_value = range_values[1]

            if enzyme in enzyme_overresults_data:
                overresults_values = enzyme_overresults_data[enzyme]
                if overresults_values.get('range') and len(overresults_values['range']) >= 2:
                    over_first_value = overresults_values['range'][0]
                    over_second_value = overresults_values['range'][1]

                    over_range = (over_first_value + over_second_value) / 2
                    wild_range = (first_value + second_value) / 2
                    
                    # Ìí¼Ó¼ì²éÌõ¼þ£¬È·±£ wild_range ²»ÎªÁã
                    if wild_range != 0:
                        range_value = over_range / wild_range
                        mean_value = (over_range + wild_range) / 2
                        range_change.append((enzyme,range_value))
                        mean_change.append((enzyme,mean_value))


    return ko_data,up_data,down_data,range_change,mean_change


def gene_reaction_map1(reaction_list,model):
    equation_dict = {}

    for reaction_id in reaction_list:
        equation = model.reactions.get_by_id(reaction_id).reaction
        equation_dict[reaction_id] = equation
    return equation_dict

def ref_e_con(inputdic,Concretemodel_Need_Data):
    # wild fix min_enz, biomass and B, maximize product, get value(model.reactions[i]) for i in reaction_list
    constr_coeff={}
    constr_coeff['fix_reactions']={}
    # constr_coeff['fix_reactions']['EX_glc__D_e_reverse']=10
    e_total, _, _ = get_constraint_params(inputdic['model'])
    Concretemodel_Need_Data['E_total']=e_total
    constr_coeff['substrate_constrain']=(get_substrate_uptake_id(inputdic),inputdic['substrate_uptake_rate'])
    obj_name=inputdic['biomass']
    obj_target='maximize'

    EcoECM_FBA_protainmodel_MIN=FBA_template2(set_obj_value=True,obj_name=obj_name,obj_target=obj_target,mode=get_fba_mode(Concretemodel_Need_Data),constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
    Model_Solve(EcoECM_FBA_protainmodel_MIN)
    print(EcoECM_FBA_protainmodel_MIN.obj())
    # wild e1
    e_ref={x: value(EcoECM_FBA_protainmodel_MIN.e1[x]) for x in EcoECM_FBA_protainmodel_MIN.e1}
    mw_dict=Concretemodel_Need_Data['mw_dict']
    reaction_dict = {x: value(EcoECM_FBA_protainmodel_MIN.reaction[x]) for x in EcoECM_FBA_protainmodel_MIN.reaction}
    reaction_dict_bio = {key: value for key, value in reaction_dict.items() if value > 1e-1}
    E_refdict = {}
    kcat_dict=Concretemodel_Need_Data['kcat_dict']
    kcat_dict={key:value for key,value in kcat_dict.items()if key in reaction_dict}
    for gene in e_ref:
        if gene in mw_dict:
            result = e_ref[gene] * mw_dict[gene]
            E_refdict[gene] = result

    E_refdict = {k: v for k, v in E_refdict.items() if v > 0}
    return E_refdict,mw_dict,reaction_dict_bio,kcat_dict

def gene_reaction_map(enzyme_list, filtered_reaction_dict,kcat_dict,Concretemodel_Need_Data,model):
    enzyme_reaction = Concretemodel_Need_Data['enzyme_rxns_dict']
    gene_reaction_mapping = {}
    gene_kcat_mapping = {}
    new_gene_reaction_mapping = {}

    for gene, reactions in enzyme_reaction.items():
        reaction_equations = []
        for reaction_name in reactions:
            if reaction_name in filtered_reaction_dict:
                reaction = model.reactions.get_by_id(reaction_name)
                reaction_equation = f'{reaction_name} ({filtered_reaction_dict[reaction_name]}): {reaction.reaction}'
                reaction_equations.append(reaction_equation)
        if reaction_equations:
            gene_reaction_mapping[gene] = ", ".join(reaction_equations)

    for gene in enzyme_list:
        if gene in gene_reaction_mapping:
            new_gene_reaction_mapping[gene] = gene_reaction_mapping[gene]

    for gene, reactions in enzyme_reaction.items():
        for reaction_name in reactions:
            if reaction_name in kcat_dict:
                kcat_value = kcat_dict[reaction_name]
                if gene in gene_kcat_mapping:
                    gene_kcat_mapping[gene].append((reaction_name, kcat_value))
                else:
                    gene_kcat_mapping[gene] = [(reaction_name, kcat_value)]
    

    return new_gene_reaction_mapping, gene_reaction_mapping,gene_kcat_mapping
# get reaction flux
def reaction_flux(inputdic,Concretemodel_Need_Data):
    # Çó×î´ó²úÆ·ËÙÂÊ
    constr_coeff={}
    constr_coeff['fix_reactions']={}
    # constr_coeff['fix_reactions']['EX_glc__D_e_reverse']=10
    Concretemodel_Need_Data['E_total']=totalE2
    if 'iML' in inputdic['model']:   
        Concretemodel_Need_Data['E_total']=totalE*1.01
    constr_coeff['substrate_constrain']=(get_substrate_uptake_id(inputdic),inputdic['substrate_uptake_rate'])
    constr_coeff['biomass_constrain']=(inputdic['biomass'],v0_biomass*0.1)
    obj_name=inputdic['product']
    obj_target='maximize'

    EcoECM_FBA_protainmodel_pro_max=FBA_template2(set_obj_value=True,obj_name=obj_name,obj_target=obj_target,mode=get_fba_mode(Concretemodel_Need_Data),constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
    Model_Solve(EcoECM_FBA_protainmodel_pro_max)
    print(EcoECM_FBA_protainmodel_pro_max.obj())

    # get reaction flux
    reaction_dict = {x: value(EcoECM_FBA_protainmodel_pro_max.reaction[x]) for x in EcoECM_FBA_protainmodel_pro_max.reaction}
    reaction_dict = {key: value for key, value in reaction_dict.items() if value > 1e-1}
    kcat_dict=Concretemodel_Need_Data['kcat_dict']
    kcat_dict={key:value for key,value in kcat_dict.items()if key in reaction_dict}
    return reaction_dict,kcat_dict,EcoECM_FBA_protainmodel_pro_max


# analysis enzyme usage_over
# e1*MW
def enzyme_usage_over(EcoECM_FBA_protainmodel_pro_max,Concretemodel_Need_Data):
    e_dict = {x: value(EcoECM_FBA_protainmodel_pro_max.e1[x]) for x in EcoECM_FBA_protainmodel_pro_max.e1}
    mw_dict=Concretemodel_Need_Data['mw_dict']
    E_dict = {}
    for gene in e_dict:
        if gene in mw_dict:
            result = e_dict[gene] * mw_dict[gene]
            E_dict[gene] = result
    # sum e
    total_sum = sum(E_dict.values())
    print("sum:", total_sum)
    # next (single enzyme usage_over)
    normalized_E_dict = {gene: (value / total_sum) * 100 for gene, value in E_dict.items()}
    return normalized_E_dict,E_dict,e_dict
def fold_change(E_dict,E_refdict):
    Fold_change={}
    for gene in E_dict:
        if gene in E_refdict:
            results = E_dict[gene]/E_refdict[gene]
            Fold_change[gene]=results
    return  Fold_change    

def must_df(enzyme_results_data,enzyme_overresults2_data,range_change,mean_change,ko_data,up_data,down_data):
    # make excel
    wild_data =[{'gene':gene,'enzyme_wild':[format(value,'.3e') for value in data['range']]}for gene,data in enzyme_results_data.items()]
    df1=pd.DataFrame(wild_data)
    over_data =[{'gene':gene,'enzyme_over':[format(value,'.3e') for value in data['range']]}for gene,data in enzyme_overresults2_data.items()]
    df2=pd.DataFrame(over_data)
    meged_df =pd.merge(df1,df2,on='gene',how='inner')
    meged_df['reaction_bio']=meged_df['gene'].map(new_gene_reaction_mapping_bio)
    meged_df['reaction_pro']=meged_df['gene'].map(new_gene_reaction_mapping)
    # enzyme_usage_dict = single_enzyme_usage_over.to_dict()
    meged_df['ref_e_con(g/gDW)']=meged_df['gene'].map(E_refdict).apply(lambda x: format(x,'.3e'))
    meged_df['over_e_con(g/gDW)']=meged_df['gene'].map(E_dict).apply(lambda x: format(x,'.3e'))
    range_change_data = [{'gene': enzyme, 'Fold_change': value} for enzyme,value in range_change]
    mean_change_data = [{'gene': enzyme, 'mean_change': value} for enzyme,value in mean_change]
    # 将 range_change_data 的 'range_change' 值处理为空或为 'full' 的情况
    for data in range_change_data:
        if data['Fold_change'] in [None, 'full']:
            data['Fold_change'] = 0
        elif isinstance(data['Fold_change'], str) and not data['Fold_change'].replace('.', '', 1).isdigit():
            # 如果值是其他非数值类型字符串，直接设置为 0
            data['Fold_change'] = 0
    for data in range_change_data:
        value = data['Fold_change']
        try:
            value = float(value)
        except (TypeError, ValueError):
            value = 0.0
        if value > 0:
            data['Fold_change'] = round(float(np.log2(value)), 3)
        else:
            data['Fold_change'] = 0.0
    range_change_df = pd.DataFrame(range_change_data)
    mean_change_df = pd.DataFrame(mean_change_data)
    meged_df = pd.merge(meged_df, range_change_df, on='gene', how='inner')
    meged_df = pd.merge(meged_df, mean_change_df, on='gene', how='inner')
    meged_df['Fold_change'] = meged_df['Fold_change'].apply(lambda x: round(x, 3))
    meged_df['mean_change'] = meged_df['mean_change'].apply(lambda x: round(x, 3))
    meged_df['enzyme_usage_over']=meged_df['gene'].map(normalized_E_dict)
    meged_df['enzyme_usage_over'] = meged_df['enzyme_usage_over'].apply(float)
    meged_df = meged_df.sort_values(by='enzyme_usage_over', ascending=False)
    meged_df['enzyme_usage_over'] = meged_df['enzyme_usage_over'].apply(lambda x: format(x, '.3e'))
    meged_df['enzyme_usage_over'] = meged_df['enzyme_usage_over'].apply(lambda x: float(x) if x != 0 else 0.0)
    meged_df.reset_index(drop=True, inplace=True)

    meged_df['manipulations'] = None

    for gene, _ in ko_data:
        meged_df.loc[meged_df['gene'] == gene, 'manipulations'] = 'ko'
    for gene, _ in up_data:
        meged_df.loc[meged_df['gene'] == gene, 'manipulations'] = 'Up'
    for gene, _ in down_data:
        meged_df.loc[meged_df['gene'] == gene, 'manipulations'] = 'down'
    # Fold_Change threshold rule on log2 scale
    meged_df.loc[(meged_df['manipulations'] == 'Up') & (meged_df['Fold_change'] <= 1), 'manipulations'] = None
    meged_df.loc[(meged_df['manipulations'] == 'down') & (meged_df['Fold_change'] >= -1), 'manipulations'] = None
    meged_df['e1(mmol/gDW)']=meged_df['gene'].map(e_dict).apply(lambda x: format(x,'.3e'))
    meged_df['kcat(1/h)'] = meged_df['gene'].map(gene_kcat_mapping)
    mw_dict=Concretemodel_Need_Data['mw_dict']
    meged_df['mw(g/mg)'] = meged_df['gene'].map(mw_dict).apply(lambda x: format(x,'.3e'))
    # meged_df['manipulations'] = meged_df.apply(lambda row: 'None' if pd.isnull(row['reaction_bio']) and pd.isnull(row['reaction_pro']) else row['manipulations'], axis=1)
    return meged_df

def get_reactions_from_expression(expression, model):
    """
    根据基因表达式获取相关反应的 ID 列表。

    参数:
    - expression: 基因逻辑表达式（如 'Cgl1271 and Cgl1272 or Cgl1273'）
    - model: COBRA 模型对象

    返回:
    - reaction_ids: 包含所有相关反应 ID 的列表
    """
    reaction_ids = set()  # 用集合来存储反应 ID，避免重复
    expression = str(expression).strip()

    if not expression:
        return []

    # 通用提取：支持 Cgl/BSU/b 等不同命名风格
    tokens = re.findall(r"[A-Za-z0-9_\-\.]+", expression)
    keywords = {"and", "or", "not", "AND", "OR", "NOT", "True", "False"}

    gene_ids = []
    for token in tokens:
        if token in keywords:
            continue
        try:
            model.genes.get_by_id(token)
            gene_ids.append(token)
        except KeyError:
            pass

    # 兜底：表达式本身就是单个基因 ID
    if not gene_ids:
        try:
            model.genes.get_by_id(expression)
            gene_ids = [expression]
        except KeyError:
            gene_ids = []

    for gene_id in set(gene_ids):
        try:
            gene = model.genes.get_by_id(gene_id)
            for reaction in gene.reactions:
                reaction_ids.add(reaction.id)
        except KeyError:
            pass

    return list(reaction_ids)
def basic(bio,pro):
    data = {
    'Metabolite': list(bio.keys()),
    'Value': list(bio.values())
    }
    bio_df = pd.DataFrame(data)
    data = {
    'Metabolite': list(pro.keys()),
    'Value': list(pro.values())
    }
    pro_df = pd.DataFrame(data)
    keyInfo=dict()
    keyInfo['wildtype']=dict()
    keyInfo['over']=dict()
    keyInfo['wildtype']['mdf']=B_value1
    keyInfo['wildtype']['growth']=v0_biomass
    keyInfo['wildtype']['min_enz']=totalE
    keyInfo['over']['mdf']=B_value2
    keyInfo['over']['product']=v1_product_max
    keyInfo['over']['min_enz']=totalE2
    KeyInfo = pd.DataFrame(keyInfo)

    return bio_df,pro_df,KeyInfo

def bottleneck_reactions(EcoECM_FBA_protainmodel_B2,model):
    if not hasattr(EcoECM_FBA_protainmodel_B2, 'Df') or not hasattr(EcoECM_FBA_protainmodel_B2, 'B'):
        print('SE mode: no Df/B variables available, skip bottleneck reaction analysis.')
        empty_df = pd.DataFrame(columns=["Reaction_ID", "Equation", "df_Value", "flux_value"])
        return [], empty_df
    Df_dict = {x: value(EcoECM_FBA_protainmodel_B2.Df[x]) for x in EcoECM_FBA_protainmodel_B2.Df}
    # get reaction and flux
    reaction={x:value(EcoECM_FBA_protainmodel_B2.reaction[x]) for x in EcoECM_FBA_protainmodel_B2.reaction}
    # get reaction>10-3
    filtered_reaction = {k: v for k, v in reaction.items() if v > 1e-3}
    # get df=B
    B_dict={x:value(EcoECM_FBA_protainmodel_B2.B[x]) for x in EcoECM_FBA_protainmodel_B2.B}
    B = next(iter(B_dict.values()))
    matching_keys = [key for key, value in Df_dict.items() if value == B]
    # find bottleneck_rxn
    matching_keys_set = set(matching_keys)
    filtered_reaction_set = set(filtered_reaction.keys())
    same_reactions = matching_keys_set.intersection(filtered_reaction_set)
    bottleneck_reactions_list = list(same_reactions)
    reaction_data = []
    for reaction_id in bottleneck_reactions_list:
        equation = model.reactions.get_by_id(reaction_id).reaction
        df_value = Df_dict.get(reaction_id, None)
        flux_value =reaction.get(reaction_id,None)
        reaction_data.append((reaction_id, equation, df_value,flux_value))

    Df = pd.DataFrame(reaction_data, columns=["Reaction_ID", "Equation", "df_Value","flux_value"])
    Df = Df.sort_values(by='flux_value', ascending=False)
    Df['df_Value'] = Df['df_Value'].apply(lambda x: format(x, '.3e') if isinstance(x, (int, float)) else x)
    Df['flux_value']=Df['flux_value'].apply(lambda x: format(x, '.3e') if isinstance(x, (int, float)) else x)
    return bottleneck_reactions_list,Df


def output(path_results,bio_df,pro_df,meged_df):
    filename=os.path.join(path_results, 'results.xlsx')
    with pd.ExcelWriter(filename) as writer:
        KeyInfo.to_excel(writer,sheet_name='basic',index=True)
        bio_df.to_excel(writer,sheet_name='EX_ref',index=True)
        pro_df.to_excel(writer,sheet_name='EX_high',index=True)
        meged_df.to_excel(writer,sheet_name='MUST',index=True)
        Df.to_excel(writer,sheet_name='bottleneck-Rxn',index=True)    



def detail_put_optforce(model_input,meged_df):
    reactions = {}
    reaction_ids = [] 
    up_df = meged_df.loc[meged_df['manipulations'] == 'Up']
    up_df = up_df.sort_values(by='enzyme_usage_over', ascending=False)
    up_reaction=up_df['gene'].tolist()
    for i in up_reaction:
        try:    
            up_reaction_name = model_input.genes.get_by_id(i).annotation.get('uniprot')
            enzyme_usage_values = up_df.loc[up_df['gene'] == i, 'enzyme_usage_over'].tolist()
            enzyme_usage_values = [float(value) for value in enzyme_usage_values if pd.notna(value)]
            if not enzyme_usage_values:
                enzyme_usage_values = [0.0]
            range_change_value = up_df.loc[up_df['gene'] == i, 'Fold_change'].values.tolist()
            range_change_value = range_change_value[0] if range_change_value else None
            reaction_ids = []  # 初始化一个空列表
            reaction_ids = get_reactions_from_expression(i, model0)

            reaction_dict = {
                    "Modification": "UP",
                    "Fold_Change":range_change_value ,
                    "Reaction_ID": reaction_ids,
                    "Uniprot_ID": up_reaction_name,
                    "Gene_ID": i,
                    "Enzyme_Cost": enzyme_usage_values,
                    "figpath":os.path.join(path_results4,'target_fig',i+'.png')
                                                    }
            reactions[i] = reaction_dict
        except KeyError:
            pass
    
    down_df = meged_df.loc[meged_df['manipulations'] == 'down']
    down_df = down_df.sort_values(by='enzyme_usage_over', ascending=False)
    down_reaction=down_df['gene'].tolist()
    for i in down_reaction:
            try:
                    down_reaction_name = model_input.genes.get_by_id(i).annotation.get('uniprot')
                    enzyme_usage_values = down_df.loc[down_df['gene'] == i, 'enzyme_usage_over'].tolist()
                    enzyme_usage_values = [float(value) for value in enzyme_usage_values if pd.notna(value)]
                    if not enzyme_usage_values:
                        enzyme_usage_values = [0.0]
                    range_change_value = down_df.loc[down_df['gene'] == i, 'Fold_change'].values.tolist()
                    range_change_value = range_change_value[0] if range_change_value else None
                    reaction_ids = []  # 初始化一个空列表
                    reaction_ids = get_reactions_from_expression(i, model0)
                    reaction_dict = {
                        "Modification": "Down",
                        "Fold_Change":range_change_value ,
                        "Reaction_ID": reaction_ids,
                        "Uniprot_ID": down_reaction_name,
                        "Gene_ID": i,
                        "Enzyme_Cost": enzyme_usage_values,
                        "figpath":os.path.join(path_results4,'target_fig',i+'.png')
                    }
                    reactions[i] = reaction_dict
            except KeyError:
                    pass
    return reactions,up_reaction,down_reaction

def output_web_optforce(reactions):
    output={}
    output['summary'] = {}
    output['yield'] = {}
    up_number = sum(1 for item in reactions.values() if str(item.get('Modification', '')).upper() == 'UP')
    down_number = sum(1 for item in reactions.values() if str(item.get('Modification', '')).upper() == 'DOWN')
    summary = {
        'UP':up_number,
        'Down':down_number,
        'KO':0,
        'Prod_Rate(mmol/gDW/h)':round(v1_product_max,3),
        'Growth(1/h)':round(v0_biomass,3),
        "Specific_Growth_Rate(1/h)":round(v0_biomass,3)
    }
    Yield = {
        'Ycm':Ycm,
        'Ycg':Ycm,
        'Yrm':round(degree_bio/degree_pro,3),
        'Yrg':degree_bio*substrate_molecular/degree_pro*product_molecular,
        'Ysm':Yield0,
        'Ysg':Mass_Yield0,
        'Yscm':Carbon_Yield0,
        'Yscg':Carbon_Yield0
    }
    output['summary'] = summary
    output['yield'] = Yield
    output['genes']=reactions
    output_file = os.path.join(path_results4, 'output.json')
    with open(output_file, 'w') as json_file:
        json.dump(output, json_file) 
    return output


def up_picture_optforce(meged_df):
    up_df = meged_df.loc[meged_df['manipulations'] == 'Up']
    first_10_columns = up_df.iloc[:, :4]  # 提取前 4 列
    result_dict = {}  # 初始化结果字典

    for index, row in first_10_columns.iterrows():
        reaction_name = row['gene']  # 提取反应名称
        flux_wild_values = row['enzyme_wild']  # 提取野生型酶通量
        flux_over_values = row['enzyme_over']  # 提取过表达型酶通量
        
        # 将数据存储到临时字典中
        temp_dict = {
            "Gene":reaction_name,
            "wild_range": flux_wild_values,
            "eng_range": flux_over_values,
        }
        # 将临时字典存储到结果字典中，键为 reaction_name
        result_dict[reaction_name] = temp_dict

    # 将字典转换为 JSON 格式的字符串
    result_json = json.dumps(result_dict, indent=2)

    # 输出 JSON 字符串
    up_file = os.path.join(path_results4, 'up.json')
    with open(up_file, 'w') as json_file:
        json_file.write(result_json)

    return result_json

def down_picture_optforce(meged_df):
    down_df = meged_df.loc[meged_df['manipulations'] == 'down']
    first_10_columns = down_df.iloc[:, :4]
    result_dict_down = {}
    for index, row in first_10_columns.iterrows():
        reaction_name = row['gene']  # 提取反应名称
        flux_wild_values = row['enzyme_wild']  # 提取野生型酶通量
        flux_over_values = row['enzyme_over']  # 提取过表达型酶通量
        
        # 将数据存储到临时字典中
        temp_dict = {
            "Gene":reaction_name,
            "wild_range": flux_wild_values,
            "eng_range": flux_over_values,
        }
        # 将临时字典存储到结果字典中，键为 reaction_name
        result_dict_down[reaction_name] = temp_dict

    # 将字典转换为 JSON 格式的字符串
    result_json = json.dumps(result_dict_down, indent=2)

    # 输出 JSON 字符串
    down_file = os.path.join(path_results4, 'down.json')
    with open(down_file, 'w') as json_file:
        json_file.write(result_json)

    return result_json


def flux(model0):
    # Step 1: 提取 up_reaction 中的基因和反应关系
    all_reactions = set()
    reactions_with_ranges = {}

    for gene_expression in up_reaction:  # 直接遍历列表
        related_reactions = get_reactions_from_expression(gene_expression, model0)
        all_reactions.update(related_reactions)
    for gene_expression in down_reaction:  # 直接遍历列表
        related_reactions = get_reactions_from_expression(gene_expression, model0)
        all_reactions.update(related_reactions)
    # Step 2: 将 meged_df 的基因和 Fold_change 对应到反应
    for index, row in meged_df.iterrows():
        gene = row['gene']
        fold_change = row['Fold_change']
        if gene in up_reaction:
            related_reactions = get_reactions_from_expression(gene, model0)
            for reaction in related_reactions:
                reactions_with_ranges[reaction] = fold_change
    for index, row in meged_df.iterrows():
        gene = row['gene']
        fold_change = row['Fold_change']
        if gene in down_reaction:
            related_reactions = get_reactions_from_expression(gene, model0)
            for reaction in related_reactions:
                reactions_with_ranges[reaction] = fold_change
    # Step 3: 生成反应列表并初始化通量
    reaction_list = [rea.id for rea in model0.reactions]
    fluxnew = {}

    for i in reaction_list:
        fluxnew[i] = 0  # 初始化通量为 0
        for j in reactions_with_ranges.keys():
            if i == j:
                fluxnew[i] += reactions_with_ranges[j]
            if i + '_num' == j[:-1]:
                fluxnew[i] += reactions_with_ranges[j]
            if i + '_reverse' == j or i + '_reverse_num' == j[:-1]:
                fluxnew[i] -= reactions_with_ranges[j]
    flux_file = os.path.join(path_results4, 'flux.json')
    with open(flux_file, 'w') as json_file:
        json.dump(fluxnew, json_file) 
    return fluxnew      


def calculate_product_fseof(Concretemodel_Need_Data,inputdic):
    model.objective = get_substrate_uptake_id(inputdic)
    objvalue2=model.optimize().objective_value
    constr_coeff={}
    constr_coeff['fix_reactions']={}
    # constr_coeff['fix_reactions']['EX_glc__D_e_reverse']=10
    constr_coeff['substrate_constrain']=(get_substrate_uptake_id(inputdic),objvalue2)
    obj_name=inputdic['product']
    obj_target='maximize'
    e_total, b_value, k_value = get_constraint_params(inputdic['model'])
    Concretemodel_Need_Data['E_total']=e_total
    Concretemodel_Need_Data['B_value']=b_value
    Concretemodel_Need_Data['K_value']=k_value
    EcoECM_FBA_protainmodel_pro=FBA_template2(set_obj_value=True,obj_name=obj_name,obj_target=obj_target,mode=get_fba_mode(Concretemodel_Need_Data),constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
    Model_Solve(EcoECM_FBA_protainmodel_pro)
    print(EcoECM_FBA_protainmodel_pro.obj())
    product=EcoECM_FBA_protainmodel_pro.obj()
    return product,objvalue2,EcoECM_FBA_protainmodel_pro


def calculate_product_fseof_bio(Concretemodel_Need_Data,inputdic):
    model.objective = get_substrate_uptake_id(inputdic)
    objvalue2=model.optimize().objective_value
    constr_coeff={}
    constr_coeff['fix_reactions']={}
    # constr_coeff['fix_reactions']['EX_glc__D_e_reverse']=10
    constr_coeff['substrate_constrain']=(get_substrate_uptake_id(inputdic),objvalue2)
    obj_name=inputdic['biomass']
    obj_target='maximize'
    e_total, b_value, k_value = get_constraint_params(inputdic['model'])
    Concretemodel_Need_Data['E_total']=e_total
    Concretemodel_Need_Data['B_value']=b_value
    Concretemodel_Need_Data['K_value']=k_value
    EcoECM_FBA_protainmodel_pro=FBA_template2(set_obj_value=True,obj_name=obj_name,obj_target=obj_target,mode=get_fba_mode(Concretemodel_Need_Data),constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
    Model_Solve(EcoECM_FBA_protainmodel_pro)
    print(EcoECM_FBA_protainmodel_pro.obj())
    biomass=EcoECM_FBA_protainmodel_pro.obj()
    return biomass,objvalue2,EcoECM_FBA_protainmodel_pro

def gene_reaction_map2(enzyme_list, filtered_reaction_dict,kcat_dict,Concretemodel_Need_Data,model):
    enzyme_reaction = Concretemodel_Need_Data['enzyme_rxns_dict']
    gene_reaction_mapping = {}
    gene_kcat_mapping = {}
    new_gene_reaction_mapping = {}

    for gene, reactions in enzyme_reaction.items():
        reaction_equations = []
        for reaction_name in reactions:
            if reaction_name in filtered_reaction_dict:
                reaction = model.reactions.get_by_id(reaction_name)
                reaction_equation = f'{reaction_name} ({filtered_reaction_dict[reaction_name]}): {reaction.reaction}'
                reaction_equations.append(reaction_equation)
        if reaction_equations:
            gene_reaction_mapping[gene] = ", ".join(reaction_equations)

    for gene in enzyme_list:
        if gene in gene_reaction_mapping:
            new_gene_reaction_mapping[gene] = gene_reaction_mapping[gene]

    for gene, reactions in enzyme_reaction.items():
        for reaction_name in reactions:
            if reaction_name in kcat_dict:
                kcat_value = kcat_dict[reaction_name]
                if gene in gene_kcat_mapping:
                    gene_kcat_mapping[gene].append((reaction_name, kcat_value))
                else:
                    gene_kcat_mapping[gene] = [(reaction_name, kcat_value)]
    return new_gene_reaction_mapping, gene_reaction_mapping,gene_kcat_mapping


def biomass(product,inputdic):
    exlist = list(np.linspace(0,product,10))
    exlistn = []
    for i in exlist:
        i = format(i,'.2f')
        exlistn.append(float(i))
    FSEOFdf = pd.DataFrame()
    reactiondf = pd.DataFrame()
    obj_name = inputdic['biomass']
    e_total, b_value, k_value = get_constraint_params(inputdic['model'])
    for i in exlistn:
        cond = i
        obj_name = inputdic['biomass']
        obj_target='maximize'
        mode = get_fba_mode(Concretemodel_Need_Data)
        constr_coeff={}
        constr_coeff['fix_reactions']={}
        Concretemodel_Need_Data['E_total']=e_total
        Concretemodel_Need_Data['B_value']=b_value
        Concretemodel_Need_Data['K_value']=k_value
        constr_coeff['substrate_constrain']=(get_substrate_uptake_id(inputdic),objvalue2)
        constr_coeff['fix_reactions'][inputdic['product']] =  [cond*0.95,np.inf]
        obj_target='maximize'

        EcoECM_FBA_protainmodel=FBA_template2(set_obj_value=True,obj_name=obj_name,obj_target=obj_target,mode=mode,constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
        Model_Solve(EcoECM_FBA_protainmodel)
        print(EcoECM_FBA_protainmodel.obj())
        constr_coeff['fix_reactions'][inputdic['biomass']]=[value(EcoECM_FBA_protainmodel.reaction[inputdic['biomass']]*0.85),np.inf]        
        if mode == 'SET':
            EcoECM_FBA_protainmodel_B2=FBA_template2(set_obj_B_value=True,obj_name=obj_name,obj_target=obj_target,mode=mode,constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
            Model_Solve(EcoECM_FBA_protainmodel_B2)
            B_value2=EcoECM_FBA_protainmodel_B2.obj()
        else:
            B_value2 = Concretemodel_Need_Data.get('B_value', 0)
            print(f"SE mode: skip B optimization in biomass loop at cond={cond}, use B_value={B_value2}")
        print(EcoECM_FBA_protainmodel.obj())
        constr_coeff={}
        constr_coeff['fix_reactions']={}
        Concretemodel_Need_Data['B_value']=B_value2
        # constr_coeff['biomass_constrain']=('EX_lys_L_e',cond*0.1)
        constr_coeff['substrate_constrain']=(get_substrate_uptake_id(inputdic),objvalue2)
        constr_coeff['fix_reactions'][inputdic['product']] =  [cond*0.99,np.inf]
        EcoECM_FBA_protainmodel=FBA_template2(set_obj_value=True,obj_name=obj_name,obj_target=obj_target,mode=mode,constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
        Model_Solve(EcoECM_FBA_protainmodel)
        biomass=EcoECM_FBA_protainmodel.obj()
        print(EcoECM_FBA_protainmodel.obj())
        constr_coeff={}
        constr_coeff['fix_reactions']={}
        constr_coeff['fix_reactions'][inputdic['biomass']]=[biomass*0.9,np.inf]
        constr_coeff['substrate_constrain']=(get_substrate_uptake_id(inputdic),objvalue2)
        constr_coeff['fix_reactions'][inputdic['product']] =  [cond*0.9,np.inf]
        obj_target='minimize'
        obj_name=inputdic['biomass']
        EcoECM_FBA_protainmodel=FBA_template2(set_obj_sum_e=True,obj_name=obj_name,obj_target=obj_target,mode=mode,constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
        Model_Solve(EcoECM_FBA_protainmodel)
        totalE2=EcoECM_FBA_protainmodel.obj()
        print(EcoECM_FBA_protainmodel.obj())

        constr_coeff={}
        constr_coeff['fix_reactions']={}
        # constr_coeff['fix_reactions']['EX_glc__D_e_reverse']=10
        Concretemodel_Need_Data['E_total']=totalE2*1.1
        Concretemodel_Need_Data['B_value']=B_value2*0.9
        Concretemodel_Need_Data['K_value']=k_value
        constr_coeff['substrate_constrain']=(get_substrate_uptake_id(inputdic),objvalue2)
        constr_coeff['fix_reactions'][inputdic['product']] = [cond*0.9,np.inf]
        # constr_coeff['biomass_constrain']=('EX_lys_L_e',cond*0.1)
        obj_target='maximize'
        EcoECM_FBA_protainmodel_max=FBA_template2(set_obj_value=True,obj_name=obj_name,obj_target=obj_target,mode=mode,constr_coeff=constr_coeff,Concretemodel_Need_Data=Concretemodel_Need_Data)
        Model_Solve(EcoECM_FBA_protainmodel_max)
        biomass_max=EcoECM_FBA_protainmodel_max.obj()
        print(EcoECM_FBA_protainmodel_max.obj())   
        e1={x: value(EcoECM_FBA_protainmodel_max.e1[x]) for x in EcoECM_FBA_protainmodel_max.e1}
        model_e1_solution = pd.DataFrame(list(e1.items()), columns=['gene', 'e1'])
        reaction_dict = conbine_flux(EcoECM_FBA_protainmodel_max,model0)
        kcat_dict=Concretemodel_Need_Data['kcat_dict']
        kcat_dict={key:value for key,value in kcat_dict.items()if key in reaction_dict}
        model_reaction_solution = pd.DataFrame(list(reaction_dict.items()), columns=['reaction', 'flux'])
        enzyme_list = list(Concretemodel_Need_Data['mw_dict'].keys()) 
        new_gene_reaction_mapping,gene_reaction_mapping,gene_kcat_mapping=gene_reaction_map2(enzyme_list, reaction_dict,kcat_dict,Concretemodel_Need_Data,model)
        model_e1_solution = pd.DataFrame(list(e1.items()), columns=['gene', 'e1'])
        FSEOFdf['gene'] = model_e1_solution['gene']
        FSEOFdf['gene = '+str(cond)] = FSEOFdf['gene'].map(new_gene_reaction_mapping)
        FSEOFdf['cond = '+str(cond)] = model_e1_solution['e1']
        reactiondf['reaction'] = model_reaction_solution['reaction']
        reactiondf['cond = '+str(cond)] = model_reaction_solution['flux']               
    return FSEOFdf,reactiondf

# def check_monotonicity(row):
#     row_values = row[2:11]  # »ñÈ¡µÚ3µ½µÚ12ÁÐµÄÊý¾Ý
#     increasing = all(row_values[i] < row_values[i + 1] for i in range(len(row_values) - 1))
#     decreasing = all(row_values[i] > row_values[i + 1] for i in range(len(row_values) - 1))

#     if increasing:
#         return 'up'
#     elif decreasing:
#         return 'down'
#     else:
#         return 'unchanged'



def check_monotonicity(row):
    try:
        # Accept either a full row or a subset (e.g., only cond columns)
        numeric_values = []
        for val in row:
            try:
                numeric_values.append(float(val))
            except (ValueError, TypeError):
                continue

        if len(numeric_values) < 2:
            return 'unchanged'

        increasing = all(numeric_values[i] < numeric_values[i + 1] for i in range(len(numeric_values) - 1))
        decreasing = all(numeric_values[i] > numeric_values[i + 1] for i in range(len(numeric_values) - 1))
        all_negative = all(value < 0 for value in numeric_values)

        if all_negative:
            if increasing:
                return 'down'
            elif decreasing:
                return 'up'
            else:
                return 'unchanged'
        else:
            if increasing:
                return 'up'
            elif decreasing:
                return 'down'
            else:
                return 'unchanged'
    except Exception as e:
        print(f'Error in check_monotonicity: {e}')
        return 'unchanged'


def get_gpr_for_reactions(reactiondf, model_input):
    gpr_dict = {}

    for idx, row in reactiondf.iterrows():
        reaction_id = row['reaction']
        gpr = model_input.reactions.get_by_id(reaction_id).gene_reaction_rule
        gpr_dict[reaction_id] = gpr

    return gpr_dict
# 处理 1 到 11 列中存在符号相反的情况
def check_opposite_signs(row):
    row_values = row[1:11]
    if any(row_values[i] * row_values[j] < 0 for i in range(len(row_values)) for j in range(i+1, len(row_values))):
        return True
    return False

def threshold(reactiondf):
    selected_columns = reactiondf.iloc[:, [1, 11]]
    mean_flux = selected_columns.mean(axis=1, numeric_only=True)
    reactiondf['mean_flux'] = mean_flux
    e_2_threshold = abs(1e-1) 
    reactiondf.loc[abs(mean_flux) <= e_2_threshold, 'manipulatios'] = None
    empty_gpr_indices = reactiondf[reactiondf['gpr'] == ''].index
    reactiondf.loc[empty_gpr_indices, 'manipulations'] = None
    up_reactions = reactiondf[reactiondf['manipulations'] == 'up']
    down_reactions = reactiondf[reactiondf['manipulations'] == 'down']

    product_flux = reactiondf.iloc[:, 1:11]
    # 检查并处理数据类型
    product_flux = product_flux.apply(pd.to_numeric, errors='coerce')
    # 计算每行的最大值和最小值的差值
    result_per_row = product_flux.apply(lambda row: round(row.max() - row.min(), 3), axis=1)
    # 将结果添加为新的一列
    reactiondf['result'] = result_per_row
    return up_reactions,down_reactions,result_per_row,reactiondf
def gene_reaction_map3(reaction_list,model_input):
    equation_dict = {}

    for reaction_id in reaction_list:
        equation = model_input.reactions.get_by_id(reaction_id).reaction
        equation_dict[reaction_id] = equation
    return equation_dict

def get_sort_key(x):
    columns = FSEOFdf.columns
    parts = x.split('=')
    if len(parts) == 2 and parts[0].strip() == 'gene':
        return float('inf'), float('inf')
    elif len(parts) == 2 and parts[0].strip() == 'cond':
        return float(parts[1]), float('inf')
    return float('inf'), float('inf')

# def check_monotonicity(row):
#     row_values = row[2:11]  # »ñÈ¡µÚ3µ½µÚ12ÁÐµÄÊý¾Ý
#     increasing = all(row_values[i] < row_values[i + 1] for i in range(len(row_values) - 1))
#     decreasing = all(row_values[i] > row_values[i + 1] for i in range(len(row_values) - 1))

#     if increasing:
#         return 'up'
#     elif decreasing:
#         return 'down'
#     else:
#         return 'unchanged'

def calculate_sum(row, column_name):
    pattern = r'\((.*?)\)'
    matches = re.findall(pattern, row)
    return sum(float(match.split(':')[0]) for match in matches)
def detail(FSEOFdf):
    columns = FSEOFdf.columns
    sorted_columns = sorted(columns, key=get_sort_key)
    # 使用排序后的列名重新构建 DataFrame
    sorted_FSEOFdf = FSEOFdf[sorted_columns]
    FSEOFdf = FSEOFdf.copy()
    cond_columns = [c for c in sorted_columns if c.strip().startswith('cond =')]
    FSEOFdf.loc[:, 'manipulations'] = FSEOFdf[cond_columns].apply(lambda row: check_monotonicity(row), axis=1)

    # FSEOFdf['manipulations'] = FSEOFdf.apply(lambda row: check_monotonicity(row), axis=1)
    # 列名列表
    column_names = FSEOFdf.columns
    # 遍历列名，为每一列创建一个新列来存储括号内数值的总和
    for column_name in column_names:
        # 使用正则表达式提取数值部分
        match = re.search(r'(\d+\.\d+)', column_name)
        
        if match:
            new_column_name = f'sum={match.group(1)}'
            FSEOFdf.loc[:, new_column_name] = FSEOFdf[column_name].astype(str).apply(lambda x: calculate_sum(x, column_name)).copy()
    return FSEOFdf
def result(FSEOFdf):
    cond_columns = [c for c in FSEOFdf.columns if c.strip().startswith('cond =')]
    cond_columns = sorted(cond_columns, key=lambda x: float(x.split('=')[1].strip()))

    FSEOFdf = FSEOFdf.copy()
    if len(cond_columns) < 3:
        FSEOFdf.loc[:, 'mean_fluxs'] = 0.0
        FSEOFdf['enzyme_usage_over'] = FSEOFdf['gene'].map(normalized_E_dict)
        FSEOFdf['enzyme_usage_over'] = pd.to_numeric(FSEOFdf['enzyme_usage_over'], errors='coerce').fillna(0.0)
        FSEOFdf['manipulations'] = None
        return FSEOFdf

    cond_df = FSEOFdf[cond_columns].apply(pd.to_numeric, errors='coerce')
    first_three_mean = cond_df.iloc[:, :3].mean(axis=1)
    last_three_mean = cond_df.iloc[:, -3:].mean(axis=1)

    ratio = last_three_mean / first_three_mean
    ratio = ratio.replace([np.inf, -np.inf], np.nan)
    ratio = ratio.where(ratio > 0)

    fold_change = np.log2(ratio)
    fold_change = fold_change.replace([np.inf, -np.inf], np.nan).fillna(0.0).round(3)

    FSEOFdf.loc[:, 'mean_fluxs'] = fold_change
    FSEOFdf['enzyme_usage_over'] = FSEOFdf['gene'].map(normalized_E_dict)
    FSEOFdf['enzyme_usage_over'] = pd.to_numeric(FSEOFdf['enzyme_usage_over'], errors='coerce').fillna(0.0)

    # Fold_Change threshold rule on log2 scale
    FSEOFdf.loc[(FSEOFdf['manipulations'] == 'up') & (FSEOFdf['mean_fluxs'] <= 1), 'manipulations'] = None
    FSEOFdf.loc[(FSEOFdf['manipulations'] == 'down') & (FSEOFdf['mean_fluxs'] >= -1), 'manipulations'] = None

    # Discard nearly-zero enzyme-cost targets
    low_cost_threshold = 1e-10
    FSEOFdf.loc[FSEOFdf['enzyme_usage_over'].abs() <= low_cost_threshold, 'manipulations'] = None

    FSEOFdf = FSEOFdf.sort_values(by='enzyme_usage_over', ascending=False)
    return FSEOFdf

def output_fseof(path_results):
    filename=os.path.join(path_results, 'results.xlsx') # Ìæ»»ÎªÄãÏë±£´æµÄÎÄ¼þÂ·¾¶
    with pd.ExcelWriter(filename) as writer:
        FSEOFdf.to_excel(writer, sheet_name='test',index=True)
        reactiondf.to_excel(writer, sheet_name='reaction',index=True)

def detail_put(model_input,meged_df):
    reactions = {}
    reaction_ids = [] 
    up_df = meged_df.loc[meged_df['manipulations'] == 'up']
    up_df = up_df.sort_values(by='enzyme_usage_over', ascending=False)
    up_reaction=up_df['gene'].tolist()
    for i in up_reaction:
        try:    
            up_reaction_name = model_input.genes.get_by_id(i).annotation.get('uniprot')
            enzyme_usage_values = up_df.loc[up_df['gene'] == i, 'enzyme_usage_over'].tolist()
            enzyme_usage_values = [float(value) for value in enzyme_usage_values if pd.notna(value)]
            if not enzyme_usage_values:
                enzyme_usage_values = [0.0]

            range_change_value = up_df.loc[up_df['gene'] == i, 'mean_fluxs'].values.tolist()
            range_change_value = round(range_change_value[0], 3) if range_change_value else 0.0
            reaction_ids = []  # 初始化一个空列表
            # reaction_list = model.genes.get_by_id(i).reactions
            # # 将每个反应的 ID 添加到列表中
            # for reaction in reaction_list:
            #     reaction_ids.append(reaction.id)
            reaction_ids = get_reactions_from_expression(i, model0)
        # ½«Ã¿¸ö·´Ó¦ºÍ¶ÔÓ¦µÄ 'range_change' Öµ´æ´¢ÔÚ×ÖµäÖÐ

            reaction_dict = {
                    "Modification": "UP",
                    "Fold_Change":range_change_value ,
                    "Reaction_ID": reaction_ids,
                    "Uniprot_ID": up_reaction_name,
                    "Gene_ID": i,
                    "Enzyme_Cost": enzyme_usage_values,
                    "figpath":os.path.join(path_results3,'target_fig',i+'.png')

                                                    }
            reactions[i] = reaction_dict
        except KeyError:
            pass
    
    down_df = meged_df.loc[meged_df['manipulations'] == 'down']
    down_df = down_df.sort_values(by='enzyme_usage_over', ascending=False)
    down_reaction=down_df['gene'].tolist()
    for i in down_reaction:
            try:
                    down_reaction_name = model_input.genes.get_by_id(i).annotation.get('uniprot')
                    enzyme_usage_values = down_df.loc[down_df['gene'] == i, 'enzyme_usage_over'].tolist()
                    enzyme_usage_values = [float(value) for value in enzyme_usage_values if pd.notna(value)]
                    if not enzyme_usage_values:
                        enzyme_usage_values = [0.0]
                    range_change_value = down_df.loc[down_df['gene'] == i, 'mean_fluxs'].values.tolist()
                    range_change_value = round(range_change_value[0], 3) if range_change_value else 0.0
                    reaction_ids = []  # 初始化一个空列表
                    reaction_ids = get_reactions_from_expression(i, model0)
                    reaction_dict = {
                        "Modification": "Down",
                        "Fold_Change":range_change_value ,
                        "Reaction_ID": reaction_ids,
                        "Uniprot_ID": down_reaction_name,
                        "Gene_ID": i,
                        "Enzyme_Cost": enzyme_usage_values,
                    "figpath":os.path.join(path_results3,'target_fig',i+'.png')
                    }
                    reactions[i] = reaction_dict
            except KeyError:
                    pass
    return reactions,up_reaction,down_reaction

def get_composition_from_formula(met_formula):
    element_re = re.compile("([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")
    composition = {}
    parsed = element_re.findall(met_formula)
    for (element, count) in parsed:
        if count == "":
            count = 1
        else:
            try:
                count = float(count)
                int_count = int(count)
                if count == int_count:
                    count = int_count
                else:
                    warn(
                            "%s is not an integer (in formula %s)"
                            % (count, self.formula)
                        )
            except ValueError:
                warn("failed to parse %s (in formula %s)" % (count, self.formula))
                self.elements = {}
        if element in composition:
            composition[element] += count
        else:
            composition[element] = count
    return (composition)

def calculated_yield(path_task,inputdic,model,prod0):
    solution=model.optimize()
    substrate=solution.fluxes[get_substrate_uptake_id(inputdic)]
    met_C_file = os.path.join(path_task,'met_contain_C_df.tsv')
    met_C = pd.read_csv(met_C_file,sep='\t', index_col=0)
    solution_select=solution.fluxes[abs(solution.fluxes)>0.0000001]
    # ¶¨ÒåÔªËØµÄÏà¶ÔÔ­×ÓÖÊÁ¿×Öµä
    atomic_mass = {'H':1,'He':4,"C":12,"N":14,"O":16,"F":19,'Ne':20,'Na':23,'Mg':24,'Al':27,'Si':28,'P':31,'S':32,'Cl':35.5,'K':39,'Ar':40,'Ca':40,'Mn':55,
                'Fe':56,'Cu':63.5,'Zn':65,'Br':80,'Ag':108,'I':127,'Ba':137,'Pt':195,'Au':197
    }

    # ¼ÆËã²¢½«Ïà¶Ô·Ö×ÓÖÊÁ¿Ìí¼Óµ½ DataFrame ÖÐ
    # ½«ÁÐ×ª»»Îª×Ö·û´®ÀàÐÍ
    met_C['formula'] = met_C['formula'].astype(str)

    # ÖØÐÂ³¢ÊÔ¼ÆËãÏà¶Ô·Ö×ÓÖÊÁ¿
    met_C['Molecular_Weight'] = met_C['formula'].apply(lambda x: sum(get_composition_from_formula(x).get(element, 0) * atomic_mass.get(element, 0) for element in ['H','He',"C","N","O","F",'Ne','Na','Mg','Al','Si','P','S','Cl','K','Ar','Ca','Mn',
                'Fe','Cu','Zn','Br','Ag','I','Ba','Pt','Au']))
    # ¼ÆËã²úÆ·
    product_c_num_r=0
    product_c_num_l=0
    product_molecular_l=0
    product_molecular_r=0
    # get product equation
    product_equ=model.reactions.get_by_id(inputdic['product'])
    # get product metabolites and coefficients
    metabolite_ids = [metabolite.id for metabolite in product_equ.metabolites.keys()]
    coefficients = [product_equ.metabolites[metabolite] for metabolite in product_equ.metabolites.keys()]
    # ½«ÏµÊýºÍ´úÐ»Îï ID ·ÅÔÚÒ»ÆðÐÎ³É×Öµä
    metabolites_coefficients = dict(zip(metabolite_ids, coefficients))
    for eachm in metabolite_ids:
        if metabolites_coefficients[eachm]<0:
            met_coef=abs(metabolites_coefficients[eachm])
            try:
                met_formula=get_composition_from_formula(model.metabolites.get_by_id(eachm).formula)
                met_molecular=met_C.loc[eachm, 'Molecular_Weight']
            except:
                try:
                    met_formula=get_composition_from_formula(met_C.loc[eachm.id,'formula'])
                    met_molecular=met_C.loc[eachm, 'Molecular_Weight']
                except:
                    print('No formula!')
                else:
                    if 'C' in met_formula.keys():
                        #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                        product_c_num_r=product_c_num_r+met_coef*met_formula['C']
                        product_molecular_r=product_molecular_r+met_molecular*met_coef
            else:
                if 'C' in met_formula.keys():
                    #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                    product_c_num_r=product_c_num_r+met_coef*met_formula['C']
                    product_molecular_r=product_molecular_r+met_molecular*met_coef
    for eachm in metabolite_ids:
        if metabolites_coefficients[eachm]>0:
            met_coef=abs(metabolites_coefficients[eachm])
            try:
                met_formula=get_composition_from_formula(model.metabolites.get_by_id(eachm).formula)
                met_molecular=met_C.loc[eachm, 'Molecular_Weight']
            except:
                try:
                    met_formula=get_composition_from_formula(met_C.loc[eachm.id,'formula'])
                    met_molecular=met_C.loc[eachm, 'Molecular_Weight']
                except:
                    print('No formula!')
                else:
                    if 'C' in met_formula.keys():
                        #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                        product_c_num_r=product_c_num_r+met_coef*met_formula['C']
                        product_molecular_r=product_molecular_r+met_molecular*met_coef
            else:
                if 'C' in met_formula.keys():
                    #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                    product_c_num_r=product_c_num_r+met_coef*met_formula['C']
                    product_molecular_r=product_molecular_r+met_molecular*met_coef
   
        product_c_num=product_c_num_r-product_c_num_l
        product_molecular=product_molecular_r-product_molecular_l
        #Ö»ÓÐ²úÆ·Ì¼Ô­×ÓÓÐ¾»Éú³ÉÊ±¼ÆËã 
        if product_c_num>0:
            #µ×ÎïÌ¼Ô­×ÓÊý
            # ¼ÆËãµ×Îï
            #µ×ÎïÌ¼Ô­×ÓÊý
            substrate_c_num_r=0
            substrate_c_num_l=0
            substrate_molecular_l=0
            substrate_molecular_r=0
            substrate_equ=model.reactions.get_by_id(get_substrate_uptake_id(inputdic))
            # get product metabolites and coefficients
            metabolite_ids_s = [metabolite.id for metabolite in substrate_equ.metabolites.keys()]
            coefficients_s = [substrate_equ.metabolites[metabolite] for metabolite in substrate_equ.metabolites.keys()]
            # ½«ÏµÊýºÍ´úÐ»Îï ID ·ÅÔÚÒ»ÆðÐÎ³É×Öµä
            metabolites_coefficients_s = dict(zip(metabolite_ids_s, coefficients_s))
            for eachm in metabolite_ids_s:
                if metabolites_coefficients_s[eachm]<0:
                    met_coef=abs(metabolites_coefficients_s[eachm])
                    try:
                        met_formula=get_composition_from_formula(model.metabolites.get_by_id(eachm).formula)
                        met_molecular=met_C.loc[eachm, 'Molecular_Weight']
                    except:
                        try:
                            met_formula=get_composition_from_formula(met_C.loc[eachm.id,'formula'])
                            met_molecular=met_C.loc[eachm, 'Molecular_Weight']
                        except:
                            print('Substrate right no formula!')
                        else:
                            if 'C' in met_formula.keys():
                                #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                                substrate_c_num_r=substrate_c_num_r+met_coef*met_formula['C']
                                substrate_molecular_r=substrate_molecular_r+met_molecular*met_coef
                    else:
                        if 'C' in met_formula.keys():
                            #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                            substrate_c_num_r=substrate_c_num_r+met_coef*met_formula['C']
                            substrate_molecular_r=substrate_molecular_r+met_molecular*met_coef
                        
            for eachm in metabolite_ids_s:
                if metabolites_coefficients_s[eachm]>0:
                    met_coef=abs(metabolites_coefficients_s[eachm])
                    try:
                        met_formula=get_composition_from_formula(model.metabolites.get_by_id(eachm).formula)
                        met_molecular=met_C.loc[eachm, 'Molecular_Weight']
                    except:
                        try:
                            met_formula=get_composition_from_formula(met_C.loc[eachm.id,'formula'])
                            met_molecular=met_C.loc[eachm, 'Molecular_Weight']
                        except:
                            print('Substrate right no formula!')
                        else:
                            if 'C' in met_formula.keys():
                                #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                                substrate_c_num_r=substrate_c_num_r+met_coef*met_formula['C']
                                substrate_molecular_r=substrate_molecular_r+met_molecular*met_coef
                    else:
                        if 'C' in met_formula.keys():
                            #print(eachm.id,met_coef,met_formula['C'],met_coef*met_formula['C'])
                            substrate_c_num_r=substrate_c_num_r+met_coef*met_formula['C']
                            substrate_molecular_r=substrate_molecular_r+met_molecular*met_coef                                     
            substrate_c_num=substrate_c_num_r-substrate_c_num_l
            substrate_molecular=substrate_molecular_r-substrate_molecular_l
            # if 'EX_glc_e_reverse' in list(flux_solution_df_select.index):
            #     result_static['C Yield']=abs(round((product_c_num*flux_solution_df_select.loc[product_id,'fluxes'])/(substrate_c_num*flux_solution_df_select.loc['EX_glc_e_reverse','fluxes']),7))
            # elif 'EX_glc_e_reverse'+"_reverse" in list(flux_solution_df_select.index):
            #     result_static['C Yield']=abs(round((product_c_num*flux_solution_df_select.loc[product_id,'fluxes'])/(substrate_c_num*flux_solution_df_select.loc['EX_glc_e_reverse'+"_reverse",'fluxes']),7))
            Ycm = substrate_c_num/product_c_num
            if get_substrate_uptake_id(inputdic) in list(list(solution_select.index)):
                Carbon_Yield0=abs(round((product_c_num*prod0)/(substrate_c_num*substrate),3)) 
                Mass_Yield0=abs(round((product_molecular*prod0)/(substrate_molecular*substrate),3))
        # #¼ÆËãµÃÂÊ
        try:
            prod0
        except:
            Yield0=prod0
        else:
            if get_substrate_uptake_id(inputdic) in list(list(solution_select.index)):
                Yield0=abs(round(prod0/solution_select[get_substrate_uptake_id(inputdic)],3))
                Mass_Yield0=abs(round((product_molecular*prod0)/(substrate_molecular*substrate),3)) 
            # elif 'EX_glc_e_reverse'+"_reverse" in list(flux_solution_df_select.index):
            #     result_static['Yield']=abs(round(flux_solution_df_select.loc[product_id,'fluxes']/flux_solution_df_select.loc['EX_glc_e_reverse'+"_reverse",'fluxes'],7))
    return Yield0,Carbon_Yield0,Ycm,Mass_Yield0,product_c_num,substrate_c_num,substrate,substrate_molecular,product_molecular
#  Yield0,Carbon_Yield0,Mass_Yield0,product_c_num,substrate_c_num,substrate = calculated_yield(path_task,inputdic,model,v1_product_max)

def reduced_degree(model,metid):
    product_equ=model.reactions.get_by_id(metid)
# get product metabolites and coefficients
    metabolite_ids = [metabolite.id for metabolite in product_equ.metabolites.keys()]
    met = model.metabolites.get_by_id(metabolite_ids[0])
    degree = 0
    for e,en in met.elements.items():
        if e == 'C':
            degree = degree  + 4*en
        if e == 'H':
            degree = degree  + 1*en
        if e == 'O':
            degree = degree  + (-2)*en
        if e == 'N':
            degree = degree  + (-3)*en
        if e == 'P':
            degree = degree  + (5)*en
        if e == 'S':
            degree = degree  + (6)*en
    # Handle missing charge by defaulting to 0
    charge = met.charge if met.charge is not None else 0
    degree = degree - charge
    degree_C = degree / met.elements['C'] 
    return degree,degree_C

def up_picture(reactiondf):
# 假设first_10_columns是包含反应数据的DataFrame
# 假设model是您的模型对象
    up_df = reactiondf.loc[reactiondf['manipulations'] == 'up']
    # up_df.drop('Unnamed: 0', axis=1, inplace=True)
    # 提取前 10 列
    first_10_columns = up_df.iloc[:, :11]

    # 将 DataFrame 转换为嵌套字典
    result_dict = {}

    for _, row in first_10_columns.iterrows():
        reaction_id = row["gene"].strip()
        reaction_ids = get_reactions_from_expression(reaction_id, model0)
        reaction_title = ";".join(reaction_ids) if reaction_ids else "NA"
        reaction_dict = {"Gene": reaction_id, "Reaction_ID": reaction_title}
        # 提取 'x' 列的值
        exlist = list(np.linspace(0,product,10))
        exlist = [round(element, 3) for element in exlist]
        reaction_dict['x'] = exlist
        # 提取 'y' 列的值
        y_values = [f"{row[col]:.3e}" for col in up_df.columns[1:11]]
        reaction_dict['y'] = y_values
        result_dict[reaction_id] = reaction_dict
    # 将字典保存为 JSON 文件
    up_file = os.path.join(path_results3, 'up.json')
    with open(up_file, 'w') as json_file:
        json.dump(result_dict, json_file, indent=2)
    
    return result_dict

def down_picture(reactiondf):
# 假设first_10_columns是包含反应数据的DataFrame
# 假设model是您的模型对象
    down_df = reactiondf.loc[reactiondf['manipulations'] == 'down']
    # 提取前 10 列
    # down_df.drop('Unnamed: 0', axis=1, inplace=True)
    first_10_columns = down_df.iloc[:, :11]
    # 将 DataFrame 转换为嵌套字典
    result_dict_down = {}

    for _, row in first_10_columns.iterrows():
        reaction_id = row["gene"].strip()
        reaction_ids = get_reactions_from_expression(reaction_id, model0)
        reaction_title = ";".join(reaction_ids) if reaction_ids else "NA"
        reaction_dict = {"Gene": reaction_id, "Reaction_ID": reaction_title}
    #     for col in first_10_columns.columns[1:]:
    #         product_name = col.split("=")[1].strip()
    #         reaction_dict[product_name] = row[col]
    #     result_dict_down[reaction_id] = reaction_dict
        # 提取 'y' 列的值
        # 提取 'x' 列的值
        exlist = list(np.linspace(0,product,10))
        exlist = [round(element, 3) for element in exlist]
        reaction_dict['x'] = exlist
        y_values = [f"{row[col]:.3e}" for col in down_df.columns[1:11]]
        reaction_dict['y'] = y_values
        result_dict_down[reaction_id] = reaction_dict
    # 将字典保存为 JSON 文件
    down_file = os.path.join(path_results3, 'down.json')
    with open(down_file, 'w') as json_file:
        json.dump(result_dict_down, json_file, indent=2)
    
    return result_dict_down
def output_web(reactions):
    output={}
    output['summary'] = {}
    output['yield'] = {}
    output['reaction'] = {}
    up_number = sum(1 for item in reactions.values() if str(item.get("Modification", "")).upper() == "UP")
    down_number = sum(1 for item in reactions.values() if str(item.get("Modification", "")).upper() == "DOWN")
    summary = {
        'UP':up_number,
        'Down':down_number,
        'KO':0,
        'Prod_Rate(mmol/gDW/h)':round(product,3),
        'Growth(1/h)':round(biomass_0,3),
        "Specific_Growth_Rate(1/h)":round(biomass_0,3)
    }
    Yield = {
        'Ycm':Ycm,
        'Ycg':Ycm,
        'Yrm':round(degree_bio/degree_pro,3),
        'Yrg':degree_bio*substrate_molecular/degree_pro*product_molecular,
        'Ysm':Yield0,
        'Ysg':Mass_Yield0,
        'Yscm':Carbon_Yield0,
        'Yscg':Carbon_Yield0
    }
    output['summary'] = summary
    output['yield'] = Yield
    output['genes']=reactions
    output_file = os.path.join(path_results3, 'output.json')
    with open(output_file, 'w') as json_file:
        json.dump(output, json_file) 
    return output


def flux_fseof(model0):
    # Step 1: 提取 up_reaction 中的基因和反应关系
    all_reactions = set()
    reactions_with_ranges = {}

    for gene_expression in up_reaction:  # 直接遍历列表
        related_reactions = get_reactions_from_expression(gene_expression, model0)
        all_reactions.update(related_reactions)
    for gene_expression in down_reaction:  # 直接遍历列表
        related_reactions = get_reactions_from_expression(gene_expression, model0)
        all_reactions.update(related_reactions)
    # Step 2: 将 meged_df 的基因和 Fold_change 对应到反应
    for index, row in FSEOFdf.iterrows():
        gene = row['gene']
        fold_change = row['mean_fluxs']
        if gene in up_reaction:
            related_reactions = get_reactions_from_expression(gene, model0)
            for reaction in related_reactions:
                reactions_with_ranges[reaction] = fold_change
    for index, row in FSEOFdf.iterrows():
        gene = row['gene']
        fold_change = row['mean_fluxs']
        if gene in down_reaction:
            related_reactions = get_reactions_from_expression(gene, model0)
            for reaction in related_reactions:
                reactions_with_ranges[reaction] = fold_change
    # Step 3: 生成反应列表并初始化通量
    reaction_list = [rea.id for rea in model0.reactions]
    fluxnew = {}

    for i in reaction_list:
        fluxnew[i] = 0  # 初始化通量为 0
        for j in reactions_with_ranges.keys():
            if i == j:
                fluxnew[i] += reactions_with_ranges[j]
            if i + '_num' == j[:-1]:
                fluxnew[i] += reactions_with_ranges[j]
            if i + '_reverse' == j or i + '_reverse_num' == j[:-1]:
                fluxnew[i] -= reactions_with_ranges[j]
    flux_file = os.path.join(path_results3, 'flux.json')
    with open(flux_file, 'w') as json_file:
        json.dump(fluxnew, json_file) 
    return fluxnew     

def drawtarget(data,mode,savepath='./'):
    for id in data:
        if mode=='E_FSEOF':
            #  提取数据
            gludy_data = data[id]
            reaction = gludy_data.get('Reaction_ID', gludy_data.get('Gene', 'NA'))
            x_values = gludy_data['x']
            y_values = gludy_data['y']
            # 绘制图形
            plt.plot(x_values, y_values, marker='o', linestyle='-')
            plt.xlabel('Product Flux mmol/gDWh')
            plt.ylabel('Enzyme concentration g/gDW')
            plt.title(id+":"+reaction)
            plt.grid(False)
        # 调整图形边距，防止坐标轴被裁剪
            plt.subplots_adjust(left=0.2)  # 调整左边距，默认是 0.125
            #plt.show()
        elif mode=='E_OptForce':
            gludy_data = data[id]
            reaction = gludy_data['Gene']
            # 定义范围的名称和宽度
            ranges = {'Reference': gludy_data['wild_range'], 'Engineered': gludy_data['eng_range']}
            # 画图
            fig, ax = plt.subplots()
            # 绘制范围
            for i, (name, (start, end)) in enumerate(ranges.items()):
                start, end = float(start), float(end) 
                width = end - start
                ax.plot([i, i], [start, end], color='blue', linewidth=6)
                ax.plot([i-0.05, i+0.05], [start, start], color='blue', linewidth=6)
                ax.plot([i-0.05, i+0.05], [end, end], color='blue', linewidth=6)
                # 标记范围宽度
                #ax.text(i, (start + end) / 2, f'{width}', ha='center', va='center')
            # 设置坐标轴标签
            ax.set_xticks(range(len(ranges)))
            ax.set_xticklabels(list(ranges.keys()))
            #ax.set_yticks(range(1, max([end for _, (_, end) in ranges.items()])+1))
            ax.set_ylabel('Value')
            ax.set_xlabel('Range')
            plt.title(id)
            plt.grid(True)
            #plt.xlabel('Product Flux mmol/gDWh')
            plt.ylabel('Enzyme concentration g/gDW')
            #plt.grid(False)
        # save the fig to the file
        plt.savefig(os.path.join(savepath,id+'.png'))
        # clear the fig
        plt.clf() 

if __name__=="__main__":
    solver_name = os.environ.get("OPTME_PYOMO_SOLVER") or os.environ.get("OPTME_COBRA_SOLVER") or "cplex"
    cobra.Configuration().solver = os.environ.get("OPTME_COBRA_SOLVER", solver_name)
    print(f"Enzyme LP solver: {solver_name}")
    if os.environ.get("OPTME_SOLVER_OPTIONS", "").strip():
        print("Enzyme LP parameters: " + os.environ["OPTME_SOLVER_OPTIONS"])
    path_model=sys.argv[1]
    path_task=sys.argv[2]
    path_map=sys.argv[3]
    path_results=sys.argv[4]
    taskname=sys.argv[5]
    path_results2=os.path.join(path_results,taskname)
    
    path_results3 = os.path.join(path_results2, "E_FSEOF")
    path_results4 = os.path.join(path_results2, "E_OptForce")
    if not os.path.exists(path_results2):
        os.makedirs(path_results2)
    if not os.path.exists(path_results3):
        os.makedirs(path_results3)
    if not os.path.exists(path_results4):
        os.makedirs(path_results4)
    Concretemodel_Need_Data,get_dictionarymodel_data2,Inc,model0,model,dictionary_model,inputdic = prepare_model(path_model,path_task,taskname)
    if inputdic['taskname'] == 'E_OptForce':    
        B_value1,v0_biomass,bio,totalE,objvalue2 = calculate_biomass(Concretemodel_Need_Data,inputdic,model)
        enzyme_list = list(Concretemodel_Need_Data['mw_dict'].keys()) 
        obj_names = enzyme_list
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

        results = pool.starmap(calculate_wildrange, [(Concretemodel_Need_Data,obj_name,inputdic) for obj_name in obj_names])

        pool.close()
        pool.join()


        final_results = {}
        for result in results:
            final_results.update(result)
        
            # 保存结果为 JSON 文件
        wild_file_path = os.path.join(path_results4, 'wild-enzyme.json')
        with open(wild_file_path, 'w') as json_file:
            json.dump(final_results, json_file)
    

        B_value2,v1_product_max,pro,totalE2,EcoECM_FBA_protainmodel_B2 = calculate_product(Concretemodel_Need_Data,inputdic)
        enzyme_list = list(Concretemodel_Need_Data['mw_dict'].keys()) 
        obj_names = enzyme_list
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

        results = pool.starmap(calculate_over, [(Concretemodel_Need_Data,obj_name,inputdic) for obj_name in obj_names])

        pool.close()
        pool.join()


        final_results = {}
        for result in results:
            final_results.update(result)
        over_file_path = os.path.join(path_results4, 'over-enzyme.json')
        with open(over_file_path, 'w') as json_file:
            json.dump(final_results, json_file)
        
        enzyme_results_data,enzyme_overresults2_data = read_file(path_results4)
        ko_data,up_data,down_data,range_change,mean_change = compare_results(enzyme_results_data, enzyme_overresults2_data)
        E_refdict,mw_dict,reaction_dict_bio,kcat_dict = ref_e_con(inputdic,Concretemodel_Need_Data)
        new_gene_reaction_mapping_bio, gene_reaction_mapping, gene_kcat_mapping=gene_reaction_map(enzyme_list, reaction_dict_bio, kcat_dict, Concretemodel_Need_Data, model)
        reaction_dict,kcat_dict,EcoECM_FBA_protainmodel_pro_max = reaction_flux(inputdic,Concretemodel_Need_Data)
        new_gene_reaction_mapping, gene_reaction_mapping,gene_kcat_mapping=gene_reaction_map(enzyme_list,reaction_dict,kcat_dict,Concretemodel_Need_Data,model)
        normalized_E_dict,E_dict,e_dict = enzyme_usage_over(EcoECM_FBA_protainmodel_pro_max,Concretemodel_Need_Data)
        meged_df = must_df(enzyme_results_data,enzyme_overresults2_data,range_change,mean_change,ko_data,up_data,down_data)
        bio_df,pro_df,KeyInfo = basic(bio,pro)
        bottleneck_reactions_list,Df=bottleneck_reactions(EcoECM_FBA_protainmodel_B2,model)
        output(path_results4,bio_df,pro_df,meged_df)
        reactions,up_reaction,down_reaction = detail_put_optforce(model,meged_df)
        Yield0,Carbon_Yield0,Ycm,Mass_Yield0,product_c_num,substrate_c_num,substrate,substrate_molecular,product_molecular = calculated_yield(path_task,inputdic,model,v1_product_max)
        degree_pro,degree_C_pro =reduced_degree(model,inputdic['product'])
        degree_bio,degree_C_bio =reduced_degree(model,get_substrate_uptake_id(inputdic))
        output = output_web_optforce(reactions)
        result_json = up_picture_optforce(meged_df)
        result_json = down_picture_optforce(meged_df)
        os.makedirs(os.path.join(path_results4,'target_fig'),exist_ok=True)
        with open(os.path.join(path_results4, 'up.json'),encoding='utf-8') as fp:
            result_json2=json.load(fp)
        if isinstance(result_json2, dict) and len(result_json2) > 0:
            drawtarget(result_json2, 'E_OptForce', os.path.join(path_results4, 'target_fig'))
        else:
            print("Warning: 'up' is empty or not a dictionary. Skipping drawtarget.")
        with open(os.path.join(path_results4, 'down.json'),encoding='utf-8') as fp:
            result_json3=json.load(fp)        
        if isinstance(result_json3, dict) and len(result_json3) > 0:
            drawtarget(result_json3, 'E_OptForce', os.path.join(path_results4, 'target_fig'))
        else:
            print("Warning: 'down' is empty or not a dictionary. Skipping drawtarget.")
        #fluxnew = flux(model0)
    if inputdic['taskname'] == 'E_FSEOF':
        product,objvalue2,EcoECM_FBA_protainmodel_pro = calculate_product_fseof(Concretemodel_Need_Data,inputdic)
        biomass_0,objvalue2,EcoECM_FBA_protainmodel_bio = calculate_product_fseof_bio(Concretemodel_Need_Data,inputdic)
        FSEOFdf,reactiondf = biomass(product,inputdic)
        columns = FSEOFdf.columns
        sorted_columns = sorted(columns, key=get_sort_key)
        sorted_FSEOFdf = FSEOFdf[sorted_columns]
        FSEOFdf = sorted_FSEOFdf[['gene'] + [col for col in sorted_columns if col != 'gene']]
        # FSEOFdf.loc[:, 'manipulations'] = FSEOFdf.apply(lambda row: check_monotonicity(row), axis=1).copy()
        # # 列名列表
        # column_names = FSEOFdf.columns

        # # 遍历列名，为每一列创建一个新列来存储括号内数值的总和
        # for column_name in column_names:
        #     # 使用正则表达式提取数值部分
        #     match = re.search(r'(\d+\.\d+)', column_name)
            
        #     if match:
        #         new_column_name = f'sum={match.group(1)}'
        #         FSEOFdf.loc[:, new_column_name] = FSEOFdf[column_name].astype(str).apply(lambda x: calculate_sum(x, column_name)).copy()
        # last_10_columns = FSEOFdf.iloc[:, -10:]
        # mean_flux = last_10_columns.mean(axis=1)
        # FSEOFdf.loc[:, 'mean_fluxs'] = mean_flux.copy()
        # e_2_threshold = 1e-1
        # FSEOFdf.loc[mean_flux <= e_2_threshold, 'manipulations'] = None
        normalized_E_dict,E_dict,e_dict = enzyme_usage_over(EcoECM_FBA_protainmodel_pro,Concretemodel_Need_Data)
        FSEOFdf =detail(FSEOFdf)
        FSEOFdf =result(FSEOFdf)
        output_fseof(path_results3)
        reactions,up_reaction,down_reaction = detail_put(model,FSEOFdf)
        Yield0,Carbon_Yield0,Ycm,Mass_Yield0,product_c_num,substrate_c_num,substrate,substrate_molecular,product_molecular = calculated_yield(path_task,inputdic,model,product)
        degree_pro,degree_C_pro =reduced_degree(model,inputdic['product'])
        degree_bio,degree_C_bio =reduced_degree(model,get_substrate_uptake_id(inputdic))
        result_dict = up_picture(FSEOFdf)
        result_dict_down = down_picture(FSEOFdf)
        os.makedirs(os.path.join(path_results3,'target_fig'),exist_ok=True)
        if len(result_dict)>0:
            drawtarget(result_dict,'E_FSEOF', os.path.join(path_results3,'target_fig'))
        if len(result_dict_down)>0:
            drawtarget(result_dict_down,'E_FSEOF', os.path.join(path_results3,'target_fig'))
        output = output_web(reactions)
        #fluxnew = flux_fseof(model0)


# python EToptme.py '/hpcfs/fhome/xuwenqi/project/OptMetarget/INPUT/model' '/hpcfs/fhome/xuwenqi/project/OptMetarget/INPUT/task' '/hpcfs/fhome/xuwenqi/project/OptMetarget/INPUT/map' '/hpcfs/fhome/xuwenqi/project/OptMetarget/OUTPUT' 'wangry@tib.cas.cn_20241128-131412_iCW773R'
# python EToptme.py '/hpcfs/fhome/xuwenqi/project/optme/ET-OptME/input/model/iCW' '/hpcfs/fhome/xuwenqi/project/optme/ET-OptME/input/task' '/hpcfs/fhome/xuwenqi/project/OptMetarget/input1/target/map' '/hpcfs/fhome/xuwenqi/project/optme/ET-OptME/output' 'ala1_iCW_F'
