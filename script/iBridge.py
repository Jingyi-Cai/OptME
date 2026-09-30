
# Adapted from iBridge （Kim 2023, Cell Systems,https://doi.org/10.1016/j.cels.2023.10.005)
import os
import glob
import time
from copy import deepcopy
import sys
import os
import re 
import warnings
import json
import pandas as pd
import numpy as np
import shutil
import cobra
from cobra.io import read_sbml_model, write_sbml_model, load_json_model
from cobra import Model, Reaction, Metabolite
# Flux sampling uses the same COBRA solver as target.py (OPTME_COBRA_SOLVER, default cplex).
from visualization2 import *
from solver_budget import install_cobra_threads, use_sequential_threads

def _load_cobra_model(path_model, model_name):
    model_file = os.path.join(path_model, model_name)
    if ".mat" in model_name:
        return cobra.io.load_matlab_model(model_file)
    if ".xml" in model_name or ".sbml" in model_name:
        return cobra.io.read_sbml_model(model_file)
    if ".json" in model_name:
        return cobra.io.load_json_model(model_file)
    return cobra.io.load_model(model_name)


def _apply_task_bounds(model, inputdic):
    """Set oxygen and CO2 bounds from the task JSON, as target.py does."""
    o2_reaction_id = inputdic.get("O2", "EX_o2_e")
    try:
        o2_reaction = model.reactions.get_by_id(o2_reaction_id)
    except KeyError:
        print("Warning: oxygen reaction %s not found, skip oxygen bound adjustment." % o2_reaction_id)
    else:
        if inputdic.get("oxygenstate") == "aerobic":
            o2_reaction.lower_bound = -1000
        elif inputdic.get("oxygenstate") == "micro_aerobic":
            o2_reaction.lower_bound = -2
        elif inputdic.get("oxygenstate") == "anaerobic":
            o2_reaction.lower_bound = 0
    try:
        co2 = model.reactions.get_by_id("EX_co2_e")
    except KeyError:
        co2 = None
    if co2 is not None:
        if str(inputdic.get("CO2_flag")) == "True":
            co2.lower_bound = -1000
        else:
            co2.lower_bound = 0


def prepare_model(path_model, path_task, taskname):
    with open(os.path.join(path_task, taskname) + ".json", encoding="utf-8") as fp:
        inputdic = json.load(fp)
    model = _load_cobra_model(path_model, inputdic["model"])
    _apply_task_bounds(model, inputdic)
    universal_model_file = os.path.join(path_task, "BiGG.xml")
    if not os.path.isfile(universal_model_file):
        raise FileNotFoundError("BiGG.xml was not found in the task directory: " + universal_model_file)
    model0 = read_sbml_model(universal_model_file)
    return model0, model, inputdic


def _flux_dict(solution):
    fluxes = {}
    for reaction_id, value in solution.fluxes.items():
        value = float(value)
        if abs(value) <= 1e-6:
            value = 0.0
        fluxes[reaction_id] = value
    return fluxes


def _optimize_cobra(model, objective, constraints, mode="max", parsimonious=False):
    """FBA, optionally followed by parsimonious FBA. Status 2 means optimal."""
    with model:
        for reaction_id, bounds in constraints.items():
            reaction = model.reactions.get_by_id(reaction_id)
            reaction.lower_bound = bounds[0]
            reaction.upper_bound = bounds[1]
        model.objective = objective
        model.objective_direction = "max" if mode == "max" else "min"
        if parsimonious:
            solution = cobra.flux_analysis.pfba(model)
        else:
            solution = model.optimize()
    if solution.status != "optimal":
        return 0, None, None
    return 2, solution, _flux_dict(solution)


def run_MetScore_simulation(model, biomass_rxn, target_rxn, constdic=None):
    """Covariance scan used by iBridge, solved with the COBRA solver."""
    if constdic is None:
        constdic = {}
    model.solver = os.environ.get("OPTME_COBRA_SOLVER", "cplex")
    stat, wild_solution, wild_opt_flux = _optimize_cobra(
        model, biomass_rxn, constdic, mode="max", parsimonious=True
    )
    if stat != 2:
        raise RuntimeError("Wild-type biomass FBA failed for iBridge.")
    print("#constraints", constdic)

    _, _, flux_dic = _optimize_cobra(
        model, target_rxn, constdic, mode="max", parsimonious=True
    )
    if flux_dic is None:
        raise RuntimeError("Product maximization failed for %s." % target_rxn)
    max_const = flux_dic[target_rxn]

    _, _, flux_dic = _optimize_cobra(
        model, target_rxn, constdic, mode="min", parsimonious=True
    )
    if flux_dic is None:
        raise RuntimeError("Product minimization failed for %s." % target_rxn)
    min_const = flux_dic[target_rxn]

    flux_dist_dic_set = {}
    count = 0
    for each_flux_const in np.linspace(min_const, max_const, 10):
        tmp_const = deepcopy(constdic)
        count += 1
        tmp_const[target_rxn] = [each_flux_const * 0.95, each_flux_const * 1.05]
        stat, _, opt_flux_dic = _optimize_cobra(
            model, biomass_rxn, tmp_const, mode="max", parsimonious=False
        )
        if stat != 2:
            print("Biomass optimization failed at count %d" % count)
            continue

        tmp_biomass_flux = opt_flux_dic[biomass_rxn]
        tmp_const[biomass_rxn] = [tmp_biomass_flux * 0.95, tmp_biomass_flux * 1.05]
        print("Count: %d\tBiomass flux: %0.6f" % (count, tmp_biomass_flux))

        with model:
            for reaction_id, bounds in tmp_const.items():
                reaction = model.reactions.get_by_id(reaction_id)
                reaction.lower_bound = bounds[0]
                reaction.upper_bound = bounds[1]
            moma_solution = cobra.flux_analysis.moma(model, solution=wild_solution, linear=True)
        if moma_solution.status == "optimal":
            flux_dic = _flux_dict(moma_solution)
            flux_dist_dic_set[each_flux_const] = flux_dic
            if abs(flux_dic.get(biomass_rxn, 0.0)) <= 1e-6:
                break

    if not flux_dist_dic_set:
        raise RuntimeError("iBridge produced no feasible flux samples.")
    flux_df = pd.DataFrame.from_dict(flux_dist_dic_set)
    flux_corr_df = flux_df.abs().T.corr()
    flux_cov_df = flux_df.abs().T.cov()
    return flux_df, flux_corr_df, flux_cov_df


def calculate_MetScore_sum(cobra_model, covariance_data):
    branch_metabolite_data = {}
    neighboring_rxns = {}

    for each_metabolite in cobra_model.metabolites:
        target_met_id = each_metabolite.id
        branch_metabolite_data[target_met_id] = 0.0
        count = 0
        for each_reaction in each_metabolite.reactions:
            reactants = [met.id for met in each_reaction.reactants]
            if target_met_id in reactants:
                if each_reaction.id in covariance_data:
                    # branch_metabolite_data[target_met_id] += abs(covariance_data[each_reaction.id])
                    try:
                        value = float(covariance_data[each_reaction.id])
                        branch_metabolite_data[target_met_id] += abs(value)
                    except ValueError:
                        # 如果转换失败，打印错误信息或选择跳过该值
                        print(f"Non-numeric data for reaction {each_reaction.id}: {covariance_data[each_reaction.id]}")
                         
    
    return branch_metabolite_data


def select_candidates(cobra_model, target_reaction, output_file, 
                      metscore_df, corr_df, cov_df, 
                      corr_threshold=0, cov_threshold=0.1):
    # Ensure that corr_df and cov_df are numeric
    corr_df = corr_df.apply(pd.to_numeric, errors='coerce')
    cov_df = cov_df.apply(pd.to_numeric, errors='coerce')

    # Filter the data based on thresholds
    pcorr_df = corr_df[corr_df > corr_threshold]
    pcov_df = cov_df[cov_df > cov_threshold]
    pcorr_df = pcorr_df.dropna()
    pcov_df = pcov_df.dropna()
    
    positive_candidate_reactions = list(set(pcorr_df.index) & set(pcov_df.index))
    
    ncorr_df = corr_df[corr_df < -corr_threshold]
    ncov_df = cov_df[cov_df < -cov_threshold]
    ncorr_df = ncorr_df.dropna()
    ncov_df = ncov_df.dropna()
    
    negative_candidate_reactions = list(set(ncorr_df.index) & set(ncov_df.index))
    
    print('Number of positive candidate reactions: %d' % (len(positive_candidate_reactions)))
    print('Number of negative candidate reactions: %d' % (len(negative_candidate_reactions)))
    
    with open(output_file, 'w') as fp:
        header = ['Metabolite', 'Score', 'Normalized score', 
                  'No. of reactions', 'No. of positive reactions', 
                  'No. of negative reactions', 'candidate reactions', 
                  'positive reactions', 'negative reactions', 
                  'Positive score', 'Negative score']
        fp.write('%s\n' % ('\t'.join(header)))
        
        for each_row, each_df in metscore_df.iterrows():
            cobra_metabolite = cobra_model.metabolites.get_by_id(each_row)
            candidate_reactions = []

            # Collect candidate reactions for the current metabolite
            for each_reaction in cobra_metabolite.reactions:
                for each_reactant in each_reaction.reactants:
                    if each_reactant.id == each_row:
                        candidate_reactions.append(each_reaction.id)
            
            candidate_reactions = list(set(candidate_reactions))
            pos_candidate_reactions = list(set(positive_candidate_reactions) & set(candidate_reactions))
            neg_candidate_reactions = list(set(negative_candidate_reactions) & set(candidate_reactions))
            
            positive_score = 0.0
            negative_score = 0.0
            for rxn in pos_candidate_reactions:
                positive_score += float(pcov_df.loc[rxn])
            for rxn in neg_candidate_reactions:
                negative_score += float(ncov_df.loc[rxn])
            
            # Compute normalized score
            if each_df[target_reaction] != 0.0 and np.sqrt(float(len(candidate_reactions))) != 0.0:
                normalized_score = each_df[target_reaction] / np.sqrt(float(len(candidate_reactions)))
            else:
                normalized_score = 0.0
            
            # Prepare the row to write to file
            contents = [each_row, each_df[target_reaction], normalized_score,
                        len(candidate_reactions), len(pos_candidate_reactions),
                        len(neg_candidate_reactions), ';'.join(candidate_reactions),
                        ';'.join(pos_candidate_reactions), ';'.join(neg_candidate_reactions),
                        positive_score, negative_score]
            
            for i, item in enumerate(contents):
                delim = '\t' if i < len(contents) - 1 else '\n'
                fp.write(str(item) + delim)
    
    return



def read_target_rxn(filename):
    fp = open(filename, 'r')
    lines = fp.read()
    lines = lines.replace('[', '')
    lines = lines.replace(']', '')
    lines = lines.replace("'", "")
    fp.close()
    target_rxn = lines.strip().split(',')
    print(len(target_rxn))
    return target_rxn
    

def write_header(filename, taget_id):
    fp = open(filename, 'r')
    lines = fp.read()
    fp.close()
    fp = open(filename, 'w')
    fp.write('%s,%s\n'%('reaction', taget_id))
    fp.write(lines.strip())
    fp.close()
    return

# Module 2

def make_candidate_reaction_sets(df, basename, input_dir, path_results3):
    ex_metabolites = []
    with open('%s/excluded_metabolites.txt'%input_dir, 'r') as fp:
        for line in fp:
            ex_metabolites.append(line.strip())
            
    score_info = {}
    for each_met, each_df in df.groupby('Metabolite'):
            
        if each_met[-1] == 'c':
            pos_reaction_num = each_df['No. of positive reactions'].values[0]
            neg_reaction_num = each_df['No. of negative reactions'].values[0]

            pos_score = each_df['Positive score'].values[0]
            neg_scroe = each_df['Negative score'].values[0]
            
            final_pos_scroe = pos_score/np.sqrt(1+pos_reaction_num)
            final_neg_scroe = neg_scroe/np.sqrt(1+neg_reaction_num)
            score_info[each_met] = [final_pos_scroe, final_neg_scroe]

    negative_score_info = {}
    positive_score_info = {}
    
    for met in score_info:
        if abs(score_info[met][0]) > abs(score_info[met][1]):
            positive_score_info[met] = score_info[met][0]
        else:
            negative_score_info[met] = score_info[met][1]
      
    if not os.path.exists('%s/application_results'%path_results3):
        os.makedirs('%s/application_results'%path_results3)
    fp = open('%s/application_results/Candidate_%s'%(path_results3, basename), 'w')
    fp.write('%s\t%s\t%s\t%s\n'\
             %('Negative metabolite', 'Positive metabolite',
               'Negative score', 'Positive score'))
    for negative_met in negative_score_info:
        negative_score = negative_score_info[negative_met]
        
        for positive_met in positive_score_info:
            
            if negative_met in ex_metabolites or positive_met in ex_metabolites:
                continue
            
            positive_score = positive_score_info[positive_met]
            fp.write('%s\t%s\t%s\t%s\n'\
                     %(negative_met, positive_met, negative_score, positive_score))
            
    fp.close()
    return
    
    
def check_reaction(filename, original_model_met_info, universal_model_met_info, flux_corr_df, flux_cov_df, path_results3):
    basename = os.path.basename(filename)
    ex_reaction_id = basename.split('Final_MetScore_')[1].strip()
    ex_reaction_id = ex_reaction_id.replace('.txt', '')
    
    df = pd.read_table(filename)
    if not os.path.exists('%s/application_result2'%path_results3):
        os.makedirs('%s/application_result2'%path_results3)
    fp = open('%s/application_result2/New_reaction_candidate_up_%s'%(path_results3, basename), 'w')
    fp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'\
             %('Target', 'Gene source', 'Negative metabolite', 
              'Positive metabolite', 'Negative score', 'Positive score', 
              'reaction', 'Equation', 'Corr', 'Cov'))
    
    for each_row, each_df in df.iterrows():
        negative_met = each_df['Negative metabolite']
        positive_met = each_df['Positive metabolite']

        negative_score = each_df['Negative score']
        positive_score = each_df['Positive score']
        
        for each_met_set in original_model_met_info:
            if negative_met in each_met_set[0] and positive_met in each_met_set[1]:
                target_reaction = each_met_set[2]
                if target_reaction not in flux_corr_df.index:
                    fp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'\
                             %(ex_reaction_id, 'Native', negative_met, positive_met, 
                              negative_score, positive_score, each_met_set[2], 
                              each_met_set[3], 'NA', 'NA'))
                else:
                    corr_val = float(flux_corr_df.loc[target_reaction])
                    cov_val = float(flux_cov_df.loc[target_reaction])
                    fp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'\
                             %(ex_reaction_id, 'Native', negative_met, positive_met, 
                              negative_score, positive_score, each_met_set[2], 
                              each_met_set[3], corr_val, cov_val))
                    
        
        for each_met_set in universal_model_met_info:
            if negative_met+'_' in each_met_set[0] and positive_met+'_' in each_met_set[1]:
                target_reaction = each_met_set[2]
                if target_reaction not in flux_corr_df.index:
                    fp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'\
                             %(ex_reaction_id, 'Universal model', negative_met, positive_met, 
                              negative_score, positive_score, each_met_set[2], 
                              each_met_set[3], 'NA', 'NA'))
    fp.close()
    return


def check_down_reaction(filename, original_model_met_info, universal_model_met_info, flux_corr_df, flux_cov_df, path_results3):
    basename = os.path.basename(filename)
    ex_reaction_id = basename.split('Final_MetScore_')[1].strip()
    ex_reaction_id = ex_reaction_id.replace('.txt', '')
    
    df = pd.read_table(filename)
    if not os.path.exists('%s/application_result2'%path_results3):
        os.makedirs('%s/application_result2'%path_results3)
    fp = open('%s/application_result2/New_reaction_candidate_down_%s'%(path_results3, basename), 'w')
    fp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'\
             %('Target', 'Gene source', 'Negative metabolite', 
              'Positive metabolite', 'Negative score', 'Positive score', 
              'reaction', 'Equation', 'Corr', 'Cov'))
    
    for each_row, each_df in df.iterrows():
        negative_met = each_df['Negative metabolite']
        positive_met = each_df['Positive metabolite']

        negative_score = each_df['Negative score']
        positive_score = each_df['Positive score']
        
        for each_met_set in original_model_met_info:
            if negative_met in each_met_set[1] and positive_met in each_met_set[0]:
                target_reaction = each_met_set[2]
                if target_reaction not in flux_corr_df.index:
                    fp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'\
                             %(ex_reaction_id, 'Native_down', negative_met, positive_met, 
                              negative_score, positive_score, each_met_set[2], 
                              each_met_set[3], 'NA', 'NA'))
                else:
                    corr_val = float(flux_corr_df.loc[target_reaction])
                    cov_val = float(flux_cov_df.loc[target_reaction])
                    fp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'\
                             %(ex_reaction_id, 'Native_down', negative_met, positive_met, 
                              negative_score, positive_score, each_met_set[2], 
                              each_met_set[3], corr_val, cov_val))
    fp.close()
    return





def metabolite_set(cobra_model):
    met_set_info = []
    for each_reaction in cobra_model.reactions:
        reactants = [met.id for met in each_reaction.reactants]
        products = [met.id for met in each_reaction.products]
        met_set_info.append([reactants, products, each_reaction.id, each_reaction.reaction])
    return met_set_info


def metabolite_set_cytosol(cobra_model):
    met_set_info = []
    for each_reaction in cobra_model.reactions:
        reactants = [met.id for met in each_reaction.reactants]
        products = [met.id for met in each_reaction.products]

        compartments = []
        for each_met in reactants+products:
            each_cmp = each_met[-2:]
            compartments.append(each_cmp)
            
        compartments = list(set(compartments))
        if len(compartments) == 1:
            if compartments[0] == 'c_':
                met_set_info.append([reactants, products, each_reaction.id, each_reaction.reaction])
        
    return met_set_info




# Module 3
def parse_reactions(reaction, ex_metabolites):
    reaction = reaction.replace('>', '')
    reaction = reaction.replace('<', '')
    
    sptlist = reaction.split('--')
    reactants = sptlist[0].strip().split(' + ')
    products = sptlist[0].strip().split(' + ')
    new_reactants = []
    new_products = []
    
    for met in reactants:
        if met not in ex_metabolites:
            new_reactants.append(met)
            
    for met in products:
        if met not in ex_metabolites:
            new_products.append(met)
            
    return new_reactants, new_products

def check_duplicate(reaction, known_reaction_info, ex_metabolites):
    flag = False
    reactants, products = parse_reactions(reaction, ex_metabolites)
    if len(reactants) > 0 and len(products) > 0: 
        for rxn in known_reaction_info:
            reactants2 = known_reaction_info[rxn][0]
            products2 = known_reaction_info[rxn][1]

            if  set(reactants).issubset(set(reactants2)) and set(products).issubset(set(products2)):
                flag = True
                return flag

            if  set(products).issubset(set(reactants2)) and set(reactants).issubset(set(products2)):
                flag = True
                return flag
        
    return flag

def detail_put(model_input):
    reactions_detail = {}
    df_up = pd.read_csv(os.path.join(path_results3,'Application_final_up_summary.txt'),delimiter='\t')
    df_down = pd.read_csv(os.path.join(path_results3,'Application_final_down_summary.txt'),delimiter='\t')
    df_up['fc'] = abs(df_up['Positive score'] - df_up['Negative score']).round(3)
    df_down['fc'] = abs(df_down['Positive score'] - df_down['Negative score']).round(3)
        # 移除 fc 小于 1 的行
    df_up = df_up[df_up['fc'] >= 1]
    df_down = df_down[df_down['fc'] >= 1]
    # df_up = meged_df.loc[meged_df['manipulation'] == 'up']
    # df_up = df_up.sort_values(by='result', ascending=False)
    df_up = df_up[['reaction','Equation','fc','Positive score']]
    up_reaction=df_up['reaction'].tolist()
    for i in up_reaction:
        try:    
            reaction = model_input.reactions.get_by_id(i)
            up_reaction_name = reaction.name
            up_reaction_equation = reaction.reaction
            gpr=model_input.reactions.get_by_id(i).gene_reaction_rule
            range_change_value = df_up.loc[df_up['reaction'] == i, 'fc'].values.tolist()
            range_change_value = range_change_value[0] if range_change_value else None
            reaction_dict = {
                    "Modification": "UP",
                    "Covariance_Difference":range_change_value ,
                    "Reaction_ID": i,
                    "Reaction_Name": up_reaction_name,
                    "Gene_ID": gpr,
                    "Equation": up_reaction_equation
                }
            reactions_detail[i] = reaction_dict
        except KeyError:
            pass
    down_df = df_down[['reaction','Equation','fc','Positive score']]
    down_reaction=down_df['reaction'].tolist()
    for i in down_reaction:
        try:
            reaction = model_input.reactions.get_by_id(i)
            down_reaction_name = reaction.name
            down_reaction_equation = reaction.reaction
            gpr=model_input.reactions.get_by_id(i).gene_reaction_rule
            range_change_value = down_df.loc[down_df['reaction'] == i, 'fc'].values.tolist()
            range_change_value = range_change_value[0] if range_change_value else None
            reaction_dict = {
                "Modification": "Down",
                "Covariance_Difference":range_change_value ,
                "Reaction_ID": i,
                "Reaction_Name": down_reaction_name,
                "Gene_ID": gpr,
                "Equation": down_reaction_equation
            }
            reactions_detail[i] = reaction_dict
        except KeyError:
                pass
    range_change = []
    # 遍历 df_up 表格中的每一行
    for index, row in df_up.iterrows():
        # 获取 reaction 和 Positive score 的值
        reaction = row['reaction']
        positive_score = row['Positive score']
        
        # 将 (reaction, Positive score) 元组添加到 range_change 列表
        range_change.append((reaction, positive_score))
    return reactions_detail,range_change,up_reaction,down_reaction,df_up,df_down
    
def calculate_prod(model, inputdic):  
    model.objective=inputdic['biomass']
    v0_biomass=model.optimize().objective_value
    model.objective=inputdic['product']
    prod = model.optimize().objective_value
    return prod,v0_biomass

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
    Yield0 = 0
    Carbon_Yield0 = 0
    Ycm = 0
    Mass_Yield0 = 0
    solution=model.optimize()
    substrate=solution.fluxes[inputdic['substrate']]
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
            substrate_equ=model.reactions.get_by_id(inputdic['substrate'])
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
            if inputdic['substrate'] in list(list(solution_select.index)):
                Carbon_Yield0=abs(round((product_c_num*prod0)/(substrate_c_num*substrate),3)) 
                Mass_Yield0=abs(round((product_molecular*prod0)/(substrate_molecular*substrate),3))
        # #¼ÆËãµÃÂÊ
        try:
            prod0
        except:
            Yield0=prod0
        else:
            if inputdic['substrate'] in list(list(solution_select.index)):
                Yield0=abs(round(prod0/solution_select[inputdic['substrate']],3))
                Mass_Yield0=abs(round((product_molecular*prod0)/(substrate_molecular*substrate),3)) 
            # elif 'EX_glc_e_reverse'+"_reverse" in list(flux_solution_df_select.index):
            #     result_static['Yield']=abs(round(flux_solution_df_select.loc[product_id,'fluxes']/flux_solution_df_select.loc['EX_glc_e_reverse'+"_reverse",'fluxes'],7))
    return Yield0,Carbon_Yield0,Ycm,Mass_Yield0,product_c_num,substrate_c_num,substrate,substrate_molecular,product_molecular
        # Yield0,Carbon_Yield0,Ycm,Mass_Yield0,product_c_num,substrate_c_num,substrate,substrate_molecular,product_molecular = calculated_yield(path_task,inputdic,model,prod0)

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
    # 处理 met.charge 为 None 的情况
    if met.charge is None:
        met.charge = 0
    degree -= met.charge

    # 避免除零错误
    if 'C' not in met.elements or met.elements['C'] == 0:
        raise ValueError("The metabolite does not contain any carbon atoms.")

    degree_C = degree / met.elements['C']
    return degree,degree_C

def output_web_ibridge(reactions_detail,prod0):
    output={}
    output['summary'] = {}
    output['reaction'] = {}
    up_number=len(up_reaction)
    down_number=len(down_reaction)
    summary = {
        'UP':up_number,
        'Down':down_number,
        'KO':0
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
    output['reaction']=reactions_detail
    output_file = os.path.join(path_results3, 'output.json')
    with open(output_file, 'w') as json_file:
        json.dump(output, json_file) 
    
    # filename=os.path.join(path_results2, 'results.xlsx') # Ìæ»»ÎªÄãÏë±£´æµÄÎÄ¼þÂ·¾¶
    # with pd.ExcelWriter(filename) as writer:
    #     meged_df.to_excel(writer, sheet_name='test',index=True)
    return output

def conbine_flux(modelid):   
    #get flux of stoichiometry from E_protain_FBA
    # # Merge up_reaction and down_reaction
    all_reactions = up_reaction + down_reaction
    result_reaction_columns = list(zip(df_up['reaction'], df_up['fc']))
    result_reaction_columns_down = list(zip(df_down['reaction'], df_down['fc']))
    reactions_with_ranges = {}
    for reaction, range_value in result_reaction_columns:
        if reaction in up_reaction:
            reactions_with_ranges[reaction] = range_value
    for reaction, range_value in result_reaction_columns_down:
        if reaction in down_reaction:
            reactions_with_ranges[reaction] = -range_value            
    flux_file = os.path.join(path_results3, 'fcflux_map.json')
    values = []
    for key, value in reactions_with_ranges.items():
        values.append(value)
    if len(values) == 0:
        mean = 0.0
        std_dev = 0.0
    else:
        mean = float(np.mean(values))
        std_dev = float(np.std(values))
    reactions_with_ranges['range_lb'] = round(mean - 2 * std_dev, 3)
    reactions_with_ranges['range_ub'] = round(mean + 2 * std_dev, 3)
    with open(flux_file, 'w') as json_file:
        json.dump(reactions_with_ranges, json_file) 
    return reactions_with_ranges          



def calculate_product(inputdic):
    model.objective=inputdic['product']
    v1_product_max=model.optimize().objective_value
    product = model.optimize()
    model_pfba_solution = cobra.flux_analysis.pfba(model)
    model_solution_frame = model_pfba_solution.to_frame()
    over_flux = model_solution_frame.fluxes.to_dict()  
    values = []
    for key, value in over_flux.items():
        values.append(value)
    # Calculate mean and standard deviation
    mean = np.mean(values)
    std_dev = np.std(values)
    over_flux['range_lb']=(mean-2*std_dev).round(3)
    over_flux['range_ub']=(mean + 2 * std_dev).round(3) 
    over_flux_path = os.path.join(path_results3, 'overflux_map.json')
    with open(over_flux_path, 'w') as file:
        json.dump(over_flux, file)
    return v1_product_max,over_flux

def calculate_biomass(inputdic):
    model.objective=inputdic['biomass']
    biomass=model.optimize().objective_value
    product = model.optimize()
    model_pfba_solution = cobra.flux_analysis.pfba(model)
    model_solution_frame = model_pfba_solution.to_frame()
    wild_flux = model_solution_frame.fluxes.to_dict()  
    values = []
    for key, value in wild_flux.items():
        values.append(value)
    # Calculate mean and standard deviation
    mean = np.mean(values)
    std_dev = np.std(values)
    wild_flux['range_lb']=(mean-2*std_dev).round(3)
    wild_flux['range_ub']=(mean + 2 * std_dev).round(3)     
    wild_flux_path = os.path.join(path_results3, 'wildflux_map.json')
    with open(wild_flux_path, 'w') as file:
        json.dump(wild_flux, file)    
    return biomass,wild_flux


def reaction_list_d3(model,inputdic):
    # model.reactions.get_by_id(inputdic['substrate']).bounds=(-inputdic['concentration'],0)
    model.reactions.get_by_id('ATPM').bounds =(0,0)
    model.objective=inputdic['product']
    model_pfba_solution = cobra.flux_analysis.pfba(model)
    model_solution_frame = model_pfba_solution.to_frame()
    Optimal_rate=model_solution_frame.loc[inputdic['product'],'fluxes']
    flux_solution_df_select=model_solution_frame[abs(model_solution_frame['fluxes'])>0.0000001]#只储存大于0的值e-6
    model_pfba_solution=pd.DataFrame() 
    if flux_solution_df_select.shape[0] and Optimal_rate > 0:
        for index, row in flux_solution_df_select.iterrows():
            model_pfba_solution.loc[index,'reaction_id']=index
            model_pfba_solution.loc[index,'reaction_name']=model.reactions.get_by_id(index).name
            # model_pfba_solution.loc[index,'fluxes']=abs(round(row['fluxes'],3))
            model_pfba_solution.loc[index,'fluxes']=round(row['fluxes'],7)
            model_pfba_solution.loc[index,'gene id']=model.reactions.get_by_id(index).gene_reaction_rule
            #model_pfba_solution.loc[index,'ec-code']=norm_model.reactions.get_by_id(index).annotation['ec-code'] 
            model_pfba_solution.loc[index,'equ']=model.reactions.get_by_id(index).reaction
            model_pfba_solution.loc[index,'equ_name']=model.reactions.get_by_id(index).build_reaction_string(True).replace('O2 O2','O2').replace('CO2 CO2','CO2').replace('H2O H2O','H2O')
            model_pfba_solution.loc[index,'abs_fluxes']=abs(round(row['fluxes'],7))
    model_pfba_solution=model_pfba_solution.sort_values(by=['abs_fluxes'],ascending=False)    
    return model_pfba_solution

def get_d3_flux(model_pfba_solution,wild_flux,over_flux,fc_flux,path_results2):
    model_pfba_solutionA = model_pfba_solution.copy()
    model_pfba_solutionB = model_pfba_solution.copy()
    model_pfba_solutionC = model_pfba_solution.copy()
    # 填充 wild_flux 值到 model_pfba_solutionA 的 fluxes 列，未找到则填充 0
    model_pfba_solutionA['fluxes'] = model_pfba_solutionA.index.map(wild_flux).fillna(0)
    # 填充 over_flux 值到 model_pfba_solutionB 的 fluxes 列，未找到则填充 0
    model_pfba_solutionB['fluxes'] = model_pfba_solutionB.index.map(over_flux).fillna(0)
    # 填充 fc_flux 值到 model_pfba_solutionC 的 fluxes 列，未找到则填充 0
    model_pfba_solutionC['fluxes'] = model_pfba_solutionC.index.map(fc_flux).fillna(0)
    # 替换 abs_fluxes 列为 fluxes 的绝对值
    model_pfba_solutionA['abs_fluxes'] = model_pfba_solutionA['fluxes'].abs()
    model_pfba_solutionB['abs_fluxes'] = model_pfba_solutionB['fluxes'].abs()
    model_pfba_solutionC['abs_fluxes'] = model_pfba_solutionC['fluxes'].abs()
    wildflux_outfile=os.path.join(path_results2, 'wild_d3flux.tsv')
    model_pfba_solutionA=model_pfba_solutionA.sort_values(by=['abs_fluxes'],ascending=False)
    model_pfba_solutionA.to_csv(wildflux_outfile,index=False,sep='\t') 
    overflux_outfile=os.path.join(path_results2, 'over_d3flux.tsv')
    model_pfba_solutionB=model_pfba_solutionB.sort_values(by=['abs_fluxes'],ascending=False)
    model_pfba_solutionB.to_csv(overflux_outfile,index=False,sep='\t') 
    fcflux_outfile=os.path.join(path_results2, 'fc_d3flux.tsv')
    model_pfba_solutionC=model_pfba_solutionC.sort_values(by=['abs_fluxes'],ascending=False)
    model_pfba_solutionC.to_csv(fcflux_outfile,index=False,sep='\t')
    
    model.objective=inputdic['biomass']
    v0_biomass=model.optimize().objective_value
    model.reactions.get_by_id(inputdic['substrate']).bounds=(-inputdic['concentration'],0)
    model_pfba_solution_wt = cobra.flux_analysis.pfba(model)
    need_fluxes_wt = model_pfba_solution_wt.fluxes[abs(model_pfba_solution_wt.fluxes)>1e-6]
    model.reactions.get_by_id(inputdic['biomass']).bounds=(v0_biomass*0.1,v0_biomass*0.1)
    model.objective=inputdic['product']
    v1_product_max=model.optimize().objective_value
    model_pfba_solution_over = cobra.flux_analysis.pfba(model)
    need_fluxes_over = model_pfba_solution_over.fluxes[abs(model_pfba_solution_over.fluxes)>1e-6]
    # pfba计算--need_fluxes
    flux_table = os.path.join(path_results2, 'wild_d3flux.tsv')
    substrateId = inputdic['substrate']
    productId = inputdic['product']
    object_rxn_id = 'wild_d3flux'
    with open(flux_table, 'w') as flux: 
        for r,v in need_fluxes_wt.items():
            rxns = model.reactions.get_by_id(r) 
            try:
                check = rxns.check_mass_balance() 
                flux.write(f"{r}\t{round(v,3)}\t{rxns.reaction}\t{rxns.build_reaction_string(use_metabolite_names=True)}\t{rxns.bounds}\t{check}\n") 
            except:
                flux.write(f"{r}\t{round(v,3)}\t{rxns.reaction}\t{rxns.build_reaction_string(use_metabolite_names=True)}\t{rxns.bounds}\tno_check_mass_balance\n") 
    visualisation_html_file_WT = pathway_tsv_visualization(flux_table,substrateId,productId,cofactor_switch = "T",flow="True",outputdir=path_results2,rxnId=object_rxn_id,model=model)
    flux_table = os.path.join(path_results2, 'fc_d3flux.tsv')
    fc_filtered = model_pfba_solutionC[['reaction_id', 'fluxes']]
    substrateId = inputdic['substrate']
    productId = inputdic['product']
    object_rxn_id = 'fc_d3flux'
    with open(flux_table, 'w') as flux: 
    # 遍历 reaction_id 和 fluxes
        for _, row in fc_filtered.iterrows():
            r = row['reaction_id']  # reaction_id
            v = row['fluxes']  
            rxns = model.reactions.get_by_id(r) 
            try:
                check = rxns.check_mass_balance() 
                flux.write(f"{r}\t{round(v,3)}\t{rxns.reaction}\t{rxns.build_reaction_string(use_metabolite_names=True)}\t{rxns.bounds}\t{check}\n") 
            except:
                flux.write(f"{r}\t{round(v,3)}\t{rxns.reaction}\t{rxns.build_reaction_string(use_metabolite_names=True)}\t{rxns.bounds}\tno_check_mass_balance\n") 
    visualisation_html_file = pathway_tsv_visualization(flux_table,substrateId,productId,cofactor_switch = "T",flow="True",outputdir=path_results2,rxnId=object_rxn_id,model=model)
    flux_table = os.path.join(path_results2, 'over_d3flux.tsv')
    substrateId = inputdic['substrate']
    productId = inputdic['product']
    object_rxn_id = 'over_d3flux'
    with open(flux_table, 'w') as flux: 
        for r,v in need_fluxes_over.items():
            rxns = model.reactions.get_by_id(r) 
            try:
                check = rxns.check_mass_balance() 
                flux.write(f"{r}\t{round(v,3)}\t{rxns.reaction}\t{rxns.build_reaction_string(use_metabolite_names=True)}\t{rxns.bounds}\t{check}\n") 
            except:
                flux.write(f"{r}\t{round(v,3)}\t{rxns.reaction}\t{rxns.build_reaction_string(use_metabolite_names=True)}\t{rxns.bounds}\tno_check_mass_balance\n") 
    visualisation_html_file_fc = pathway_tsv_visualization(flux_table,substrateId,productId,cofactor_switch = "T",flow="True",outputdir=path_results2,rxnId=object_rxn_id,model=model)
    

    return model_pfba_solutionA,model_pfba_solutionB,model_pfba_solutionC,visualisation_html_file,visualisation_html_file_WT,visualisation_html_file_fc      
        

if __name__=="__main__":
    path_model=sys.argv[1]
    path_task=sys.argv[2]
    path_map=sys.argv[3]
    path_results=sys.argv[4]
    taskname=sys.argv[5]
    path_results2=os.path.join(path_results,taskname)

    if not os.path.exists(path_results2):
        os.makedirs(path_results2)

    path_results3 = os.path.join(path_results2, "iBridge")
    if not os.path.exists(path_results3):
        os.makedirs(path_results3)

    use_sequential_threads()
    install_cobra_threads()
    cobra.Configuration().solver = os.environ.get("OPTME_COBRA_SOLVER", "cplex")
    model0,model,inputdic = prepare_model(path_model,path_task,taskname)
    if inputdic['taskname'] == 'iBridge':
        target_rxn_list = [inputdic['product']]
        biomass_rxn = inputdic['biomass']
        model_file = os.path.join(path_model,inputdic['model'])
        constrict = {}
        for target_reaction in target_rxn_list:
            print('Target reaction: %s'%target_reaction)
            model_reactions = [reaction.id for reaction in model.reactions]

            corr_output = '%s/corr_%s.csv'%(path_results3,target_reaction)
            cov_output = '%s/cov_%s.csv'%(path_results3, target_reaction)

            flux_df, flux_corr_df, flux_cov_df = run_MetScore_simulation(model, biomass_rxn,
                                                                        target_reaction, constrict)
            flux_corr_df = flux_corr_df[target_reaction].loc[model_reactions]
            flux_corr_df.to_csv(corr_output)
            flux_cov_df = flux_cov_df[target_reaction].loc[model_reactions]
            flux_cov_df.to_csv(cov_output)

            write_header(corr_output, target_reaction)
            write_header(cov_output, target_reaction)

            df = pd.read_csv(cov_output, index_col=0, header=0)

            final_dic = {}
            for each_col in df.columns:
                final_dic[each_col]={}
                fluxsum_dic = calculate_MetScore_sum(model,  dict(df[each_col]))
                final_dic[each_col]=fluxsum_dic

            final_df = pd.DataFrame.from_dict(final_dic)
            final_df.to_csv('%s/MetScore_%s.csv'%(path_results3, target_reaction))

            flux_corr_df = pd.read_csv(corr_output, index_col=0)
            flux_cov_df = pd.read_csv(cov_output, index_col=0)
            output_file = '%s/Final_MetScore_%s.txt'%(path_results3, target_reaction)
            select_candidates(model, target_reaction, output_file, 
                            final_df, flux_corr_df, flux_cov_df, 0, 0.0)
            
    files = glob.glob('%s/Final_*.txt'%path_results3)
    for filename in files:
        basename = os.path.basename(filename)
        df = pd.read_table(filename)
        make_candidate_reaction_sets(df, basename, path_task, path_results3)
    
    original_model_met_set_info = metabolite_set(model)
    universal_model_met_set_info = metabolite_set_cytosol(model0)
    print('Number of metabolite set for original model: %d'\
            %len(original_model_met_set_info))
    print('Number of metabolite set in cytosol for universal model: %d'\
            %len(universal_model_met_set_info))

    files = glob.glob('%s/application_results/*.txt'%path_results3)

    cnt = 1
    for each_file in files:
        s = time.time()
        basename = os.path.basename(each_file)
        ex_reaction_id = basename.split('Final_MetScore_')[1].strip()
        ex_reaction_id = ex_reaction_id.replace('.txt', '')

        corr_file = glob.glob('%s/corr*%s*.csv'%(path_results3, ex_reaction_id))[0]
        cov_file = glob.glob('%s/cov*%s*.csv'%(path_results3, ex_reaction_id))[0]

        flux_corr_df = pd.read_csv(corr_file, index_col=0)
        flux_cov_df = pd.read_csv(cov_file, index_col=0)

        print('Analysis number: %d\tTarget reaction: %s'%(cnt, ex_reaction_id))

        check_reaction(each_file, original_model_met_set_info, universal_model_met_set_info, flux_corr_df, flux_cov_df, path_results3)
        check_down_reaction(each_file, original_model_met_set_info, universal_model_met_set_info, flux_corr_df, flux_cov_df, path_results3)
        cnt+=1
        e = time.time()
        print('Elapsed time: %fs'%(e-s))
    

    ex_metabolites = []
    with open('%s/excluded_metabolites.txt'%path_task, 'r') as fp:
        for line in fp:
            ex_metabolites.append(line.strip())

    known_reaction_info = {}
    for each_reaction in model.reactions:
        reactants, products = parse_reactions(each_reaction.reaction, ex_metabolites)
        known_reaction_info[each_reaction] = [reactants, products]
        
            

    os.makedirs(os.path.join(path_results3, "application_result2"), exist_ok=True)
    files = [item for item in os.listdir('%s/application_result2/'%path_results3) if '_up_' in item]
    fp = open('%s/Application_final_up_summary.txt'%path_results3, 'w')
    fp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'\
            %('Target', 'Gene source', 'Negative metabolite', 
            'Positive metabolite', 'Negative score', 'Positive score', 
            'reaction', 'Equation', 'Corr', 'Cov'))

    for each_file in files:
        with open('%s/application_result2/%s'%(path_results3, each_file), 'r') as fp2:
            fp2.readline()
            for line in fp2:
                sptlist = line.strip().split('\t')
                equation = sptlist[7].strip()
                met1 = sptlist[2].strip()
                met2 = sptlist[3].strip()
                if met1 not in ex_metabolites and met2 not in ex_metabolites:
                    fp.write(line.strip() + '\n')
    fp.close()


    files = [item for item in os.listdir('%s/application_result2/'%path_results3) if '_down_' in item]
    fp = open('%s/Application_final_down_summary.txt'%path_results3, 'w')
    fp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'\
            %('Target', 'Gene source', 'Negative metabolite', 
            'Positive metabolite', 'Negative score', 'Positive score', 
            'reaction', 'Equation', 'Corr', 'Cov'))

    for each_file in files:
        with open('%s/application_result2/%s'%(path_results3, each_file), 'r') as fp2:
            fp2.readline()
            for line in fp2:
                sptlist = line.strip().split('\t')
                equation = sptlist[7].strip()
                met1 = sptlist[2].strip()
                met2 = sptlist[3].strip()
                if met1 not in ex_metabolites and met2 not in ex_metabolites:
                    fp.write(line.strip() + '\n')
    fp.close()
    
    reactions_detail,range_change,up_reaction,down_reaction,df_up,df_down = detail_put(model)
    prod0,v0_biomass = calculate_prod(model, inputdic)
    Yield0,Carbon_Yield0,Ycm,Mass_Yield0,product_c_num,substrate_c_num,substrate,substrate_molecular,product_molecular = calculated_yield(path_task,inputdic,model,prod0)
    degree_pro,degree_C_pro =reduced_degree(model,inputdic['product'])
    degree_bio,degree_C_bio =reduced_degree(model,inputdic['substrate'])
    output = output_web_ibridge(reactions_detail,prod0)
    reactions_with_ranges = conbine_flux(model)
    product,over_flux_fseof =calculate_product(inputdic)
    v_biomass,wild_flux_fseof = calculate_biomass(inputdic)    
    model_pfba_solution = reaction_list_d3(model,inputdic)
    # model_pfba_solutionA,model_pfba_solutionB,model_pfba_solutionC,visualisation_html_file,visualisation_html_file_WT,visualisation_html_file_fc =  get_d3_flux(model_pfba_solution,wild_flux_fseof,over_flux_fseof,reactions_with_ranges,path_results3)
    source_file = os.path.join(path_map, 'iJO1366.Central_metabolism.json')
    destination_file = os.path.join(path_results3, 'map.json')
    shutil.copy2(source_file, destination_file)
    os.rename(destination_file, os.path.join(path_results3, 'map.json'))

