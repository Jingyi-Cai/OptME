# coding:utf-8
"""
Author: Jingyi Cai, Wenqi Xu (2024-2026)
Function: Predict metabolic engineering targets with FSEOF or loopless OptForce (MUST). The task field taskname selects which method runs.
Input: python target.py <model_dir> <task_dir> <map_dir> <results_dir> <task_id>. Reads the task JSON, the model, INPUT/map/iJO1366.Central_metabolism.json, and the carbon table. The COBRA solver is OPTME_COBRA_SOLVER, default cplex.
Output: <results_dir>/<task_id>/FSEOF/ or OptForce/, including output.json, results.xlsx, flux maps, d3flux HTML, and target figures.
"""
import sys
import cobra
import cProfile
from modelbuilder import *
import copy
from cobra import Model, Reaction, Metabolite
from sympy import subsets
from cobra.flux_analysis.variability import flux_variability_analysis
from cobra import Reaction, Metabolite, Model
from cobra.flux_analysis.loopless import add_loopless, loopless_solution
from cobra.flux_analysis import pfba
import json
import numpy as np
import os
import pandas as pd
import re
import time
import shutil
import uuid
import boto3
import re
from visualization2 import *
import matplotlib.pyplot as plt
from visualcode import *




def prepare_model_optforce(path_model,path_task,taskname):
    # model_file0="/home/sun/ETGEMS-10.20/data/iML1515_new.json"
    # model_file="/home/sun/ETGEMS-10.20/data/iML1515_irr_enz_constraint_adj_round2.json"
    with open(os.path.join(path_task,taskname)+'.json',encoding='utf-8') as fp:
        inputdic=json.load(fp)
    if  '.mat' in inputdic['model']: 
        model =cobra.io.load_matlab_model(os.path.join(path_model,inputdic['model']))
    elif '.xml' in inputdic['model']:
        model = cobra.io.read_sbml_model(os.path.join(path_model,inputdic['model']))
    elif '.sbml' in inputdic['model']:
        model = cobra.io.read_sbml_model(os.path.join(path_model,inputdic['model']))   #  need add judge
    elif '.json' in inputdic['model']:
        model = cobra.io.load_json_model(os.path.join(path_model,inputdic['model']))
    elif '' in inputdic['model']:
        model = cobra.io.load_model( inputdic['model'])
    o2_reaction_id = inputdic.get('O2', 'EX_o2_e')
    try:
        o2_reaction = model.reactions.get_by_id(o2_reaction_id)
    except KeyError:
        print(f"Warning: oxygen reaction {o2_reaction_id} not found in model {inputdic['model']}, skip oxygen bound adjustment.")
    else:
        if inputdic['oxygenstate']=='aerobic':
            o2_reaction.lower_bound = -1000
        if inputdic['oxygenstate']=='micro_aerobic': 
            o2_reaction.lower_bound = -2 
        if inputdic['oxygenstate']=='anaerobic': 
            o2_reaction.lower_bound = 0
    if inputdic['CO2_flag']=='True':
        model.reactions.get_by_id('EX_co2_e').lower_bound  = -1000
    else:
        model.reactions.get_by_id('EX_co2_e').lower_bound  = 0
    return model,inputdic      
# model0,model,inputdic =prepare_model(path_model,path_task,taskname)

def calculate_biomass_optforce(inputdic,model,path_results):
    model.objective=inputdic['biomass']
    v0_biomass=model.optimize().objective_value
    fva_wild = flux_variability_analysis(model, loopless=True,fraction_of_optimum=0.9)
        # 创建空字典
    result_dict = {}
    # model.reactions.get_by_id(inputdic['substrate']).bounds=(-inputdic['substrate_uptake_rate'],0)
    model_pfba_solution = cobra.flux_analysis.pfba(model)
    WTbiomass = model_pfba_solution.fluxes[inputdic['biomass']]
    model_solution_frame_wt = model_pfba_solution.to_frame()
    wild_flux = model_solution_frame_wt.fluxes.to_dict()
    values = []
    for key, value in wild_flux.items():
        values.append(value)
    # Calculate mean and standard deviation
    mean = np.mean(values)
    std_dev = np.std(values)
    wild_flux['range_lb']=(mean-2*std_dev).round(3)
    wild_flux['range_ub']=(mean + 2 * std_dev).round(3)  
    # 遍历 DataFrame 中的每一行
    for index, row in fva_wild.iterrows():
        result_dict[index] = {'range': [row['minimum'], row['maximum']]}

    # 将字典写入 JSON 文件
    wild_file_path = os.path.join(path_results, 'wild-enzyme.json')
    with open(wild_file_path, 'w') as file:
        json.dump(result_dict, file)

    wild_flux_path = os.path.join(path_results, 'wildflux_map.json')
    with open(wild_flux_path, 'w') as file:
        json.dump(wild_flux, file)
    return v0_biomass,fva_wild,wild_flux
# v0_biomass,fva_wild=calculate_biomass(inputdic,model,path_results)

def calculate_product_optforce(inputdic,v0_biomass,model,path_results):
    model.reactions.get_by_id(inputdic['biomass']).bounds=(v0_biomass*0.1,v0_biomass*0.1)
    model.objective=inputdic['product']
    v1_product_max=model.optimize().objective_value
    fva_over = flux_variability_analysis(model, loopless=True,fraction_of_optimum=0.9)
    # 创建空字典
    result_dict = {}
    model_pfba_solution = cobra.flux_analysis.pfba(model)
    product_max = model_pfba_solution.fluxes[inputdic['product']]
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
    # 遍历 DataFrame 中的每一行
    for index, row in fva_over.iterrows():
        result_dict[index] = {'range': [row['minimum'], row['maximum']]}
    over_file_path = os.path.join(path_results, 'over-enzyme.json')
    # 将字典写入 JSON 文件
    with open(over_file_path, 'w') as file:
        json.dump(result_dict, file)

    over_flux_path = os.path.join(path_results, 'overflux_map.json')
    with open(over_flux_path, 'w') as file:
        json.dump(over_flux, file)   
    return v1_product_max,fva_over,over_flux
# v1_product_max,fva_over = calculate_product(inputdic,v0_biomass,model,path_results)

def read_file(path_results):
    wild_file_path = os.path.join(path_results, 'wild-enzyme.json')
    over_file_path = os.path.join(path_results, 'over-enzyme.json') 
    with open(wild_file_path, 'r') as enzyme_results_file:
        enzyme_results_data = json.load(enzyme_results_file)
    with open(over_file_path, 'r') as enzyme_overresults2_file:
        enzyme_overresults2_data = json.load(enzyme_overresults2_file)
    return enzyme_results_data,enzyme_overresults2_data
# enzyme_results_data,enzyme_overresults2_data = read_file(path_results)


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
    direction = {}
    for enzyme, values in enzyme_results_data.items():
        if values.get('range') and len(values['range']) >= 2:
            range_values = values['range']
            first_value = range_values[0]

            if first_value > 1e-2:  # 满足第一个条件
                if enzyme in enzyme_overresults_data:
                    overresults_values = enzyme_overresults_data[enzyme]
                    if overresults_values.get('range') and len(overresults_values['range']) >= 2:
                        over_first_value = overresults_values['range'][0]
                        over_second_value = overresults_values['range'][1]
                        if abs(over_first_value) < 0.001 and abs(over_second_value) < 0.001:
                            over_first_value = 0
                            over_second_value = 0
                        if over_first_value==over_second_value == 0:  # 满足第二个条件
                            ko_data.append((enzyme, range_values))
                            direction[enzyme] = 0   # 添加direction正
            elif range_values[1] < -1e-2: 
                if enzyme in enzyme_overresults_data:
                    overresults_values = enzyme_overresults_data[enzyme]
                    if overresults_values.get('range') and len(overresults_values['range']) >= 2:
                        over_first_value = overresults_values['range'][0]
                        over_second_value = overresults_values['range'][1]    
                        if abs(over_first_value) < 0.001 and abs(over_second_value) < 0.001:
                            over_first_value = 0
                            over_second_value = 0
                        if over_first_value==over_second_value == 0:  # 满足第三个条件
                            ko_data.append((enzyme, range_values))
                            direction[enzyme] = 1  # 添加direction负

    # up flux
    up_data = []

    for enzyme, values in enzyme_results_data.items():
        if values.get('range') and len(values['range']) >= 2:
            range_values = values['range']
            second_value = range_values[1]
            
            if enzyme in enzyme_overresults_data:
                overresults_values = enzyme_overresults_data[enzyme]
                
                if overresults_values.get('range') and len(overresults_values['range']) >= 2:
                    over_range_values = overresults_values['range']
                    over_first_value = over_range_values[0]
                    over_second_value = over_range_values[1]
                    
                    num_zeros = sum(1 for value in range_values if value == 0)
                    if num_zeros != 3 and second_value <= over_first_value and (range_values != [0, 0] and over_range_values != [0, 0]):
                        avg_values = sum(range_values) / len(range_values)
                        avg_over_values = sum(over_range_values) / len(over_range_values)
                        
                        if avg_values > 1e-1 or avg_over_values > 1e-1:
                            up_data.append((enzyme, range_values))
                            direction[enzyme] = 0 # 添加direction正                            
# 2.上调野生型、过表达型4个值均小于0，野生型最小值大于过表达型最大值
    for enzyme, values in enzyme_results_data.items():
        if values.get('range') and len(values['range']) >= 2:
            range_values = values['range']
            second_value = range_values[1]
            
            if enzyme in enzyme_overresults2_data:
                overresults_values = enzyme_overresults2_data[enzyme]
                
                if overresults_values.get('range') and len(overresults_values['range']) >= 2:
                    over_range_values = overresults_values['range']
                    over_first_value = over_range_values[0]
                    over_second_value = over_range_values[1]
                    
                    # 检查是否所有值都小于0
                    if all(value < 0 for value in range_values) and all(value < 0 for value in over_range_values):
                        # 检查野生型最小值是否大于过表达型最大值
                        if min(range_values) > max(over_range_values):
                            # 检查是否野生型和过表达型均不包含0
                            if not (0 in range_values or 0 in over_range_values):
                                up_data.append((enzyme, range_values))
                                direction[enzyme] = 1  # 添加direction负
    for enzyme, values in enzyme_results_data.items():
        if values.get('range') and len(values['range']) >= 2:
            range_values = values['range']
            
            # 检查野生型的2个值是否一正一负
            if range_values[0] * range_values[1] < 0:
            
                if enzyme in enzyme_overresults2_data:
                    overresults_values = enzyme_overresults2_data[enzyme]
                    
                    if overresults_values.get('range') and len(overresults_values['range']) >= 2:
                        over_range_values = overresults_values['range']
                        
                        # # 检查过表达型的两个值是否与野生型的范围不重合
                        # if (min(over_range_values) > max(range_values)) or (max(over_range_values) < min(range_values)):
                        #     up_data.append((enzyme, range_values))
                        # 检查过表达型的两个值是否与野生型的范围不重合
                        if min(over_range_values) > max(range_values):
                            up_data.append((enzyme, range_values))
                            direction[enzyme] = 0  # 添加direction正
                        elif max(over_range_values) < min(range_values):
                            up_data.append((enzyme, range_values))
                            direction[enzyme] = 1  #负    

    # down_flux
    # down_flux
    down_data = []

    for enzyme, values in enzyme_results_data.items():
        if values.get('range') and len(values['range']) >= 2:
            range_values = values['range']
            first_value = range_values[0]

            if enzyme in enzyme_overresults2_data:
                overresults_values = enzyme_overresults2_data[enzyme]
                if overresults_values.get('range') and len(overresults_values['range']) >= 2:
                    over_range_values = overresults_values['range']
                    over_first_value = overresults_values['range'][0]
                    over_second_value = overresults_values['range'][1]

                    num_zeros = sum(1 for value in range_values if value == 0)
                    # ¼ì²é·¶Î§ÊÇ·ñÂú×ãÌõ¼þ£¨²»ÊÇ&#8203;``¡¾oaicite:2¡¿``&#8203;£©
                    if num_zeros != 3 and first_value >= over_second_value and first_value >= over_first_value and (range_values != [0, 0] and over_range_values != [0, 0]):
                        avg_values = sum(range_values) / len(range_values)
                        avg_over_values = sum(over_range_values) / len(over_range_values) 

                        if avg_values > 1e-1 or avg_over_values > 1e-1:                       
                            down_data.append((enzyme, range_values))
                            direction[enzyme] = 0  # 添加direction正
    for enzyme, values in enzyme_results_data.items():
        if values.get('range') and len(values['range']) >= 2:
            range_values = values['range']
            
            # 检查野生型的所有值是否均小于0
            if all(value < 0 for value in range_values):
            
                if enzyme in enzyme_overresults_data:
                    overresults_values = enzyme_overresults_data[enzyme]
                    
                    if overresults_values.get('range') and len(overresults_values['range']) >= 2:
                        over_range_values = overresults_values['range']
                        
                        # 检查过表达型的最小值是否大于野生型的最大值
                        if min(over_range_values) > max(range_values):
                            down_data.append((enzyme, range_values))
                            direction[enzyme] =1  # 添加direction负
    for enzyme, values in enzyme_results_data.items():
        if values.get('range') and len(values['range']) >= 2:
            range_values = values['range']
            
            # 检查野生型的两个值是否一个为正，一个为负
        if enzyme in enzyme_overresults_data:
            overresults_values = enzyme_overresults_data[enzyme]
            
            if overresults_values.get('range') and len(overresults_values['range']) >= 2:
                over_range_values = overresults_values['range']
                
                # 检查过表达型的两个值是否一个为正，一个为负
                if over_range_values[0] * over_range_values[1] < 0:
                
                    # # 检查野生型的范围和过表达型的范围是否不重合
                    # if (min(over_range_values) > max(range_values)) or (max(over_range_values) < min(range_values)):
                    #     down_data.append((enzyme, range_values))
                                            # 检查野生型的范围和过表达型的范围是否不重合
                    if (min(over_range_values) > max(range_values)) :
                        down_data.append((enzyme, range_values))
                        direction[enzyme] = 1  # 添加direction 负
                    elif (max(over_range_values) < min(range_values)):
                        down_data.append((enzyme, range_values))
                        direction[enzyme] = 0 # 添加direction正

    range_change = []
    mean_change = []
    for enzyme, values in enzyme_results_data.items():
        if values.get('range') and len(values['range']) >= 2:
            range_values = values['range']
            first_value = range_values[0]
            second_value = range_values[1]

            if enzyme in enzyme_overresults2_data:
                overresults_values = enzyme_overresults2_data[enzyme]
                if overresults_values.get('range') and len(overresults_values['range']) >= 2:
                    over_first_value = overresults_values['range'][0]
                    over_second_value = overresults_values['range'][1]

                    over_range = (over_first_value + over_second_value) / 2
                    wild_range = (first_value + second_value) / 2
                    
                    # Ìí¼Ó¼ì²éÌõ¼þ£¬È·±£ wild_range ²»ÎªÁã
                    if wild_range >=1e-5:
                        range_value = over_range / wild_range
                        mean_value = (over_range + wild_range) / 2
                        range_change.append((enzyme,range_value))
                        mean_change.append((enzyme,mean_value))
    
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
    return ko_data,up_data,down_data,range_change,mean_change,direction
# ko_data,up_data,down_data,range_change,mean_change,direction=compare_results(enzyme_results_data, enzyme_overresults_data)


def gene_reaction_map(reaction_list,model_input):
    equation_dict = {}

    for reaction_id in reaction_list:
        equation = model_input.reactions.get_by_id(reaction_id).reaction
        equation_dict[reaction_id] = equation
    return equation_dict
# equation_dict = gene_reaction_map(reaction_list,model_input)
def get_gpr_for_reactions_optforce(reaction_list, model_input):
    gpr_dict = {}
    for reaction_id in reaction_list:
        gpr = model_input.reactions.get_by_id(reaction_id).gene_reaction_rule
        gpr_dict[reaction_id] = gpr

    return gpr_dict
# output
# output
# output
def must_df(enzyme_results_data,enzyme_overresults2_data,range_change,ko_data,up_data,down_data):
    wild_data = [{'reaction': reaction, 'flux_wild': [format(value, '.3f') for value in data['range']]} for reaction, data in enzyme_results_data.items()]
    df1 = pd.DataFrame(wild_data)
    over_data = [{'reaction': reaction, 'flux_over': [format(value, '.3f') for value in data['range']]} for reaction, data in enzyme_overresults2_data.items()]
    df2 = pd.DataFrame(over_data)
    meged_df = pd.merge(df1, df2, on='reaction', how='inner')
    meged_df['equation'] = meged_df['reaction'].map(equation_dict)
    meged_df['gpr'] = meged_df['reaction'].map(gpr_dict)
    # ?? 'manipulations' Â¨Â¢D3?Â¨Âº??Â¡Â¥?a None
    meged_df['manipulations'] = None
    for reaction, _ in ko_data:
        meged_df.loc[meged_df['reaction'] == reaction, 'manipulations'] = 'ko'
    for reaction, _ in up_data:
        meged_df.loc[meged_df['reaction'] == reaction, 'manipulations'] = 'Up'
    for reaction, _ in down_data:
        meged_df.loc[meged_df['reaction'] == reaction, 'manipulations'] = 'down'
    range_change_data = [{'reaction': enzyme, 'range_change': value} for enzyme,value in range_change]
    mean_change_data = [{'reaction': enzyme, 'mean_change': value} for enzyme,value in mean_change]
    # 更新 range_change 的值
    # 将 range_change_data 的 'range_change' 值处理为空或为 'full' 的情况
    for data in range_change_data:
        if data['range_change'] in [None, 'full']:
            data['range_change'] = 0
        elif isinstance(data['range_change'], str) and not data['range_change'].replace('.', '', 1).isdigit():
            # 如果值是其他非数值类型字符串，直接设置为 0
            data['range_change'] = 0
    for data in range_change_data:
        value = data['range_change']
        if value > 0:
            data['range_change'] = round(np.log2(value), 3)
        else:
            data['range_change'] = 0
    range_change_df = pd.DataFrame(range_change_data)

    mean_change_df = pd.DataFrame(mean_change_data)
    # ºÏ²¢µ½ megad_df
    meged_df = pd.merge(meged_df, range_change_df, on='reaction', how='inner')
    meged_df = pd.merge(meged_df, mean_change_df, on='reaction', how='inner')
    meged_df['range_change'] = meged_df['range_change'].apply(lambda x: round(x, 3))
    meged_df['mean_change'] = meged_df['mean_change'].apply(lambda x: round(x, 3))
    empty_gpr_indices = meged_df[meged_df['gpr'] == ''].index
    meged_df.loc[(meged_df['manipulations'] == 'Up') & (meged_df['range_change'] <= 1), 'manipulations'] = None
    meged_df.loc[(meged_df['manipulations'] == 'down') & (meged_df['range_change'] >= -1), 'manipulations'] = None
# 将这些行的 'manipulation' 列改为 None
    meged_df.loc[empty_gpr_indices, 'manipulations'] = None
    transport_rxns = []
    exchange_rxns = []
    for rxn in model.reactions:
        reactants_mets = [str(m).rsplit('_', 1)[0] for m in rxn.reactants]
        products_mets = [str(m).rsplit('_', 1)[0] for m in rxn.products]
        
        if sorted(reactants_mets) == sorted(products_mets):
            transport_rxns.append(rxn.id)
        if len(products_mets) == 0:
            exchange_rxns.append(rxn.id)
    transport_rxns_pattern = '|'.join(transport_rxns)
    # 检查 transport_rxns_pattern 是否为有效正则表达式
    try:
        re.compile(transport_rxns_pattern)
        is_valid_regex = True
    except re.error:
        is_valid_regex = False

    if is_valid_regex:
        # 如果 transport_rxns_pattern 是有效的正则表达式，执行过滤
        meged_df.loc[
            meged_df['reaction'].str.contains(transport_rxns_pattern, regex=True, na=False),
            'manipulations'
        ] = None
    else:
        # 如果无效，则跳过该步骤
        print("Warning: Invalid regex pattern in transport_rxns_pattern, skipping this step.")
    # meged_df.loc[meged_df['reaction'].str.contains(transport_rxns_pattern), 'manipulations'] = None
    return meged_df
# meged_df = must_df(enzyme_results_data,enzyme_overresults2_data,range_change,ko_data,up_data,down_data)
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

def calculated_yield_optforce(path_task,inputdic,model,prod0):
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
    # 处理 met.charge 为 None 的情况
    if met.charge is None:
        met.charge = 0
    degree -= met.charge

    # 避免除零错误
    if 'C' not in met.elements or met.elements['C'] == 0:
        raise ValueError("The metabolite does not contain any carbon atoms.")

    degree_C = degree / met.elements['C']
    return degree,degree_C


def detail_put_optforce(model_input,meged_df,path_results4):
    reactions = {}
    ko_df = meged_df.loc[meged_df['manipulations'] == 'ko']
    ko_df = ko_df.sort_values(by='range_change', ascending=False)
    ko_reaction=ko_df['reaction'].tolist()
    for i in ko_reaction:
        try:    
            reaction = model_input.reactions.get_by_id(i)
            ko_reaction_name = reaction.name
            ko_reaction_equation = reaction.reaction
            gpr=model_input.reactions.get_by_id(i).gpr.to_string()
            range_change_value = ko_df.loc[ko_df['reaction'] == i, 'range_change'].values.tolist()
            range_change_value = range_change_value[0] if range_change_value else None
            mean_value = ko_df.loc[ko_df['reaction'] == i,'mean_change'].values.tolist()
            mean_value = round(mean_value[0], 3) if mean_value else None
            if isinstance(mean_value, list):
                mean_value = round(mean_value[0], 3) if mean_value else None
        # ½«Ã¿¸ö·´Ó¦ºÍ¶ÔÓ¦µÄ 'range_change' Öµ´æ´¢ÔÚ×ÖµäÖÐ

            reaction_dict = {
                    "Modification": "KO",
                    "Fold_Change":range_change_value ,
                    "Average_Flux":mean_value,
                    "Reaction_ID": i,
                    "Enzyme_Name": ko_reaction_name,
                    "Gene_ID": gpr,
                    "Equation": ko_reaction_equation,
                    "figpath":os.path.join(path_results4,'target_fig',i+'.png')
                }
            reactions[i] = reaction_dict
        except KeyError:
            pass
    up_df = meged_df.loc[meged_df['manipulations'] == 'Up']
    up_df = up_df.sort_values(by='range_change', ascending=False)
    up_reaction=up_df['reaction'].tolist()
    for i in up_reaction:
        try:    
            reaction = model_input.reactions.get_by_id(i)
            up_reaction_name = reaction.name
            up_reaction_equation = reaction.reaction
            gpr=model_input.reactions.get_by_id(i).gpr.to_string()
            range_change_value = up_df.loc[up_df['reaction'] == i, 'range_change'].values.tolist()
            range_change_value = range_change_value[0] if range_change_value else None
            mean_value = up_df.loc[up_df['reaction'] == i,'mean_change'].values.tolist()
            mean_value = round(mean_value[0], 3) if mean_value else None
            if isinstance(mean_value, list):
                mean_value = round(mean_value[0], 3) if mean_value else None
        # ½«Ã¿¸ö·´Ó¦ºÍ¶ÔÓ¦µÄ 'range_change' Öµ´æ´¢ÔÚ×ÖµäÖÐ

            reaction_dict = {
                    "Modification": "UP",
                    "Fold_Change":range_change_value ,
                    "Average_Flux":mean_value,
                    "Reaction_ID": i,
                    "Enzyme_Name": up_reaction_name,
                    "Gene_ID": gpr,
                    "Equation": up_reaction_equation,
                    "figpath":os.path.join(path_results4,'target_fig',i+'.png')
                }
            reactions[i] = reaction_dict
        except KeyError:
            pass
    
    down_df = meged_df.loc[meged_df['manipulations'] == 'down']
    down_df = down_df.sort_values(by='range_change', ascending=False)
    down_reaction=down_df['reaction'].tolist()
    for i in down_reaction:
            try:
                    reaction = model_input.reactions.get_by_id(i)
                    down_reaction_name = reaction.name
                    down_reaction_equation = reaction.reaction
                    gpr=model_input.reactions.get_by_id(i).gpr.to_string()
                    range_change_value = down_df.loc[down_df['reaction'] == i, 'range_change'].values.tolist()
                    range_change_value = range_change_value[0] if range_change_value else None
                    mean_value = down_df.loc[down_df['reaction'] == i,['mean_change']].values.tolist()
                    mean_value = [round(value[0], 3) for value in mean_value] if mean_value else None
                    # mean_value = down_df.loc[down_df['reaction'] == i, 'mean_change'].iloc[0]
                    # mean_value = round(mean_value[0], 3) if mean_value else None
                    if isinstance(mean_value, list):
                        mean_value = round(mean_value[0], 3) if mean_value else None
                    reaction_dict = {
                        "Modification": "Down",
                        "Fold_Change":range_change_value ,
                        "Average_Flux":mean_value,
                        "Reaction_ID": i,
                        "Enzyme_Name": down_reaction_name,
                        "Gene_ID": gpr,
                        "Equation": down_reaction_equation,
                        "figpath":os.path.join(path_results4,'target_fig',i+'.png')
                    }
                    reactions[i] = reaction_dict
            except KeyError:
                    pass
    return reactions,up_reaction,down_reaction,ko_reaction
def output_web_optforce(reactions):
    output={}
    output['summary'] = {}
    output['reaction'] = {}
    ko_number=len(ko_reaction)
    up_number=len(up_reaction)
    down_number=len(down_reaction)
    summary = {
        'UP':up_number,
        'Down':down_number,
        'KO':ko_number,
    }
    output['summary'] = summary
    ko_number=len(ko_reaction)
    up_number=len(up_reaction)
    down_number=len(down_reaction)
    output['reaction']=reactions
    output_file = os.path.join(path_results4, 'output.json')
    with open(output_file, 'w') as json_file:
        json.dump(output, json_file) 
    
    filename=os.path.join(path_results4, 'results.xlsx') # Ìæ»»ÎªÄãÏë±£´æµÄÎÄ¼þÂ·¾¶
    with pd.ExcelWriter(filename) as writer:
        meged_df.to_excel(writer, sheet_name='test',index=True)
    return output

def ko_picture_optforce(meged_df):
    ko_df = meged_df.loc[meged_df['manipulations'] == 'ko']
    first_10_columns = ko_df.iloc[:, :4]
    result_dict = {}

    for index, row in first_10_columns.iterrows():
        reaction_name = row['reaction']
        equation = row['equation']
        flux_wild_values = [float(value) for value in row['flux_wild']] 
        flux_over_values = [float(value) for value in row['flux_over']]

        temp_dict = {"reaction": equation, "wild_range": flux_wild_values, "eng_range": flux_over_values}

        result_dict[reaction_name] = temp_dict

    # 将字典转换为 JSON 格式的字符串
    result_json = json.dumps(result_dict, indent=2)

    # 输出 JSON 字符串
    ko_file = os.path.join(path_results4, 'ko.json')
    with open(ko_file, 'w') as json_file:
        json_file.write(result_json)

    return result_json

def up_picture_optforce(meged_df):
    up_df = meged_df.loc[meged_df['manipulations'] == 'Up']
    first_10_columns = up_df.iloc[:, :4]
    result_dict = {}

    for index, row in first_10_columns.iterrows():
        reaction_name = row['reaction']
        equation = row['equation']
        flux_wild_values = [float(value) for value in row['flux_wild']] 
        flux_over_values = [float(value) for value in row['flux_over']]

        temp_dict = {"reaction": equation, "wild_range": flux_wild_values, "eng_range": flux_over_values}

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
        reaction_name = row['reaction']
        equation = row['equation']
        flux_wild_values = [float(value) for value in row['flux_wild']] 
        flux_over_values = [float(value) for value in row['flux_over']]

        temp_dict = {"reaction": equation}
        
        # 将 flux_wild 和 flux_over 合并为列表形式
        temp_dict = {"reaction": equation, "wild_range": flux_wild_values, "eng_range": flux_over_values}

        result_dict_down[reaction_name] = temp_dict

    # 将字典转换为 JSON 格式的字符串
    result_json = json.dumps(result_dict_down, indent=2)

    # 输出 JSON 字符串
    down_file = os.path.join(path_results4, 'down.json')
    with open(down_file, 'w') as json_file:
        json_file.write(result_json)

    return result_json



def conbine_flux_optforce(modelid):   
    # 获取 stoichiometry flux
    # 合并 up_reaction 和 down_reaction
    all_reactions = up_reaction + down_reaction
    # 找到所有反应的范围值
    reactions_with_ranges = {}
    result_reaction_columns = list(zip(meged_df['reaction'], meged_df['range_change']))
    
    for reaction, range_value in result_reaction_columns:
        if reaction in up_reaction:
            # 处理 range_value 的取 log10
            if range_value > 0:  # log10 只适用于正值
                range_value = np.log10(range_value)
            else:
                range_value = 0  # 对于非正值，可以设为 0 或其他默认值
            # 将处理后的值放入字典
            reactions_with_ranges[reaction] = range_value
    for reaction, range_value in result_reaction_columns:
        if reaction in down_reaction:
            # 处理 range_value 的取 log10
            if range_value > 0:  # log10 只适用于正值
                range_value = -np.log10(range_value)
            else:
                range_value = 0  # 对于非正值，可以设为 0 或其他默认值
            # 将处理后的值放入字典
            reactions_with_ranges[reaction] = range_value
    flux_file = os.path.join(path_results4, 'fcflux_map.json')
    values = []
    for key, value in reactions_with_ranges.items():
        values.append(value)
    
    # 计算均值和标准差
    mean = np.mean(values)
    std_dev = np.std(values)
    
    # 计算下限和上限
    reactions_with_ranges['range_lb'] = round(mean - 2 * std_dev, 3)
    reactions_with_ranges['range_ub'] = round(mean + 2 * std_dev, 3)
    
    # 将结果写入 JSON 文件
    with open(flux_file, 'w') as json_file:
        json.dump(reactions_with_ranges, json_file) 
    
    # 返回结果
    return reactions_with_ranges


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
    # model.reactions.get_by_id(inputdic['substrate']).bounds=(-inputdic['substrate_uptake_rate'],0)
    model_pfba_solution_wt = cobra.flux_analysis.pfba(model)
    need_fluxes_wt = model_pfba_solution_wt.fluxes[abs(model_pfba_solution_wt.fluxes)>1e-1]
    model.reactions.get_by_id(inputdic['biomass']).bounds=(v0_biomass*0.1,v0_biomass*0.1)
    model.objective=inputdic['product']
    v1_product_max=model.optimize().objective_value
    model_pfba_solution_over = cobra.flux_analysis.pfba(model)
    need_fluxes_over = model_pfba_solution_over.fluxes[abs(model_pfba_solution_over.fluxes)>1e-1]
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
    visualisation_html_file = pathway_tsv_visualization(flux_table,substrateId,productId,cofactor_switch = "T",flow="True",outputdir=path_results2,rxnId=object_rxn_id,model=model)
    
    return model_pfba_solutionA,model_pfba_solutionB,model_pfba_solutionC,visualisation_html_file,visualisation_html_file_WT
        


def prepare_model(path_model,path_task,taskname):
    # model_file0="/home/sun/ETGEMS-10.20/data/iML1515_new.json"
    # model_file="/home/sun/ETGEMS-10.20/data/iML1515_irr_enz_constraint_adj_round2.json"
    with open(os.path.join(path_task,taskname)+'.json',encoding='utf-8') as fp:
        inputdic=json.load(fp)
    if  '.mat' in inputdic['model']: 
        model =cobra.io.load_matlab_model(os.path.join(path_model,inputdic['model']))
    elif '.xml' in inputdic['model']:
        model = cobra.io.read_sbml_model(os.path.join(path_model,inputdic['model']))
    elif '.sbml' in inputdic['model']:
        model = cobra.io.read_sbml_model(os.path.join(path_model,inputdic['model']))   #  need add judge
    elif '.json' in inputdic['model']:
        model = cobra.io.load_json_model(os.path.join(path_model,inputdic['model']))
    elif '' in inputdic['model']:
        model = cobra.io.load_model( inputdic['model'])
    o2_reaction_id = inputdic.get('O2', 'EX_o2_e')
    try:
        o2_reaction = model.reactions.get_by_id(o2_reaction_id)
    except KeyError:
        print(f"Warning: oxygen reaction {o2_reaction_id} not found in model {inputdic['model']}, skip oxygen bound adjustment.")
    else:
        if inputdic['oxygenstate']=='aerobic':
            o2_reaction.lower_bound = -1000
        if inputdic['oxygenstate']=='micro_aerobic': 
            o2_reaction.lower_bound = -2 
        if inputdic['oxygenstate']=='anaerobic': 
            o2_reaction.lower_bound = 0
    if inputdic['CO2_flag']=='True':
        model.reactions.get_by_id('EX_co2_e').lower_bound  = -1000
    else:
        model.reactions.get_by_id('EX_co2_e').lower_bound  = 0
    return model,inputdic     

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

def biomass(v1_product_max,inputdic):    # 以最大biomass为目标，计算每份所有反应的通量分布情况
    exlist = list(np.linspace(0,v1_product_max,10))
    exlistn = []
    for i in exlist:
        i = format(i,'.2f')
        exlistn.append(float(i))
    reactiondf = pd.DataFrame()
    for i in exlistn:
        cond = i
        model.reactions.get_by_id(inputdic['product']).bounds=(cond*0.99,cond*0.99)
        model.objective=inputdic['biomass']
        BIOMASS=model.optimize().objective_value
        BIOMASSFLUX = model.optimize()
        loopless = loopless_solution(model)
        temp_df = pd.DataFrame({'reaction': loopless.fluxes.index, 'product = '+str(cond): loopless.fluxes.values})
        if reactiondf.empty:
            reactiondf = temp_df
        else:
            reactiondf = pd.merge(reactiondf, temp_df, on='reaction', how='outer')
  
    return reactiondf


def reaction_list_d3(model,inputdic):
    # model.reactions.get_by_id(inputdic['substrate']).bounds=(-inputdic['substrate_uptake_rate'],0)
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
    # model_pfba_solution=model_pfba_solution.sort_values(by=['abs_fluxes'],ascending=False)    
    output_path =  os.path.join(path_results3, 'pfba.csv')
    model_pfba_solution.to_csv(output_path, index=False)
    return model_pfba_solution
# model_pfba_solution = reaction_list_d3(model,inputdic)

def check_monotonicity(row):
    row_values = row[2:11]  # 获取第3到第12列的数据
    increasing = all(row_values[i] < row_values[i + 1] for i in range(len(row_values) - 1))
    decreasing = all(row_values[i] > row_values[i + 1] for i in range(len(row_values) - 1))
    all_negative = all(value < 0 for value in row_values)
    
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
        # 对于所有的 `manipulation` 列为 'up' 或 'down' 且所有值均为负数的反应，将 direction 的值设为 1
    direction = {}
    for index, row in reactiondf.iterrows():
        row_values = row[2:11]  # 获取第3到第12列的数据
        has_positive = any(value > 0 for value in row_values)  # 判断是否存在正值
        has_negative = any(value < 0 for value in row_values)  # 判断是否存在负值
        if has_positive and has_negative:
            if row['manipulation'] in ['up', 'down']:  # 检查 manipulation 列是否为 'up' 或 'down'
                direction[row['reaction']] = 1
        elif has_negative:
            direction[row['reaction']] = -1
        elif has_positive:
            direction[row['reaction']] = 2
    reactiondf['direction'] = reactiondf.apply(
        lambda row: 1 if any(value > 0 for value in row[2:11]) and any(value < 0 for value in row[2:11]) and row['manipulation'] in ['up', 'down'] else 0,axis=1)    
    selected_columns = reactiondf.iloc[:, [2, 11]]
    mean_flux = selected_columns.mean(axis=1, numeric_only=True)
    reactiondf['mean_flux'] = mean_flux
    e_2_threshold = 1e-1
    reactiondf.loc[abs(mean_flux)  <= e_2_threshold, 'manipulation'] = None
    up_reactions = reactiondf[reactiondf['manipulation'] == 'up']
    down_reactions = reactiondf[reactiondf['manipulation'] == 'down']
    # product_flux = reactiondf.iloc[:, 1:11]
    # # 检查并处理数据类型
    # product_flux = product_flux.apply(pd.to_numeric, errors='coerce')
    product_flux = reactiondf.iloc[:, 1:11]
    product_flux = product_flux.apply(pd.to_numeric, errors='coerce')
    # 计算最后三列和前三列均值，使用 log2(last3/first3) 作为 Fold_Change
    last_three_cols_mean = product_flux.iloc[:, -3:].mean(axis=1)
    first_three_cols_mean = product_flux.iloc[:, :3].mean(axis=1)

    ratio = last_three_cols_mean / first_three_cols_mean
    ratio = ratio.replace([np.inf, -np.inf], np.nan)
    ratio = ratio.where(ratio > 0)

    reactiondf['result'] = np.log2(ratio)
    reactiondf['result'] = reactiondf['result'].replace([np.inf, -np.inf], np.nan).fillna(0).round(3)

    # Fold_Change threshold rule:
    # Up targets must satisfy log2FC > 1, Down targets must satisfy log2FC < -1
    reactiondf.loc[(reactiondf['manipulation'] == 'up') & (reactiondf['result'] <= 1), 'manipulation'] = None
    reactiondf.loc[(reactiondf['manipulation'] == 'down') & (reactiondf['result'] >= -1), 'manipulation'] = None
    result_per_row = {}
    for index, row in reactiondf.iterrows():
        result_per_row[row['reaction']] = row['result']

        # 找到 GPR 列为空字符串的行索引
    empty_gpr_indices = reactiondf[reactiondf['gpr'] == ''].index
# 将这些行的 'manipulation' 列改为 None
    reactiondf.loc[empty_gpr_indices, 'manipulation'] = None
    transport_rxns = []
    exchange_rxns = []
    for rxn in model.reactions:
        reactants_mets = [str(m).rsplit('_', 1)[0] for m in rxn.reactants]
        products_mets = [str(m).rsplit('_', 1)[0] for m in rxn.products]
        
        if sorted(reactants_mets) == sorted(products_mets):
            transport_rxns.append(rxn.id)
        if len(products_mets) == 0:
            exchange_rxns.append(rxn.id)
    transport_rxns_pattern = '|'.join(transport_rxns)
    # 检查 transport_rxns_pattern 是否为有效正则表达式
    try:
        re.compile(transport_rxns_pattern)
        is_valid_regex = True
    except re.error:
        is_valid_regex = False

    if is_valid_regex:
        # 如果 transport_rxns_pattern 是有效的正则表达式，执行过滤
        reactiondf.loc[
            reactiondf['reaction'].str.contains(transport_rxns_pattern, regex=True, na=False),
            'manipulation'
        ] = None
    else:
        # 如果无效，则跳过该步骤
        print("Warning: Invalid regex pattern in transport_rxns_pattern, skipping this step.")
    # reactiondf.loc[reactiondf['reaction'].str.contains(transport_rxns_pattern), 'manipulation'] = None
    # 将 direction 为 1 的行的 manipulation 列设置为 None
    reactiondf.loc[reactiondf['direction'] == 1, 'manipulation'] = None

    # Cap up and down targets to MAX_TARGETS each, ranked by abs(mean_flux) * result
    MAX_TARGETS = 30
    reactiondf['_rank_score'] = abs(reactiondf['mean_flux']) * reactiondf['result']
    for _direction in ['up', 'down']:
        _mask = reactiondf['manipulation'] == _direction
        _candidates = reactiondf[_mask].sort_values('_rank_score', ascending=False)
        if len(_candidates) > MAX_TARGETS:
            _drop_indices = _candidates.index[MAX_TARGETS:]
            reactiondf.loc[_drop_indices, 'manipulation'] = None
    reactiondf.drop(columns=['_rank_score'], inplace=True)
    up_reactions = reactiondf[reactiondf['manipulation'] == 'up']
    down_reactions = reactiondf[reactiondf['manipulation'] == 'down']

    return up_reactions,down_reactions,result_per_row,reactiondf

def gene_reaction_map(reaction_list,model_input):
    equation_dict = {}

    for reaction_id in reaction_list:
        equation = model_input.reactions.get_by_id(reaction_id).reaction
        equation_dict[reaction_id] = equation
    return equation_dict

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



def detail_put(model_input,meged_df,path_results3):
    reactions = {}
    up_df = meged_df.loc[meged_df['manipulation'] == 'up']
    up_df = up_df.sort_values(by='result', ascending=False)
    up_reaction=up_df['reaction'].tolist()
    for i in up_reaction:
        try:    
            reaction = model_input.reactions.get_by_id(i)
            up_reaction_name = reaction.name
            up_reaction_equation = reaction.reaction
            gpr=model_input.reactions.get_by_id(i).gpr.to_string()
            range_change_value = up_df.loc[up_df['reaction'] == i, 'result'].values.tolist()
            range_change_value = range_change_value[0] if range_change_value else None
            mean_value = up_df.loc[up_df['reaction'] == i,'mean_flux'].values.tolist()
            mean_value = round(mean_value[0], 3) if mean_value else None
            # mean_value = up_df.loc[up_df['reaction'] == i, 'mean_change'].iloc[0]
            # mean_value = round(mean_value[0], 3) if mean_value else None
            if isinstance(mean_value, list):
                mean_value = round(mean_value[0], 3) if mean_value else None
        # ½«Ã¿¸ö·´Ó¦ºÍ¶ÔÓ¦µÄ 'range_change' Öµ´æ´¢ÔÚ×ÖµäÖÐ

            reaction_dict = {
                    "Modification": "UP",
                    "Fold_Change":range_change_value ,
                    "Average_Flux":mean_value,
                    "Reaction_ID": i,
                    "Enzyme_Name": up_reaction_name,
                    "Gene_ID": gpr,
                    "Equation": up_reaction_equation,
                    "figpath":os.path.join(path_results3,'target_fig',i+'.png')
                }
            reactions[i] = reaction_dict
        except KeyError:
            pass
    
    down_df = meged_df.loc[meged_df['manipulation'] == 'down']
    down_df = down_df.sort_values(by='result', ascending=True)
    down_reaction=down_df['reaction'].tolist()
    for i in down_reaction:
            try:
                    reaction = model_input.reactions.get_by_id(i)
                    down_reaction_name = reaction.name
                    down_reaction_equation = reaction.reaction
                    gpr=model_input.reactions.get_by_id(i).gpr.to_string()
                    range_change_value = down_df.loc[down_df['reaction'] == i, 'result'].values.tolist()
                    range_change_value = range_change_value[0] if range_change_value else None
                    mean_value = down_df.loc[down_df['reaction'] == i,['mean_flux']].values.tolist()
                    mean_value = [round(value[0], 3) for value in mean_value] if mean_value else None
                    # mean_value = down_df.loc[down_df['reaction'] == i, 'mean_change'].iloc[0]
                    # mean_value = round(mean_value[0], 3) if mean_value else None
                    if isinstance(mean_value, list):
                        mean_value = round(mean_value[0], 3) if mean_value else None                    

                    reaction_dict = {
                        "Modification": "Down",
                        "Fold_Change":range_change_value ,
                        "Average_Flux":mean_value,
                        "Reaction_ID": i,
                        "Enzyme_Name": down_reaction_name,
                        "Gene_ID": gpr,
                        "Equation": down_reaction_equation,
                        "figpath":os.path.join(path_results3,'target_fig',i+'.png')
                    }
                    reactions[i] = reaction_dict
            except KeyError:
                    pass
    return reactions,up_reaction,down_reaction

def up_picture(reactiondf):
    # 确保数据中只包含 'up' 类型的操作
    up_df = reactiondf.loc[reactiondf['manipulation'] == 'up']
    
    # 将 DataFrame 转换为嵌套字典
    result_dict = {}

    for _, row in up_df.iterrows():
        reaction_id = row["reaction"].strip()
        
        # 使用 PySCeS 获取反应方程式
        reaction_equation = model.reactions.get_by_id(reaction_id).reaction
        
        reaction_dict = {"reaction": reaction_equation}
        
        
        # 提取 'x' 列的值
        exlist = list(np.linspace(0,product,10))
        exlist = [round(element, 3) for element in exlist]
        reaction_dict['x'] = exlist
        # 提取 'y' 列的值
        y_values = [round(float(row[col]),3) for col in up_df.columns[1:11]]
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
    down_df = reactiondf.loc[reactiondf['manipulation'] == 'down']
    # 提取前 10 列
    first_10_columns = down_df.iloc[:, :11]

    # 将 DataFrame 转换为嵌套字典
    result_dict_down = {}

    for _, row in first_10_columns.iterrows():
        reaction_id = row["reaction"].strip()
        
        # 使用 PySCeS 获取反应方程式
        reaction_equation = model.reactions.get_by_id(reaction_id).reaction
        
        reaction_dict = {"reaction": reaction_equation}
                # 提取 'y' 列的值
        # 提取 'x' 列的值
        exlist = list(np.linspace(0,product,10))
        exlist = [round(element, 3) for element in exlist]
        reaction_dict['x'] = exlist
        y_values = [round(float(row[col]),3) for col in down_df.columns[1:11]]
        reaction_dict['y'] = y_values
        
        # for col in first_10_columns.columns[1:]:
        #     product_name = col.split("=")[1].strip()
        #     reaction_dict[product_name] = row[col]
        
        result_dict_down[reaction_id] = reaction_dict

    # # 将字典保存为 JSON 文件
    down_file = os.path.join(path_results3, 'down.json')
    with open(down_file, 'w') as json_file:
        json.dump(result_dict_down, json_file, indent=2)
    
    return result_dict_down

def output_web(reactions):
    output={}
    output['summary'] = {}
    output['reaction'] = {}
    up_number=len(up_reaction)
    down_number=len(down_reaction)
    summary = {
        'UP':up_number,
        'Down':down_number,
        'KO':0,
    }
    output['summary'] = summary
    output['reaction']=reactions
    output_file = os.path.join(path_results3, 'output.json')


    
    with open(output_file, 'w') as json_file:
        json.dump(output, json_file) 
    filename=os.path.join(path_results3, 'results.xlsx') # Ìæ»»ÎªÄãÏë±£´æµÄÎÄ¼þÂ·¾¶
    with pd.ExcelWriter(filename) as writer:
        reactiondf.to_excel(writer, sheet_name='reaction',index=True)
    return output


def conbine_flux(modelid):   
    # 获取 stoichiometry flux
    # 合并 up_reaction 和 down_reaction
    all_reactions = up_reaction + down_reaction
    # 找到所有反应的范围值
    reactions_with_ranges = {}
    result_reaction_columns = list(zip(reactiondf['reaction'], reactiondf['result']))
    for reaction, range_value in result_reaction_columns:
        if reaction in up_reaction:
            # 处理 range_value 的取 log10
            if range_value > 0:  # log10 只适用于正值
                range_value = np.log10(range_value)
            else:
                range_value = 0  # 对于非正值，可以设为 0 或其他默认值
            # 将处理后的值放入字典
            reactions_with_ranges[reaction] = range_value
    for reaction, range_value in result_reaction_columns:
        if reaction in down_reaction:
            # 处理 range_value 的取 log10
            if range_value > 0:  # log10 只适用于正值
                range_value = -np.log10(range_value)
            else:
                range_value = 0  # 对于非正值，可以设为 0 或其他默认值
            # 将处理后的值放入字典
            reactions_with_ranges[reaction] = range_value    
    # for reaction, range_value in result_reaction_columns:
    #     if reaction in all_reactions:
    #         # 处理 range_value 的取 log10
    #         if range_value > 0:  # log10 只适用于正值
    #             range_value = np.log10(range_value)
    #         else:
    #             range_value = 0  # 对于非正值，可以设为 0 或其他默认值
            
    #         # 将处理后的值放入字典
            reactions_with_ranges[reaction] = range_value

    flux_file = os.path.join(path_results3, 'fcflux_map.json')

    values = []
    for key, value in reactions_with_ranges.items():
        values.append(value)
    
    # 计算均值和标准差
    mean = np.mean(values)
    std_dev = np.std(values)
    
    # 计算下限和上限
    reactions_with_ranges['range_lb'] = round(mean - 2 * std_dev, 3)
    reactions_with_ranges['range_ub'] = round(mean + 2 * std_dev, 3)
    
    # 将结果写入 JSON 文件
    with open(flux_file, 'w') as json_file:
        json.dump(reactions_with_ranges, json_file) 
    
    # 返回结果
    return reactions_with_ranges


def drawtarget(data,mode,savepath='./'):
    for id in data:
        if mode=='FSEOF':
            #  提取数据
            gludy_data = data[id]
            reaction = gludy_data['reaction']
            x_values = gludy_data['x']
            y_values = gludy_data['y']
            # 绘制图形
            plt.plot(x_values, y_values, marker='o', linestyle='-')
            plt.xlabel('Product Flux mmol/gDWh')
            plt.ylabel('reaction flux mmol/gDWh')
            plt.title(id+":"+reaction)
            plt.grid(False)
            #plt.show()
        elif mode=='OptForce':
            gludy_data = data[id]
            reaction = gludy_data['reaction']
            # 定义范围的名称和宽度
            ranges = {'Reference': gludy_data['wild_range'], 'Engineered': gludy_data['eng_range']}
            # 画图
            fig, ax = plt.subplots()
            # 绘制范围
            for i, (name, (start, end)) in enumerate(ranges.items()):
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
            plt.title(id+":"+reaction)
            plt.grid(True)
            #plt.xlabel('Product Flux mmol/gDWh')
            plt.ylabel('reaction flux mmol/gDWh')
            #plt.grid(False)
        # save the fig to the file
        plt.savefig(os.path.join(savepath,id+'.png'))
        # clear the fig
        plt.clf() 


def fc_d3fulx(reactiondf):
    #modelname=taskname.split('_')[-1]
    bigg_models_metabolites_file = os.path.join(path_task,'bigg_models_metabolites.txt')
    met_C_file = os.path.join(path_task,'met_contain_C_df.tsv')
    substrate_id = inputdic['substrate']
    product_id = inputdic['product']
    reactiondf.to_csv(os.path.join(path_results3,'reaction_df.csv'))
    # substrate_id2=[substrate_id.replace('SK_','')]
    # product_id2=[product_id.replace('user_add_DM_','')]
    # load model
    # model_substrate = pd.read_csv(model_substrate_file,sep='\t',index_col='substrate_model_correspondind')
    # model_product = pd.read_csv(model_product_file, sep='\t', index_col="obj_model_correspondind")  
    met_C = pd.read_csv(met_C_file,sep='\t', index_col=0)
    #biomass=get_initial_obj_info( model)
    rxn_switch="T"
    coordinates =""
    cofactor_switch="F"
    if'FSEOF' in inputdic['taskname']:
        # show_factor_file2 = os.path.join(path_results3, 'fc_d3flux2.html')
        # show_factor_file1 = os.path.join(path_results3, 'fc_d3flux1.html')
        show_factor_file = os.path.join(path_results3, 'fc_d3flux.html')
                # 提取 reaction, result, manipulation 列
        meged_df_extracted = reactiondf[['reaction', 'result', 'manipulation']]

        # 根据 manipulation 列调整 result 列的值
        meged_df_extracted['result'] = meged_df_extracted.apply(
            lambda row: -row['result'] if row['manipulation'] == 'down' else row['result'], axis=1
        )

        # 生成新表格 fc
        fc = meged_df_extracted[['reaction', 'result', 'manipulation']].copy()

        # 提取 model_pfba_solution 的 reaction_id
        model_reactions = model_pfba_solution[['reaction_id']].rename(columns={'reaction_id': 'reaction'})

        # 与 meged_df 匹配 reaction，提取对应的 result 和 manipulation
        matched_results = pd.merge(model_reactions, fc, on='reaction', how='left')

        # 生成最终表格，并将无穷值替换为 0
        final_table = matched_results[['reaction', 'result', 'manipulation']]
        final_table['result'] = final_table['result'].replace([np.inf, -np.inf], 0)
        # 保存结果
        #reactiondf.to_csv(os.path.join(path_results3, 'reactiondf.csv'), index=False)


        # 更改列名
        final_table.columns = ['reaction_id', 'FC', 'manipulation']
        final_table.to_csv(os.path.join(path_results3,'fc.csv'))
        model_comparison_results=getinputtable(final_table,model)
        
        bigg_models_metabolites = pd.read_csv(bigg_models_metabolites_file,sep='\t',index_col=0)
        typelist=[ 'Up', 'Down','none']
        # replace NaN in model_comparison_results with 0
        #model_comparison_results.fillna(0, inplace=True)
        model_comparison_results.to_csv(os.path.join(path_results3,'model_comparison_results.csv'))
        #typelist = sorted(typelist,key = lambda i:len(i),reverse=False)
        exclude_rxn = ['EX_co2_e','H2Otex','O2tpp','CO2tpp','O2tex','NH4tex','CO2tex','NH4tpp','H2Otpp','EX_o2_e','EX_h2o_e','EX_nh4_e','EX_o2_e','ADD_h_c','ADD_h_p','ADD_respiratory1']
        total_common_factor_list =  ['hco3_c','ni2_c','cdp_c', 'ag_c', 'dctp_c', 'dutp_c', 'ctp_c', 'gdp_c', 'gtp_c', 'ump_c', 'ca2_c', \
                            'h2o_c', 'datp_c', 'co2_c', 'no2_c', 'no_c', 'k_c', 'zn2_c', 'no3_c', 'o2_c', 'cl_c', 'udp_c', 'damp_c',\
                            'ditp_c', 'dump_c', 'q8h2_c', 'pppi_c', 'idp_c', 'dimp_c', 'pi_c', 'dttp_c', 'so4_c', 'adp_c', 'xtp_c',\
                            'dgtp_c', 'dadp_c', 'coa_c', 'ppi_c', 'h2_c', 'cmp_c', 'fe2_c', 'o2s_c', 'h_c', 'gmp_c', 'itp_c', 'q8_c', \
                            'cobalt2_c', 'n2o_c', 'xmp_c', 'xdp_c', 'nadph_c', 'cu_c', 'cu2_c', 'atp_c', 'dgmp_c', 'imp_c', 'h2s_c', 'utp_c',\
                            'dtmp_c', 'fadh2_c', 'so3_c', 'fad_c', 'cd2_c', 'dgdp_c', 'nad_c', 'nadh_c', 'hg2_c', 'dcmp_c', 'dudp_c', 'dtdp_c',\
                            'didp_c', 'mn2_c', 'dcdp_c', 'nh4_c', 'amp_c', 'fe3_c', 'nadp_c', 'so2_c', 'h2o2_c', 'mg2_c',\
                                'hco3_p','ni2_p','cdp_p', 'ag_p', 'dctp_p', 'dutp_p', 'ctp_p', 'gdp_p', 'gtp_p', 'ump_p', 'ca2_p', \
                            'h2o_p', 'datp_p', 'co2_p', 'no2_p', 'no_p', 'k_p', 'zn2_p', 'no3_p', 'o2_p', 'cl_p', 'udp_p', 'damp_p',\
                            'ditp_p', 'dump_p', 'q8h2_p', 'pppi_p', 'idp_p', 'dimp_p', 'pi_p', 'dttp_p', 'so4_p', 'adp_p', 'xtp_p',\
                            'dgtp_p', 'dadp_p', 'coa_p', 'ppi_p', 'h2_p', 'cmp_p', 'fe2_p', 'o2s_p', 'h_p', 'gmp_p', 'itp_p', 'q8_p', \
                            'cobalt2_p', 'n2o_p', 'xmp_p', 'xdp_p', 'nadph_p', 'cu_p', 'cu2_p', 'atp_p', 'dgmp_p', 'imp_p', 'h2s_p', 'utp_p',\
                            'dtmp_p', 'fadh2_p', 'so3_p', 'fad_p', 'cd2_p', 'dgdp_p', 'nad_p', 'nadh_p', 'hg2_p', 'dcmp_p', 'dudp_p', 'dtdp_p',\
                            'didp_p', 'mn2_p', 'dcdp_p', 'nh4_p', 'amp_p', 'fe3_p', 'nadp_p', 'so2_p', 'h2o2_p', 'mg2_p',\
                                'hco3_e','ni2_e','cdp_e', 'ag_e', 'dctp_e', 'dutp_e', 'ctp_e', 'gdp_e', 'gtp_e', 'ump_e', 'ca2_e', \
                            'h2o_e', 'datp_e', 'co2_e', 'no2_e', 'no_e', 'k_e', 'zn2_e', 'no3_e', 'o2_e', 'cl_e', 'udp_e', 'damp_e',\
                            'ditp_e', 'dump_e', 'q8h2_e', 'pppi_e', 'idp_e', 'dimp_e', 'pi_e', 'dttp_e', 'so4_e', 'adp_e', 'xtp_e',\
                            'dgtp_e', 'dadp_e', 'coa_e', 'ppi_e', 'h2_e', 'cmp_e', 'fe2_e', 'o2s_e', 'h_e', 'gmp_e', 'itp_e', 'q8_e', \
                            'cobalt2_e', 'n2o_e', 'xmp_e', 'xdp_e', 'nadph_e', 'cu_e', 'cu2_e', 'atp_e', 'dgmp_e', 'imp_e', 'h2s_e', 'utp_e',\
                            'dtmp_e', 'fadh2_e', 'so3_e', 'fad_e', 'cd2_e', 'dgdp_e', 'nad_e', 'nadh_e', 'hg2_e', 'dcmp_e', 'dudp_e', 'dtdp_e',\
                            'didp_e', 'mn2_e', 'dcdp_e', 'nh4_e', 'amp_e', 'fe3_e', 'nadp_e', 'so2_e', 'h2o2_e', 'mg2_e','flxso_c','flxr_c','WATER_c','PROTON_c','NADH_c','NAD_c','NADPH_c','NADP_c','ATP_c','ADP_c','Pi_c']  
        metlink_model= build_model_from_FBA_result_cb_plus2(model_comparison_results,typelist,rxn_switch,bigg_models_metabolites,coordinates,substrate_id,product_id)
        id_or_name = True
        flow = True
        All_met, All_rxn= get_model_information(metlink_model)
        common_factor_list=list(set(total_common_factor_list).intersection(set(All_met)))          
        if cofactor_switch == "T":
            d3flux.update_cofactors(metlink_model, common_factor_list)
        html=d3flux.flux_map(metlink_model, overwrite_reversibility=False,
                                excluded_metabolites = common_factor_list,excluded_reactions =exclude_rxn,
                                figsize=(1500, 1500), display_name_format=id_or_name,
                                hide_unused="true", hide_unused_cofactors="true",
                                flowLayout = flow)  # 是否要隐藏flux为0的==
        # save the html file
        # with open((show_factor_file1), "w")as wf:
        #     wf.write(html.data)
        
        color_list=["'#d62728'","'#32CD32'","'#696969'"]
        legend_list=[]
        for i in range(len(typelist)): 
            if typelist[i]=='Down':
                outstr='{ '+"'name':'"+typelist[i]+"', 'color':"+color_list[1]+'}'
                legend_list.append(outstr)
            if typelist[i]=='Up':
                outstr='{ '+"'name':'"+typelist[i]+"', 'color':"+color_list[0]+'}'
                legend_list.append(outstr)
            if typelist[i]=='none':
                outstr='{ '+"'name':'"+typelist[i]+"', 'color':"+color_list[2]+'}'
                legend_list.append(outstr)
        legend_list_out=str(legend_list).replace("\"","")    
        newstr=html.data.replace("'d3flux_legend_color_substitution'",legend_list_out)
        try: 
            with open((show_factor_file), "w")as wf:
                wf.write(newstr) 
        except: # if there is error write the error to the file 
            show_factor_file = os.path.join(path_results3, 'fc_d3flux_error.txt')
            with open((show_factor_file), "w")as wf:
                wf.write(str(sys.exc_info()))
    return newstr

def fc_d3flux_optforce(model_pfba_solution):
    modelname=taskname.split('_')[-1]
    # fc_data_file = './datafile/fc_data.csv'
    # model= cobra.io.load_json_model('/hpcfs/fhome/xuwenqi/project/OptMetarget/INPUT/model/wangry@tib.cas.cn_20241128-131413_iCW773R.json')
    # model_substrate_file='/hpcfs/bdc/Web_Backend/OPTME/new_product_substrate/substrate2/'+modelname+'_substrate_use_list.tsv'
    # model_product_file='/hpcfs/bdc/Web_Backend/OPTME/new_product_substrate/product3/'+modelname+'_product_use_list.tsv'
    bigg_models_metabolites_file = os.path.join(path_task,'bigg_models_metabolites.txt')
    met_C_file = os.path.join(path_task,'met_contain_C_df.tsv')
    substrate_id = inputdic['substrate']
    product_id = inputdic['product']
    # substrate_id2=[substrate_id.replace('SK_','')]
    # product_id2=[product_id.replace('user_add_DM_','')]
    # load model
    # model_substrate = pd.read_csv(model_substrate_file,sep='\t',index_col='substrate_model_correspondind')
    # model_product = pd.read_csv(model_product_file, sep='\t', index_col="obj_model_correspondind")  
    met_C = pd.read_csv(met_C_file,sep='\t', index_col=0)
    biomass=get_initial_obj_info( model)
    rxn_switch="T"
    coordinates =""
    cofactor_switch="F"    
    show_factor_file = os.path.join(path_results4, 'fc_d3flux.html')
    # 提取 reaction, result, manipulation 列
    meged = pd.read_excel(os.path.join(path_results4,'results.xlsx'))
    meged_df_extracted1 = meged[['reaction', 'range_change', 'manipulations']]

    # 根据 manipulation 列调整 result 列的值
    meged_df_extracted1['range_change'] = meged_df_extracted1.apply(
        lambda row: -row['range_change'] if row['manipulations'] == 'down' else row['range_change'], axis=1
    )

    # 生成新表格 fc
    fc1 = meged_df_extracted1[['reaction', 'range_change', 'manipulations']].copy()
    model_reactions1 = model_pfba_solution[['reaction_id']].rename(columns={'reaction_id': 'reaction'})
    # 与 meged_df 匹配 reaction，提取对应的 result 和 manipulation
    matched_result = pd.merge(model_reactions1, fc1, on='reaction', how='left')

    # 生成最终表格，并将无穷值替换为 0
    final_table1 = matched_result[['reaction', 'range_change', 'manipulations']]
    final_table1['range_change'] = final_table1['range_change'].replace([np.inf, -np.inf], 0)

    # 更改列名
    final_table1.columns = ['reaction_id', 'FC', 'manipulation']
    final_table1['manipulation'] = final_table1['manipulation'].replace('Up', 'up')
    final_table1.to_csv(os.path.join(path_results4,'fc.csv'))
    fc_data_file = (os.path.join(path_results4,'fc.csv'))
    fc_data = pd.read_csv(fc_data_file)        
    model_comparison_results=getinputtable(final_table1,model)
    bigg_models_metabolites = pd.read_csv(bigg_models_metabolites_file,sep='\t',index_col=0)
    typelist=[ 'Up', 'Down','none']
    exclude_rxn = ['EX_co2_e','H2Otex','O2tpp','CO2tpp','O2tex','NH4tex','CO2tex','NH4tpp','H2Otpp','EX_o2_e','EX_h2o_e','EX_nh4_e','EX_o2_e','ADD_h_c','ADD_h_p','ADD_respiratory1']
    total_common_factor_list =  ['hco3_c','ni2_c','cdp_c', 'ag_c', 'dctp_c', 'dutp_c', 'ctp_c', 'gdp_c', 'gtp_c', 'ump_c', 'ca2_c', \
                        'h2o_c', 'datp_c', 'co2_c', 'no2_c', 'no_c', 'k_c', 'zn2_c', 'no3_c', 'o2_c', 'cl_c', 'udp_c', 'damp_c',\
                        'ditp_c', 'dump_c', 'q8h2_c', 'pppi_c', 'idp_c', 'dimp_c', 'pi_c', 'dttp_c', 'so4_c', 'adp_c', 'xtp_c',\
                        'dgtp_c', 'dadp_c', 'coa_c', 'ppi_c', 'h2_c', 'cmp_c', 'fe2_c', 'o2s_c', 'h_c', 'gmp_c', 'itp_c', 'q8_c', \
                        'cobalt2_c', 'n2o_c', 'xmp_c', 'xdp_c', 'nadph_c', 'cu_c', 'cu2_c', 'atp_c', 'dgmp_c', 'imp_c', 'h2s_c', 'utp_c',\
                        'dtmp_c', 'fadh2_c', 'so3_c', 'fad_c', 'cd2_c', 'dgdp_c', 'nad_c', 'nadh_c', 'hg2_c', 'dcmp_c', 'dudp_c', 'dtdp_c',\
                        'didp_c', 'mn2_c', 'dcdp_c', 'nh4_c', 'amp_c', 'fe3_c', 'nadp_c', 'so2_c', 'h2o2_c', 'mg2_c',\
                            'hco3_p','ni2_p','cdp_p', 'ag_p', 'dctp_p', 'dutp_p', 'ctp_p', 'gdp_p', 'gtp_p', 'ump_p', 'ca2_p', \
                        'h2o_p', 'datp_p', 'co2_p', 'no2_p', 'no_p', 'k_p', 'zn2_p', 'no3_p', 'o2_p', 'cl_p', 'udp_p', 'damp_p',\
                        'ditp_p', 'dump_p', 'q8h2_p', 'pppi_p', 'idp_p', 'dimp_p', 'pi_p', 'dttp_p', 'so4_p', 'adp_p', 'xtp_p',\
                        'dgtp_p', 'dadp_p', 'coa_p', 'ppi_p', 'h2_p', 'cmp_p', 'fe2_p', 'o2s_p', 'h_p', 'gmp_p', 'itp_p', 'q8_p', \
                        'cobalt2_p', 'n2o_p', 'xmp_p', 'xdp_p', 'nadph_p', 'cu_p', 'cu2_p', 'atp_p', 'dgmp_p', 'imp_p', 'h2s_p', 'utp_p',\
                        'dtmp_p', 'fadh2_p', 'so3_p', 'fad_p', 'cd2_p', 'dgdp_p', 'nad_p', 'nadh_p', 'hg2_p', 'dcmp_p', 'dudp_p', 'dtdp_p',\
                        'didp_p', 'mn2_p', 'dcdp_p', 'nh4_p', 'amp_p', 'fe3_p', 'nadp_p', 'so2_p', 'h2o2_p', 'mg2_p',\
                            'hco3_e','ni2_e','cdp_e', 'ag_e', 'dctp_e', 'dutp_e', 'ctp_e', 'gdp_e', 'gtp_e', 'ump_e', 'ca2_e', \
                        'h2o_e', 'datp_e', 'co2_e', 'no2_e', 'no_e', 'k_e', 'zn2_e', 'no3_e', 'o2_e', 'cl_e', 'udp_e', 'damp_e',\
                        'ditp_e', 'dump_e', 'q8h2_e', 'pppi_e', 'idp_e', 'dimp_e', 'pi_e', 'dttp_e', 'so4_e', 'adp_e', 'xtp_e',\
                        'dgtp_e', 'dadp_e', 'coa_e', 'ppi_e', 'h2_e', 'cmp_e', 'fe2_e', 'o2s_e', 'h_e', 'gmp_e', 'itp_e', 'q8_e', \
                        'cobalt2_e', 'n2o_e', 'xmp_e', 'xdp_e', 'nadph_e', 'cu_e', 'cu2_e', 'atp_e', 'dgmp_e', 'imp_e', 'h2s_e', 'utp_e',\
                        'dtmp_e', 'fadh2_e', 'so3_e', 'fad_e', 'cd2_e', 'dgdp_e', 'nad_e', 'nadh_e', 'hg2_e', 'dcmp_e', 'dudp_e', 'dtdp_e',\
                        'didp_e', 'mn2_e', 'dcdp_e', 'nh4_e', 'amp_e', 'fe3_e', 'nadp_e', 'so2_e', 'h2o2_e', 'mg2_e','flxso_c','flxr_c','WATER_c','PROTON_c','NADH_c','NAD_c','NADPH_c','NADP_c','ATP_c','ADP_c','Pi_c']  
    metlink_model= build_model_from_FBA_result_cb_plus2(model_comparison_results,typelist,rxn_switch,bigg_models_metabolites,coordinates,substrate_id,product_id)
    id_or_name = True
    flow = True
    All_met, All_rxn= get_model_information(metlink_model)
    common_factor_list=list(set(total_common_factor_list).intersection(set(All_met)))          
    if cofactor_switch == "T":
        d3flux.update_cofactors(metlink_model, common_factor_list)
    html = d3flux.flux_map(metlink_model, overwrite_reversibility=False,
                            excluded_metabolites = common_factor_list,excluded_reactions =exclude_rxn,
                            figsize=(1500, 1500), display_name_format=id_or_name,
                            hide_unused="true", hide_unused_cofactors="true",
                            flowLayout = flow)  # 是否要隐藏flux为0的
    color_list=["'#d62728'","'#32CD32'","'#696969'"]
    legend_list=[]
    for i in range(len(typelist)): 
        if typelist[i]=='Down':
            outstr='{ '+"'name':'"+typelist[i]+"', 'color':"+color_list[1]+'}'
            legend_list.append(outstr)
        if typelist[i]=='Up':
            outstr='{ '+"'name':'"+typelist[i]+"', 'color':"+color_list[0]+'}'
            legend_list.append(outstr)
        if typelist[i]=='none':
            outstr='{ '+"'name':'"+typelist[i]+"', 'color':"+color_list[2]+'}'
            legend_list.append(outstr)
    legend_list_out=str(legend_list).replace("\"","")    
    newstr=html.data.replace("'d3flux_legend_color_substitution'",legend_list_out)
    with open((show_factor_file), "w")as wf:
        wf.write(newstr)                
    return newstr    

if __name__=="__main__":
    path_model=sys.argv[1]
    path_task=sys.argv[2]
    path_map=sys.argv[3]
    path_results=sys.argv[4]
    taskname=sys.argv[5]
    path_results2=os.path.join(path_results,taskname)
    path_results3 = os.path.join(path_results2, "FSEOF")
    path_results4 = os.path.join(path_results2, "OptForce")
    if not os.path.exists(path_results2):
        os.makedirs(path_results2)
    if not os.path.exists(path_results3):
        os.makedirs(path_results3)
    if not os.path.exists(path_results4):
        os.makedirs(path_results4)
    cobra_config = cobra.Configuration()
    cobra_config.solver = os.environ.get("OPTME_COBRA_SOLVER", "cplex")
    model,inputdic = prepare_model(path_model,path_task,taskname)
    model_pfba_solution = reaction_list_d3(model,inputdic)
    if inputdic['taskname'] == 'FSEOF':
        product,over_flux_fseof =calculate_product(inputdic)
        v_biomass,wild_flux_fseof = calculate_biomass(inputdic)
        reactiondf = biomass(product,inputdic)
        reactiondf['manipulation'] = reactiondf.apply(lambda row: check_monotonicity(row), axis=1)
        gpr_dict = get_gpr_for_reactions(reactiondf, model)
        reactiondf['gpr'] = reactiondf['reaction'].map(gpr_dict)
        reactiondf['manipulation'] = reactiondf.apply(lambda row: None if check_opposite_signs(row) else row['manipulation'], axis=1)
        up_reactions,down_reactions,result_per_row,reactiondf = threshold(reactiondf)
        reaction_list=[]
        for rea in model.reactions:
            reaction_list.append(rea.id)
        equation_dict =gene_reaction_map(reaction_list,model)
        reactiondf['equation'] = reactiondf['reaction'].map(equation_dict)
        Yield0,Carbon_Yield0,Ycm,Mass_Yield0,product_c_num,substrate_c_num,substrate,substrate_molecular,product_molecular=calculated_yield(path_task,inputdic,model,product)
        degree_pro,degree_C_pro =reduced_degree(model,inputdic['product'])
        degree_bio,degree_C_bio =reduced_degree(model,inputdic['substrate'])
        reactions,up_reaction,down_reaction =detail_put(model,reactiondf,path_results3)
        result_dict = up_picture(reactiondf)
        result_dict_down = down_picture(reactiondf)

        os.makedirs(os.path.join(path_results3,'target_fig'),exist_ok=True)
        if len(result_dict)>0:
            drawtarget(result_dict,'FSEOF', os.path.join(path_results3,'target_fig'))
        if len(result_dict_down)>0:
            drawtarget(result_dict_down,'FSEOF', os.path.join(path_results3,'target_fig'))

        output =output_web(reactions)
        reactions_with_ranges = conbine_flux(model)
        source_file = os.path.join(path_map, 'iJO1366.Central_metabolism.json')
        destination_file = os.path.join(path_results3, 'map.json')
        shutil.copy2(source_file, destination_file)
        os.rename(destination_file, os.path.join(path_results3, 'map.json'))
        try:
            model_pfba_solutionA_fseof,model_pfba_solutionB_fseof,model_pfba_solutionC_fseof,visualisation_html_file,visualisation_html_file_WT    =  get_d3_flux(model_pfba_solution,wild_flux_fseof,over_flux_fseof,reactions_with_ranges,path_results3)
            reactiondf.to_csv(os.path.join(path_results3, 'reactiondf.csv'), index=False)
            newstr = fc_d3fulx(reactiondf)
        except Exception as e:
            print(f"An error occurred: {e}")

    if inputdic['taskname'] == 'loopless_optforce_MUST':
        model,inputdic =prepare_model_optforce(path_model,path_task,taskname)
        v0_biomass,fva_wild,wild_flux=calculate_biomass_optforce(inputdic,model,path_results4)
        v1_product_max,fva_over,over_flux = calculate_product_optforce(inputdic,v0_biomass,model,path_results4)
        enzyme_results_data,enzyme_overresults2_data = read_file(path_results4)
        ko_data,up_data,down_data,range_change,mean_change,direction=compare_results(enzyme_results_data, enzyme_overresults2_data)
        reaction_list=[]
        for rea in model.reactions:
            reaction_list.append(rea.id)
        equation_dict = gene_reaction_map(reaction_list,model)
        gpr_dict = get_gpr_for_reactions_optforce(reaction_list, model)
        meged_df = must_df(enzyme_results_data,enzyme_overresults2_data,range_change,ko_data,up_data,down_data)
        Yield0,Carbon_Yield0,Ycm,Mass_Yield0,product_c_num,substrate_c_num,substrate,substrate_molecular,product_molecular = calculated_yield_optforce(path_task,inputdic,model,v1_product_max)
        degree_pro,degree_C_pro =reduced_degree(model,inputdic['product'])
        degree_bio,degree_C_bio =reduced_degree(model,inputdic['substrate'])
        reactions,up_reaction,down_reaction,ko_reaction =detail_put_optforce(model,meged_df,path_results4)
        output = output_web_optforce(reactions)
        reactions_with_ranges = conbine_flux_optforce(model)
        result_json1 = ko_picture_optforce(meged_df)
        result_json2 = up_picture_optforce(meged_df)
        result_json3 = down_picture_optforce(meged_df)

        os.makedirs(os.path.join(path_results4,'target_fig'),exist_ok=True)
        with open(os.path.join(path_results4, 'ko.json'),encoding='utf-8') as fp:
            result_json1=json.load(fp)
        # 检查 result_json1 是否为空
        if isinstance(result_json1, dict) and len(result_json1) > 0:
            drawtarget(result_json1, 'OptForce', os.path.join(path_results4, 'target_fig'))
        else:
            print("Warning: 'ko' is empty or not a dictionary. Skipping drawtarget.")
        with open(os.path.join(path_results4, 'up.json'),encoding='utf-8') as fp:
            result_json2=json.load(fp)
        if isinstance(result_json2, dict) and len(result_json2) > 0:
            drawtarget(result_json2, 'OptForce', os.path.join(path_results4, 'target_fig'))
        else:
            print("Warning: 'up' is empty or not a dictionary. Skipping drawtarget.")
        with open(os.path.join(path_results4, 'down.json'),encoding='utf-8') as fp:
            result_json3=json.load(fp)        
        if isinstance(result_json3, dict) and len(result_json3) > 0:
            drawtarget(result_json3, 'OptForce', os.path.join(path_results4, 'target_fig'))
        else:
            print("Warning: 'down' is empty or not a dictionary. Skipping drawtarget.")
        meged = pd.read_excel(os.path.join(path_results4,'results.xlsx'))
        model_pfba_solution.to_csv(os.path.join(path_results4,'pfbatest.csv'))
        source_file = os.path.join(path_map, 'iJO1366.Central_metabolism.json')
        destination_file = os.path.join(path_results4, 'map.json')
        shutil.copy2(source_file, destination_file)
        os.rename(destination_file, os.path.join(path_results4, 'map.json'))
        try:
            model_pfba_solutionA_optforce,model_pfba_solutionB_optforce,model_pfba_solutionC_optforce,visualisation_html_file,visualisation_html_file_WT    =  get_d3_flux(model_pfba_solution,wild_flux,over_flux,reactions_with_ranges,path_results4)
            newstr = fc_d3flux_optforce(model_pfba_solution)
        except Exception as e:
            print(f"An error occurred: {e}")
