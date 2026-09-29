# coding:utf-8
"""
Author: Jingyi Cai, Wenqi Xu (2024-2026)
Function: Run each method named in the task and combine their target lists. FSEOF and loopless OptForce call target.py. E_FSEOF and E_OptForce call EToptme.py. llm calls llmoptme.py.
Input: python summary.py <model_dir> <task_dir> <map_dir> <results_dir> <task_id>. Reads <task_dir>/<task_id>.json. The selected solver comes from OPTME_COBRA_SOLVER and OPTME_PYOMO_SOLVER.
Output: One folder per method under <results_dir>/<task_id>/, plus a combined <results_dir>/<task_id>/output.json. Prints METHOD_RUNTIME lines for the log.
"""
import sys
import cobra
import cProfile
# sys.path.append('/home/sun/ETGEMS-10.20')
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
import shutil
import os
import json
import copy
import subprocess
import time
from visualization2 import *
from task_config import OPTME_PYTHON, IBRIDGE_PYTHON, TARGET_SCRIPT, IBRIDGE_SCRIPT, ETOPTME_SCRIPT, LLM_PYTHON, LLMOPTME_SCRIPT
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
    # model_file0 = os.path.join(path_task, 'iCW773_uniprot_modification_del.json')
    # model0=cobra.io.load_json_model(model_file0)
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
    return model,inputdic      
def biomass(inputdic,model):
    model.objective=inputdic['biomass']
    v0_biomass=model.optimize().objective_value    
    return v0_biomass
    # v0_biomass =  biomass(inputdic,model)
def product(inputdic,model):
    model.reactions.get_by_id(inputdic['biomass']).bounds=(v0_biomass*0.1,v0_biomass*0.1)
    model.objective=inputdic['product']
    v1_product_max=model.optimize().objective_value    
    return v1_product_max


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
    product_c_num = 0
    substrate_c_num = 0
    substrate_molecular = 0
    product_molecular = 0
    substrate = 0
    if inputdic['substrate'].endswith('_reverse'):
        inputdic['substrate'] = inputdic['substrate'][:-len('_reverse')]
    solution=model.optimize()
    substrate=solution.fluxes.get(inputdic['substrate'], 0)
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
        print(f"product_c_number: {product_c_num}")
        product_molecular=product_molecular_r-product_molecular_l
        if product_c_num <= 0:
            print(f"Skip carbon-yield calculation for {inputdic['product']}: product_c_num={product_c_num}")
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
            if inputdic['substrate'] in list(list(solution_select.index)) and substrate:
                Yield0=abs(round(prod0/solution_select[inputdic['substrate']],3))
                if substrate_molecular:
                    Mass_Yield0=abs(round((product_molecular*prod0)/(substrate_molecular*substrate),3)) 
            # elif 'EX_glc_e_reverse'+"_reverse" in list(flux_solution_df_select.index):
            #     result_static['Yield']=abs(round(flux_solution_df_select.loc[product_id,'fluxes']/flux_solution_df_select.loc['EX_glc_e_reverse'+"_reverse",'fluxes'],7))
    return Yield0,Carbon_Yield0,Ycm,Mass_Yield0,product_c_num,substrate_c_num,substrate,substrate_molecular,product_molecular
#  Yield0,Carbon_Yield0,Mass_Yield0,product_c_num,substrate_c_



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
    model_pfba_solution=model_pfba_solution.sort_values(by=['abs_fluxes'],ascending=False)    
    return model_pfba_solution

def get_d3_flux(model_pfba_solution,fc_flux,path_results2):
    model_pfba_solutionA = model_pfba_solution.copy()
    model_pfba_solutionB = model_pfba_solution.copy()
    model_pfba_solutionC = model_pfba_solution.copy()
    # # 填充 wild_flux 值到 model_pfba_solutionA 的 fluxes 列，未找到则填充 0
    # model_pfba_solutionA['fluxes'] = model_pfba_solutionA.index.map(wild_flux).fillna(0)
    # # 填充 over_flux 值到 model_pfba_solutionB 的 fluxes 列，未找到则填充 0
    # model_pfba_solutionB['fluxes'] = model_pfba_solutionB.index.map(over_flux).fillna(0)
    # 填充 fc_flux 值到 model_pfba_solutionC 的 fluxes 列，未找到则填充 0
    model_pfba_solutionC['fluxes'] = model_pfba_solutionC.index.map(fc_flux).fillna(0)
    # 替换 abs_fluxes 列为 fluxes 的绝对值
    # model_pfba_solutionA['abs_fluxes'] = model_pfba_solutionA['fluxes'].abs()
    # model_pfba_solutionB['abs_fluxes'] = model_pfba_solutionB['fluxes'].abs()
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
    

    return model_pfba_solutionC,visualisation_html_file,visualisation_html_file_WT,visualisation_html_file_fc      
        



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

def output(v0_biomass,v1_product_max):
    output={}
    output['summary'] = {}
    output['yield'] = {}
    summary = {
        'Prod_Rate(mmol/gDW/h)':round(v1_product_max,3),
        'Growth(1/h)':round(v0_biomass,3),
        'Yield': {
        'Ycm':Ycm,
        'Ycg':Ycm,
        'Yrm': round(degree_bio / degree_pro, 3) if degree_pro else 'N/A',
        'Yrg': round((round(degree_bio / degree_pro, 3) * product_molecular / substrate_molecular), 3) if degree_pro and substrate_molecular else 'N/A',
        'Ysm':Yield0,
        'Ysg':Mass_Yield0,
        'Yscm':Carbon_Yield0,
        'Yscg':Carbon_Yield0
    }    
    }
    output['summary'] = summary
    output_file = os.path.join(path_results2, 'summary.json')
    with open(output_file, 'w') as json_file:
        json.dump(output, json_file)     
    return output

def remove_path(output, method,type):
    if type=="d3flux":
        str1=method+"_d3flux"+"_"+"wild"
        str2=method+"_d3flux"+"_"+"over"
        str3=method+"_d3flux"+"_"+"fc"
        output["data"]["Visualization"][str1]=''
        output["data"]["Visualization"][str2]=''
        output["data"]["Visualization"][str3]=''
    if type=="escher":
        str1=method+"_escher"+"_"+"wild"
        str2=method+"_escher"+"_"+"over"
        str3=method+"_escher"+"_"+"fc"
        output["data"]["Visualization"][str1]=''
        output["data"]["Visualization"][str2]=''
        output["data"]["Visualization"][str3]=''
    return output

def _format_duration(seconds):
    total = int(max(0.0, float(seconds)))
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f"{hours}:{minutes:02d}:{secs:02d}"
    return f"{minutes:02d}:{secs:02d}"


def _print_method_runtime(method_times):
    """Print one line per method so the pipeline log can collect them."""
    print("=== Method runtime ===")
    total = 0.0
    for name, elapsed, status in method_times:
        total += elapsed
        print(f"  {name:<36} {_format_duration(elapsed):>10}  {elapsed:10.1f} s  {status}")
    print(f"  {'methods total':<36} {_format_duration(total):>10}  {total:10.1f} s")


def run_metabolic_tasks(path_results, path_task, path_model, path_map, taskname):
    """
    Executes specified metabolic tasks, aggregates results, and generates a combined output.
    
    Args:
        path_results (str): Results file path.
        path_task (str): Task file path.
        path_model (str): Model file path.
        path_map (str): Map file path.
        taskname (str): Task name.

    Returns:
        None
    """
    #print('hello world')
    # Load the input JSON file
    with open(os.path.join(path_task, taskname) + '.json', encoding='utf-8') as fp:
        inputdic = json.load(fp)
    all_method = inputdic['taskname']
    # Create a deepcopy of the input dictionary for modifications
    methods = inputdic['taskname']
    if isinstance(methods, str):
        methods = [methods]
    inputdic_tmp = copy.deepcopy(inputdic)

    # Initialize the final output dictionary
    if inputdic.get('yield_na', False):
        yield_summary = {
            'Ycm': 'N/A',
            'Ycg': 'N/A',
            'Yrm': 'N/A',
            'Yrg': 'N/A',
            'Ysm': 'N/A',
            'Ysg': 'N/A',
            'Yscm': 'N/A',
            'Yscg': 'N/A'
        }
    else:
        yrm = round(degree_bio / degree_pro, 3) if degree_pro else 'N/A'
        yrg = round(yrm * product_molecular / substrate_molecular, 3) if yrm != 'N/A' and substrate_molecular else 'N/A'
        yield_summary = {
            'Ycm': Ycm,
            'Ycg': Ycm,
            'Yrm': yrm,
            'Yrg': yrg,
            'Ysm': Yield0,
            'Ysg': Mass_Yield0,
            'Yscm': Carbon_Yield0,
            'Yscg': Carbon_Yield0
        }

    summary = {
        'Prod_Rate(mmol/gDW/h)': round(v1_product_max, 3),
        'Growth(1/h)': round(v0_biomass, 3),
        'Yield': yield_summary
    }
    final_output = {
        "summary":summary,
        "data": {}
    }

    # Initialize combined summary counters
    combined_summary = {"UP": 0, "Down": 0, "KO": 0}
    combined_Enzyme_summary = {"UP": 0, "Down": 0, "KO": 0}
    combined_up = []
    combined_down = []
    combined_Enzyme_up =[]
    combined_Enzyme_down = []

    def _to_number(value, default=0.0):
        try:
            if value is None:
                return default
            return float(value)
        except (TypeError, ValueError):
            return default

    def _extract_enzyme_cost(value):
        if isinstance(value, list):
            if not value:
                return 0.0
            return _to_number(value[0], 0.0)
        return _to_number(value, 0.0)

    def _is_e_series_method(method_name):
        """ET-OptME enzyme methods only (E_FSEOF, E_OptForce, ...). Excludes llm."""
        return str(method_name).startswith('E_')

    def _summary_mod_counts(summary):
        return {
            "UP": int(summary.get("UP", 0)),
            "Down": int(summary.get("Down", 0)),
            "KO": int(summary.get("KO", 0)),
        }

    # Iterate through taskname methods
    print(inputdic)
    method_times = []
    for method in inputdic['taskname']:

        # Set taskname to the current method
        inputdic_tmp['taskname'] = method
        tmppath = os.path.join(path_task, taskname) + '.json'

        # Save the modified JSON file
        with open(tmppath, 'w', encoding='utf-8') as file:
            json.dump(inputdic_tmp, file)
        # Construct the command based on the method
        print(method)
        if method in ["FSEOF", "loopless_optforce_MUST"]:
            command = [
                OPTME_PYTHON, TARGET_SCRIPT,
                path_model,
                path_task,
                path_map,
                path_results,
                taskname
            ]
        elif method == "iBridge":
            command = [
                IBRIDGE_PYTHON, IBRIDGE_SCRIPT,
                path_model,
                path_task,
                path_map,
                path_results,
                taskname
            ]
        elif method in ["E_FSEOF", "E_OptForce"]:
            command = [
                OPTME_PYTHON, ETOPTME_SCRIPT,
                path_model,
                path_task,
                path_map,
                path_results,
                taskname
            ]
        elif method in ["llm"]:
            command = [
                LLM_PYTHON, LLMOPTME_SCRIPT,
                path_model,
                path_task,
                path_map,
                path_results,
                taskname
            ]
        print(command)
        # Execute the command and capture output
        method_name = method
        started = time.perf_counter()
        result = subprocess.run(command, capture_output=True, text=True)
        elapsed = time.perf_counter() - started
        status = "ok" if result.returncode == 0 else "failed"
        method_times.append((method_name, elapsed, status))
        print(f"METHOD_RUNTIME\t{method_name}\t{elapsed:.3f}\t{status}")
        print(f"Method {method_name} runtime: {elapsed:.1f} s ({status})")

        # Check if the command executed successfully
        if result.returncode != 0:
            print(f"Method {method} failed: {result.stderr}")
            continue

                # Replace loopless_optforce_MUST with OptForce
        if method == "loopless_optforce_MUST":
            method = "OptForce"   
        # Read the output JSON file for the current method
        output_file_path = os.path.join(path_results, taskname, method, "output.json")
        try:
            with open(output_file_path, encoding='utf-8') as fp:
                method_data = json.load(fp)
            if 'E_' not in method and method != 'llm':
                reactions = method_data.get("reaction", {})
                for _, details in reactions.items():
                    if not isinstance(details, dict):
                        continue
                    fold_change = _to_number(details.get("Fold_Change"), 0.0)
                    average_flux = abs(_to_number(details.get("Average_Flux"), 0.0))
                    details["score"] = round(abs(fold_change) * average_flux, 6)
                _raw_scores = [d["score"] for d in reactions.values() if isinstance(d, dict) and "score" in d]
                _max_score = max(_raw_scores) if _raw_scores else 0.0
                if _max_score > 0:
                    for _, details in reactions.items():
                        if isinstance(details, dict) and "score" in details:
                            details["score"] = round(details["score"] / _max_score, 6)

                # Add the method data to the final output
                final_output["data"][method] = {
                    "summary": method_data.get("summary", {}),
                    "reaction": reactions
                }

                # Update combined summary counters
                method_summary = method_data.get("summary", {})
                combined_summary["UP"] += method_summary.get("UP", 0)
                combined_summary["Down"] += method_summary.get("Down", 0)
                combined_summary["KO"] += method_summary.get("KO", 0)

                # Append combined up and down reaction data
                for reaction, details in reactions.items():
                    if details["Modification"] == "UP":
                        combined_up.append(details)
                    elif details["Modification"] == "Down":
                        combined_down.append(details)
            else:
                genes = method_data.get("genes", {})
                if method in ["E_FSEOF", "E_OptForce"]:
                    for _, details in genes.items():
                        if not isinstance(details, dict):
                            continue
                        fold_change = _to_number(details.get("Fold_Change"), 0.0)
                        enzyme_cost = abs(_extract_enzyme_cost(details.get("Enzyme_Cost", 0.0)))
                        details["score"] = round(abs(fold_change) * enzyme_cost, 6)
                    _raw_scores = [d["score"] for d in genes.values() if isinstance(d, dict) and "score" in d]
                    _max_score = max(_raw_scores) if _raw_scores else 0.0
                    if _max_score > 0:
                        for _, details in genes.items():
                            if isinstance(details, dict) and "score" in details:
                                details["score"] = round(details["score"] / _max_score, 6)
                elif method == "llm":
                    for _, details in genes.items():
                        if isinstance(details, dict):
                            details["Fold_change"] = ""
                            details["score"] = ""

                final_output["data"][method] = {
                    "summary": method_data.get("summary", {}),
                    "genes": genes 
                }
                # combined_Enzyme = sum of E_* method summaries only (not llm)
                if _is_e_series_method(method):
                    method_summary = _summary_mod_counts(method_data.get("summary", {}))
                    combined_Enzyme_summary["UP"] += method_summary["UP"]
                    combined_Enzyme_summary["Down"] += method_summary["Down"]
                    combined_Enzyme_summary["KO"] += method_summary["KO"]

                    for reaction, details in genes.items():
                        if not isinstance(details, dict):
                            continue
                        mod = str(details.get("Modification", "")).upper()
                        if mod == "UP":
                            combined_Enzyme_up.append(details)
                        elif mod == "DOWN":
                            combined_Enzyme_down.append(details)

        except FileNotFoundError:
            print(f"Warning: Output file not found: {output_file_path}")
        except json.JSONDecodeError:
            print(f"Error: Unable to parse output file: {output_file_path}")

    # Add combined summary and file paths to the final output
    final_output["data"]["combined"] = {
        "summary": combined_summary,
        "combined_up": os.path.join(path_results2, "UP_modification.tsv"),
        "combined_down": os.path.join(path_results2, "down_modification.tsv")
    }
    has_e_series_method = any(_is_e_series_method(m) for m in all_method) if isinstance(all_method, list) else _is_e_series_method(all_method)
    if has_e_series_method:
        final_output["data"]["combined_Enzyme"] = {
            "summary": combined_Enzyme_summary,
            "combined_Enzyme_up": os.path.join(path_results2, "UP_modification_E.tsv"),
            "combined_Enzyme_down": os.path.join(path_results2, "down_modification_E.tsv")
        }

    # --- combined_AI from llm results ---
    _all_methods_list = all_method if isinstance(all_method, list) else [all_method]
    if "llm" in _all_methods_list:
        llm_block = final_output["data"].get("llm", {})
        llm_genes = llm_block.get("genes", {})
        combined_AI_up = []
        combined_AI_down = []
        combined_AI_ko = []
        for _, ginfo in llm_genes.items():
            if not isinstance(ginfo, dict):
                continue
            mod = str(ginfo.get("Modification", "")).upper()
            if mod == "UP":
                combined_AI_up.append(ginfo)
            elif mod == "DOWN":
                combined_AI_down.append(ginfo)
            elif mod in ("KO", "DEL"):
                combined_AI_ko.append(ginfo)
        ai_columns = ["Gene_Name", "Enzyme_Name", "Modification", "Fold_change",
                       "Reaction_Name", "Reaction Formula", "Type", "Rationale"]
        ai_up_path = os.path.join(path_results2, "UP_modification_AI.tsv")
        ai_down_path = os.path.join(path_results2, "down_modification_AI.tsv")
        ai_ko_path = os.path.join(path_results2, "KO_modification_AI.tsv")
        for fpath, rows in [
            (ai_up_path, combined_AI_up),
            (ai_down_path, combined_AI_down),
            (ai_ko_path, combined_AI_ko),
        ]:
            with open(fpath, "w", encoding="utf-8") as fout:
                fout.write("\t".join(ai_columns) + "\n")
                for item in rows:
                    fout.write("\t".join(str(item.get(c, "")) for c in ai_columns) + "\n")
        # combined_AI mirrors llm summary (single AI method)
        combined_AI_summary = _summary_mod_counts(llm_block.get("summary", {}))
        final_output["data"]["combined_AI"] = {
            "summary": combined_AI_summary,
            "combined_AI_up": ai_up_path,
            "combined_AI_down": ai_down_path,
            "combined_AI_ko": ai_ko_path,
        }

    # --- Build Visualization dict dynamically based on task methods ---
    _viz_methods = {"FSEOF", "loopless_optforce_MUST", "OptForce", "iBridge"}
    _folder_map = {"loopless_optforce_MUST": "OptForce"}
    visualization = {}
    for m in _all_methods_list:
        if m not in _viz_methods:
            continue
        folder = _folder_map.get(m, m)
        visualization[f"{folder}_escher_wild"] = os.path.join(path_results2, folder, "wildflux_map.json")
        visualization[f"{folder}_escher_over"] = os.path.join(path_results2, folder, "overflux_map.json")
        visualization[f"{folder}_escher_fc"]   = os.path.join(path_results2, folder, "fcflux_map.json")
        visualization[f"{folder}_d3flux_wild"] = os.path.join(path_results2, folder, "wild_d3flux.html")
        visualization[f"{folder}_d3flux_over"] = os.path.join(path_results2, folder, "over_d3flux.html")
        visualization[f"{folder}_d3flux_fc"]   = os.path.join(path_results2, folder, "fc_d3flux.html")
    final_output["data"]["Visualization"] = visualization
    # Save the final output to a JSON file
    combined_output_path = os.path.join(path_results, taskname, "output.json")
    inputdic['taskname']= all_method
    with open(combined_output_path, 'w', encoding='utf-8') as fp:
        json.dump(final_output, fp, indent=4, ensure_ascii=False)
    # check if the taskname is already in the path
    tmppath=os.path.join(path_task,taskname)+'.json'
    with open(tmppath, 'w') as file:
        json.dump(inputdic, file)
    print(f"Combined results saved to {combined_output_path}.")
    _print_method_runtime(method_times)


if __name__=="__main__":
    path_model=sys.argv[1]
    path_task=sys.argv[2]
    path_map=sys.argv[3]
    path_results=sys.argv[4]
    taskname=sys.argv[5]
    path_results2=os.path.join(path_results,taskname)
    path_results3=os.path.join(path_results2, "iBridge")
    if not os.path.exists(path_results2):
        os.makedirs(path_results2)
    cobra.Configuration().solver = os.environ.get("OPTME_COBRA_SOLVER", "cplex")
    model,inputdic = prepare_model(path_model,path_task,taskname)
    v0_biomass =  biomass(inputdic,model)
    v1_product_max = product(inputdic,model)

    if inputdic.get('yield_na', False):
        print(f"Skip yield calculation: {inputdic.get('yield_na_reason', '')}")
        Yield0 = 'N/A'
        Carbon_Yield0 = 'N/A'
        Ycm = 'N/A'
        Mass_Yield0 = 'N/A'
        product_c_num = 'N/A'
        substrate_c_num = 'N/A'
        substrate = 'N/A'
        substrate_molecular = 'N/A'
        product_molecular = 'N/A'
        degree_pro = 'N/A'
        degree_C_pro = 'N/A'
        degree_bio = 'N/A'
        degree_C_bio = 'N/A'
    else:
        Yield0,Carbon_Yield0,Ycm,Mass_Yield0,product_c_num,substrate_c_num,substrate,substrate_molecular,product_molecular = calculated_yield_optforce(path_task,inputdic,model,v1_product_max)
        degree_pro,degree_C_pro =reduced_degree(model,inputdic['product'])
        degree_bio,degree_C_bio =reduced_degree(model,inputdic['substrate'])
    run_metabolic_tasks(path_results, path_task, path_model, path_map, taskname)

    try:
        if "iBridge" in inputdic['taskname'] :
            model_pfba_solution = reaction_list_d3(model,inputdic)
            with open(os.path.join(path_results3, 'fcflux_map.json'),encoding='utf-8') as fp:
                reactions_with_ranges=json.load(fp)    
            model_pfba_solutionC,visualisation_html_file,visualisation_html_file_WT,visualisation_html_file_fc =  get_d3_flux(model_pfba_solution,reactions_with_ranges,path_results3)
    except Exception as e:
        print(f"An error occurred: {e}")

    # Generate fc_d3flux.html for FSEOF and OptForce
    _fc_gen = {"FSEOF": "FSEOF", "loopless_optforce_MUST": "OptForce"}
    for _task_m, _folder in _fc_gen.items():
        if _task_m not in inputdic['taskname']:
            continue
        _method_dir = os.path.join(path_results2, _folder)
        _fc_json = os.path.join(_method_dir, "fcflux_map.json")
        _fc_html = os.path.join(_method_dir, "fc_d3flux.html")
        if not os.path.isfile(_fc_json) or os.path.isfile(_fc_html):
            continue
        try:
            with open(_fc_json, encoding="utf-8") as _fp:
                _fc_flux = json.load(_fp)
            _fc_tsv = os.path.join(_method_dir, "fc_d3flux_viz.tsv")
            with open(_fc_tsv, "w") as _fout:
                for _rxn_id, _flux_val in _fc_flux.items():
                    try:
                        _rxn = model.reactions.get_by_id(_rxn_id)
                        _check = _rxn.check_mass_balance()
                        _fout.write(f"{_rxn_id}\t{round(_flux_val,3)}\t{_rxn.reaction}\t{_rxn.build_reaction_string(use_metabolite_names=True)}\t{_rxn.bounds}\t{_check}\n")
                    except Exception:
                        pass
            pathway_tsv_visualization(
                _fc_tsv, inputdic["substrate"], inputdic["product"],
                cofactor_switch="T", flow="True",
                outputdir=_method_dir, rxnId="fc_d3flux", model=model
            )
            print(f"Generated fc_d3flux.html for {_folder}")
        except Exception as _e:
            print(f"Warning: fc_d3flux generation for {_folder} failed: {_e}")
    methods = inputdic['taskname']

    # --- Post-run validation ---
    taskidmap = {
        'FSEOF': 'FSEOF',
        'loopless_optforce_MUST': 'OptForce',
        'iBridge': 'iBridge',
        'E_FSEOF': 'E_FSEOF',
        'E_OptForce': 'E_OptForce',
        'llm': 'llm',
    }
    d3flux_methods = {'FSEOF', 'loopless_optforce_MUST', 'iBridge'}

    checkdic = pd.DataFrame({
        'target_prediction': [0] * len(methods),
        'd3flux': [0] * len(methods),
        'errorstring': [''] * len(methods),
    }, index=methods)

    for m in methods:
        folder_name = taskidmap.get(m, m)
        out_path = os.path.join(path_results2, folder_name, 'output.json')
        if os.path.isfile(out_path) and os.path.getsize(out_path) > 0:
            checkdic.loc[m, 'target_prediction'] = 1
        else:
            checkdic.loc[m, 'errorstring'] += f' {folder_name} target prediction failed;'

    for m in methods:
        if m in d3flux_methods:
            folder_name = taskidmap.get(m, m)
            html_path = os.path.join(path_results2, folder_name, 'fc_d3flux.html')
            if os.path.isfile(html_path) and os.path.getsize(html_path) > 0:
                checkdic.loc[m, 'd3flux'] = 1
            else:
                checkdic.loc[m, 'errorstring'] += f' {folder_name} d3flux pathway visualization failed;'

    d3flux_required = [m for m in d3flux_methods if m in checkdic.index]

    if all(checkdic['target_prediction'] == 1):
        if len(d3flux_required) == 0 or all(checkdic.loc[d3flux_required, 'd3flux'] == 1):
            overall_status = 'success'
            errorstring = ''
        else:
            overall_status = 'partial success'
            errorstring = checkdic.loc[d3flux_required, 'errorstring'].sum()
    elif all(checkdic['target_prediction'] == 0):
        overall_status = 'failed'
        errorstring = checkdic['errorstring'].sum()
    else:
        overall_status = 'partial success'
        errorstring = checkdic['errorstring'].sum()

    checkdic2 = checkdic.to_dict(orient='index')

    combined_output_path = os.path.join(path_results2, 'output.json')
    if os.path.isfile(combined_output_path):
        with open(combined_output_path, 'r') as f:
            output_data = json.load(f)

        for m_key, m_info in checkdic2.items():
            folder_name = taskidmap.get(m_key, m_key)
            if m_info['target_prediction'] == 0:
                output_data = remove_path(output_data, folder_name, "escher")
            if m_key in d3flux_methods and m_info['d3flux'] == 0:
                output_data = remove_path(output_data, folder_name, "d3flux")

        output_data['overall_status'] = overall_status
        output_data['errorstring'] = errorstring
        output_data['error_detail'] = checkdic2
        with open(combined_output_path, 'w') as f:
            json.dump(output_data, f, indent=4, ensure_ascii=False)


# python summary.py '/hpcfs/fhome/xuwenqi/project/OptMetarget/INPUT/model' '/hpcfs/fhome/xuwenqi/project/OptMetarget/INPUT/task' '/hpcfs/fhome/xuwenqi/project/OptMetarget/INPUT/map' '/hpcfs/fhome/xuwenqi/project/OptMetarget/OUTPUT' 'wangry@tib.cas.cn_20241203-092527_iML1515'
