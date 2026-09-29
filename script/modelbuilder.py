# coding:utf-8
"""
Author: Jingyi Cai (2024-2026)
Function: Build irreversible and isoenzyme-split models, assemble enzyme-constrained Pyomo models, and solve them.
Input: Called by target.py, EToptme.py, and summary.py. Receives a COBRA model plus reaction ΔG, metabolite lnC, and kcat/MW tables. Model_Solve reads OPTME_PYOMO_SOLVER and OPTME_SOLVER_OPTIONS.
Output: Returns solved models and flux dictionaries to the caller. Writes a temporary irreversible JSON under output/tmp/ while building the enzyme model.
"""
import cobra
import json
import os
import re
import numpy as np
import pandas as pd
import pyomo.environ as pyo
from pyomo.environ import *
from cobra.core import Reaction
from cobra.util.solver import set_objective


def get_dictionarymodel_data2(dictionarymodel,Concretemodel_Need_Data,kinetics_reactions):
    kcat_dict={}
    mw_dict={}
    enzyme_rxns_dict={}
    e0=[]
    totle_E=[]
    for enz in dictionarymodel['enzyme']:
        flag_dict=dictionarymodel['enzyme'][enz]
        if flag_dict['MW'] !='null':
            mw_dict[enz]=flag_dict['MW']
        enzyme_rxns_dict[enz]=list(flag_dict['reactions'].keys())
        for rxn in flag_dict['reactions']:
            if flag_dict['reactions'][rxn]['kcat'] != {}:
                kcat_dict[rxn]=flag_dict['reactions'][rxn]['kcat']
    e0 = list(kcat_dict.keys())
    totle_E = list(enzyme_rxns_dict.keys())   
    kmapp={}
    for enz in dictionarymodel['enzyme']:
        for rea in dictionarymodel['enzyme'][enz]['reactions']:
            if dictionarymodel['enzyme'][enz]['reactions'][rea]['km'] !={}:
                kmapp[rea]=dictionarymodel['enzyme'][enz]['reactions'][rea]['km'] 
    Concretemodel_Need_Data['kcat_dict']=kcat_dict
    Concretemodel_Need_Data['mw_dict']=mw_dict
    Concretemodel_Need_Data['enzyme_rxns_dict']=enzyme_rxns_dict
    Concretemodel_Need_Data['e0']=e0
    Concretemodel_Need_Data['totle_E']=totle_E
    kmapp2={}
    allsubstr=[]
    allsubstr2=[]
    P_lb={}
    P_ub={}
    for i in kinetics_reactions:
        if i in dictionarymodel['reactions'] and i in kmapp and i in kcat_dict:
            kmapp2[i]=kmapp[i]
        else:
            print(i)
    for i in kmapp2:
        P_lb[i]=0
        P_ub[i]=1
        if len(kmapp2[i])>2:
            print(i,kmapp2[i])
        if len(kmapp2[i])==2:
            ab =list(kmapp2[i].keys())
            allsubstr2.append(ab[0]+'_'+ab[1])
        for j in kmapp2[i]:
            allsubstr.append(j)
    allsubstr = list(set(allsubstr))
    allsubstr2 = list(set(allsubstr2))    
    Concretemodel_Need_Data['kmapp']=kmapp2
    Concretemodel_Need_Data['allsubstr']=allsubstr
    Concretemodel_Need_Data['allsubstr2']=allsubstr2
    Concretemodel_Need_Data['P_lb']=P_lb
    Concretemodel_Need_Data['P_ub']=P_ub


def get_data(model):
    reaction_list = []
    metabolite_list = []
    lb_list = {}
    ub_list = {}
    coef_matrix = {}

    # Pyomo S/SE/SET modes use NonNegativeReals. Reactions with ub<=0 and lb<0 only
    # allow negative flux in COBRA; re-orient to forward (v>=0): new_lb=-ub, new_ub=-lb,
    # and flip stoichiometric coefficients.
    for rea in model.reactions:
        reaction_list.append(rea.id)
        lb = rea.lower_bound
        ub = rea.upper_bound
        flip_negative_range = ub <= 0 and lb < 0
        if flip_negative_range:
            lb, ub = -ub, -lb

        lb_list[rea.id] = lb
        ub_list[rea.id] = ub

        for met in rea.metabolites:
            metabolite_list.append(met.id)
            if (met.id, rea.id) not in coef_matrix:
                coeff = rea.get_coefficient(met.id)
                if flip_negative_range:
                    coeff = -coeff
                coef_matrix[(met.id, rea.id)] = coeff

    metabolite_list = list(set(metabolite_list))

    return reaction_list, metabolite_list, lb_list, ub_list, coef_matrix


def json_load(path: str) :
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: str ~ The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary



# def Get_Concretemodel_Need_Data(reaction_g0_file,metabolites_lnC_file,model_file,reaction_kcat_MW_file):
#     Concretemodel_Need_Data={}
#     reaction_g0=pd.read_csv(reaction_g0_file,index_col=0,sep='\t')
#     Concretemodel_Need_Data['reaction_g0']=reaction_g0
#     metabolites_lnC = pd.read_csv(metabolites_lnC_file, index_col=0, sep='\t')
#     Concretemodel_Need_Data['metabolites_lnC']=metabolites_lnC
#     if re.search('\.xml',model_file):
#         model = cobra.io.read_sbml_model(model_file)
#     elif re.search('\.json',model_file):
#         model = cobra.io.json.load_json_model(model_file)
#         dicmodel=json_load(model_file)
#     #convert_to_irreversible(model)
#     reaction_kcat_MW=pd.read_csv(reaction_kcat_MW_file,index_col=0)
#     Concretemodel_Need_Data['model']=model
#     Concretemodel_Need_Data['reaction_kcat_MW']=reaction_kcat_MW
#     [reaction_list,metabolite_list,lb_list,ub_list,coef_matrix]=Get_Model_Data(model)
#     Concretemodel_Need_Data['reaction_list']=reaction_list
#     Concretemodel_Need_Data['metabolite_list']=metabolite_list
#     Concretemodel_Need_Data['lb_list']=lb_list
#     Concretemodel_Need_Data['ub_list']=ub_list
#     Concretemodel_Need_Data['coef_matrix']=coef_matrix

#     def get_gene_reaction_dict(dicmodel):
#         reaction_gene_dict = {}
#         for reaction in dicmodel["reactions"]:
#             if reaction["gene_reaction_rule"] != "":
#                 reaction_gene_dict[reaction["id"]] = reaction["gene_reaction_rule"].replace(" ","_")
#         gene_reaction_dict = {}
#         for key, value in reaction_gene_dict.items():
#                 if value not in gene_reaction_dict:
#                     gene_reaction_dict[value] = [key]
#                 else:
#                     gene_reaction_dict[value].append(key)
#         return gene_reaction_dict

#     gene_reaction_dict=get_gene_reaction_dict(dicmodel)
#     totle_E=list(gene_reaction_dict.keys())
#     e0=list(reaction_kcat_MW.index)
#     Concretemodel_Need_Data['totle_E']=totle_E
#     Concretemodel_Need_Data['gene_reaction_dict']=gene_reaction_dict
#     Concretemodel_Need_Data['e0']=e0
    
#     theMW={}

#     for i in dicmodel['reactions']:
#         # print(type(i['kcat']))
#         if i['kcat']!='' and i['kcat_MW']!='':
#             theMW[i['gene_reaction_rule']]=(i['kcat']/i['kcat_MW'])

#     Concretemodel_Need_Data['theMW']=theMW
    
#     return Concretemodel_Need_Data,model,dicmodel  #,e0


def remove_unused_metabolites(model):
    metabolite_used_list = []
    for reaction in model.reactions:
        for metabolite in reaction.metabolites:
            metabolite_used_list.append(metabolite.id)
            metabolite_used_list = list(set(metabolite_used_list))
    metabolite_unused_num = 1
    while(metabolite_unused_num > 0):
        metabolite_unused_num = 0
        for metabolite in model.metabolites:
            if metabolite.id not in metabolite_used_list:
                metabolite_unused_num += 1
                model.remove_metabolites([metabolite])
    return model


def isoenzyme_split(model):
    """Split isoenzyme reaction to mutiple reaction

    Arguments
    ----------
    * model: cobra.Model.
    
    :return: new cobra.Model.
    """  
    for r in model.reactions:
        if re.search(" or ", r.gene_reaction_rule):
            rea = r.copy()
            gene = r.gene_reaction_rule.split(" or ")
            for index, value in enumerate(gene):
                if index == 0:
                    r.id = r.id + "_num1"
                    r.gene_reaction_rule = value
                else:
                    r_add = rea.copy()
                    r_add.id = rea.id + "_num" + str(index+1)
                    r_add.gene_reaction_rule = value
                    model.add_reactions([r_add])
    for r in model.reactions:
        r.gene_reaction_rule = r.gene_reaction_rule.strip("( )")
    return model

# dictionary_model:  json_load(model_file)


def trans_model2standard_json_etgem2(model,model_file):
    model_name = model_file.split('/')[-1].split('.')[0]
    os.makedirs(os.path.join("output", "tmp"), exist_ok=True)
    json_path = os.path.join("output", "tmp", "%s_irreversible.json" % model_name)
    print(json_path)
    cobra.io.save_json_model(model, json_path)
    dictionary_model = json_load(json_path)
    for eachreaction in dictionary_model['reactions']:
        eachreaction["g0"] = -99999
    for eachmetabolite in dictionary_model['metabolites']:
        if eachmetabolite["id"][-1] != 'c':
            eachmetabolite["concentration_ub"] = 0.0
            eachmetabolite["concentration_lb"] = 0.0
        elif eachmetabolite["id"][-1] == 'c':
            eachmetabolite["concentration_ub"] = -3.912023005
            eachmetabolite["concentration_lb"] = -14.50865774

    reaction_gene_dict = {}
    for reaction in dictionary_model["reactions"]:
        if reaction["gene_reaction_rule"] != "":
            flag1=reaction["gene_reaction_rule"].split(' and ')
            flag1.sort(key=None, reverse=False)
            flag2=' and '.join(flag1)
            reaction_gene_dict[reaction["id"]] = flag2
    gene_reaction_dict = {}
    for key, value in reaction_gene_dict.items():
            if value not in gene_reaction_dict:
                gene_reaction_dict[value] = [key]
            else:
                gene_reaction_dict[value].append(key)
    enzyme_dict = {}
    for key,value in gene_reaction_dict.items():
        new_enzyme = {}
        new_enzyme["gpr"] = key
        gene_list = key.split("_and_")
        enzyme_gene_dict = {}
        for gene in gene_list:
            enzyme_gene_dict[gene] = 1
        new_enzyme["genes"] = enzyme_gene_dict
        new_enzyme["MW"] = "null"
        reaction_dict = {}
        for reaction in value:
            
            single_reaction_dict = {}
            single_reaction_dict["kcat"] = {}
            single_reaction_dict["km"] = {}
            # single_reaction_dict["ki"] = "null"
            # single_reaction_dict["dGr0"] = get_g0_by_id(dictionary_model,reaction)
            # if "annotation" in get_reaction_by_id(dictionary_model,reaction).keys():
            #     if "ec-code" in get_reaction_by_id(dictionary_model,reaction)["annotation"].keys():
            #         single_reaction_dict["EC"] = get_reaction_by_id(dictionary_model,reaction)["annotation"]["ec-code"]
            # single_reaction_dict["SA"] = "null"
            # single_reaction_dict["cofactor"] = "null"
            # single_reaction_dict["inhibitor"] = "null"
            # single_reaction_dict["ACTIVATING_COMPOUND"] = "null"
            reaction_dict[reaction] = single_reaction_dict
        new_enzyme["reactions"] = reaction_dict
        enzyme_dict[key] = new_enzyme
    dictionary_model["enzyme"] = enzyme_dict
    dictionary_model['metabolites']
    
    dict_met={}
    for eachmetabolite in dictionary_model['metabolites']:
        dict_met[eachmetabolite["id"]]=eachmetabolite
    dictionary_model['metabolites']=dict_met 
    dict_rea={}
    for eachmetabolite in dictionary_model['reactions']:
        dict_rea[eachmetabolite["id"]]=eachmetabolite
    dictionary_model['reactions']=dict_rea 
    dict_gene={}
    for eachmetabolite in dictionary_model['genes']:
        dict_gene[eachmetabolite["id"]]=eachmetabolite
    dictionary_model['genes']=dict_gene     
    reactions_reaction={}
    for eachmetabolite in dictionary_model['reactions']:
        dictionary_model['reactions'][eachmetabolite]['reaction']=model.reactions.get_by_id(eachmetabolite).reaction
    return dictionary_model,model


def convert_to_irreversible(model):
    """Split reversible reactions into two irreversible reactions

    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.

    Arguments
    ----------
    * model: cobra.Model ~ A Model object which will be modified in place.

    """
    #warn("deprecated, not applicable for optlang solvers", DeprecationWarning)
    reactions_to_add = []
    coefficients = {}
    for reaction in model.reactions:
        # If a reaction is reverse only, the forward reaction (which
        # will be constrained to 0) will be left in the model.
        if reaction.lower_bound < 0 and reaction.upper_bound > 0:
            reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction.name=reaction.name
            reverse_reaction.lower_bound = max(0, -reaction.upper_bound)
            reverse_reaction.upper_bound = -reaction.lower_bound
            coefficients[
                reverse_reaction] = reaction.objective_coefficient * -1
            reaction.lower_bound = max(0, reaction.lower_bound)
            reaction.upper_bound = max(0, reaction.upper_bound)
            # Make the directions aware of each other
            reverse_reaction.notes = reaction.notes
            reverse_reaction.annotation = reaction.annotation
            reaction.notes["reflection"] = reverse_reaction.id
            reverse_reaction.notes["reflection"] = reaction.id
            reaction_dict = {k: v * -1
                             for k, v in reaction._metabolites.items()}
            reverse_reaction.add_metabolites(reaction_dict)
            reverse_reaction._model = reaction._model
            reverse_reaction._genes = reaction._genes
            for gene in reaction._genes:
                gene._reaction.add(reverse_reaction)
            reverse_reaction.subsystem = reaction.subsystem
            reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
            
    model.add_reactions(reactions_to_add)
    set_objective(model, coefficients, additive=True)
    return model
    
#Encapsulating parameters used in Concretemodel


def Get_Concretemodel_Need_Data2(model):
    Concretemodel_Need_Data={}
    Concretemodel_Need_Data['model']=model
    [reaction_list,metabolite_list,lb_list,ub_list,coef_matrix]=get_data(model)
    Concretemodel_Need_Data['reaction_list']=reaction_list
    Concretemodel_Need_Data['metabolite_list']=metabolite_list
    Concretemodel_Need_Data['lb_list']=lb_list
    Concretemodel_Need_Data['ub_list']=ub_list
    Concretemodel_Need_Data['coef_matrix']=coef_matrix
    return (Concretemodel_Need_Data)
#set_obj_value,set_metabolite,set_Df: only 'True' and 'False'
    #  kcat_dict:   reaction:value    
    #  mw_dict  :   enzyme:value
    #  enzyme_rxns_dict:     enzyme:rxns_list


def FBA_template2(coef_matrix=None,metabolites_lnC=None,reaction_g0=None,lb_list=None,ub_list=None,reaction_kcat_MW=None,B_value=None,kcat_dict=None,mw_dict=None,K_value=None,\
            product_value=None,biomass_value=None,substrate_value=None,E_total=None,\
            product_name=None,reaction_list=None,metabolite_list=None,enzyme_rxns_dict=None,e0=None,totle_E=None,obj_name=None,obj_target=None,\
            biomass_id=None,substrate_name=None,Bottleneck_reaction_list=None,biomass_name=None,\
            set_obj_B_value=False,set_obj_sum_e=False,set_obj_value_e=False,set_obj_E_value=False,set_obj_single_E_value=False,set_obj_X_value=False,\
            set_obj_value=False,set_obj_V_value=False,set_obj_TM_value=False,set_obj_Met_value=False,set_obj_Two_value=False,\
            set_metabolite=False,set_Df=False,set_thermodynamics0=False,set_thermodynamics=False,set_integer=False,set_fix_sum_e=False,\
            set_enzyme_constraint=False,set_substrate_ini=False, set_biomass_ini=False,set_product_ini=False,set_reactions_ini=False,\
            set_metabolite_ratio=False,set_Bottleneck_reaction=False,set_Df_value=False,set_enzyme_value=False,set_max_sum_e=False,\
            mode=None,constr_coeff=None,Concretemodel_Need_Data=None):
    
    """According to the parameter conditions provided by the user, the specific pyomo model is returned.

    Notes
    ----------
    * reaction_list: List of reaction IDs for the model.
    It is obtained by analyzing SBML model with get_data(model) function.
    
    * metabolite_list: List of metabolite IDs for the model.
    It is obtained by analyzing SBML model with get_data(model) function.
    
    * coef_matrix: The model coefficient matrix.
    It is obtained by analyzing SBML model with get_data(model) function.
    
    * metabolites_lnC:Upper and lower bound information of metabolite concentration (Natural logarithm).
    The format is as follows (in .txt file):
    id	lnClb	lnCub
    2pg_c	-14.50865774	-3.912023005
    13dpg_c	-14.50865774	-3.912023005
    ...
    
    * reaction_g0: Thermodynamic parameter of reactions (drG'¡ã) .
    The format is as follows (in .txt file):
    reaction	g0
    13PPDH2	-21.3
    13PPDH2_reverse	21.3
    ...

    * reaction_kcat_MW: The enzymatic data of kcat (turnover number, /h) divided by MW (molecular weight, kDa).
    The format is as follows (in .csv file):
    kcat,MW,kcat_MW
    AADDGT	231887.7737	40.6396	5705.956104
    AAMYL	23490.41652	56.63940007	414.736323
    ...
    
    * lb_list: Lower bound of reaction flux. It is obtained by analyzing SBML model with get_data(model) function.
    * ub_list: Upper bound of reaction flux. It is obtained by analyzing SBML model with get_data(model) function.
    * obj_name: Name of object, such as set_obj_value, set_obj_single_E_value, set_obj_TM_value and set_obj_Met_value.    
    * K_value: the maximum value minus the minimum value of reaction thermodynamic driving force (1249, "max(DFi,max)-min(DFi,min)").
    * obj_target: Type of object function (maximize or minimize).
    * set_obj_value: Set the flux as the object function (True or False).
    * set_substrate_ini: Adding substrate amount constraints (True or False).
    * substrate_name: Substrate input reaction ID in the model (such as EX_glc__D_e_reverse).
    * substrate_value: Set the upper bound for substrate input reaction flux (10 mmol/h/gDW).
    * set_biomass_ini: Set the lower bound for biomass synthesis reaction flux (True or False).
    * biomass_value: The lower bound of biomass synthesis reaction flux.
    * biomass_id: Biomass synthesis reaction ID in the model (BIOMASS_Ec_iML1515_core_75p37M).
    * set_metabolite: Set the upper and lower bounds of metabolite concentration (True or False).
    * set_Df: Adding thermodynamic driving force expression for reactions (True or False).
    * set_obj_B_value: The object function is the maximizing thermodynamic driving force of a pathway (True or False)
    * set_stoi_matrix: Adding flux balance constraints (True or False).
    * set_bound: Set the upper and lower bounds of reaction flux (True or False).
    * set_enzyme_constraint: Adding enzymamic constraints (True or False).
    * set_integer: Adding binary variables constraints (True or False).
    * set_metabolite_ratio: Adding concentration ratio constraints for metabolites (True or False).
    * set_thermodynamics: Adding thermodynamic constraints (True or False).
    * B_value: The value of maximizing the minimum thermodynamic driving force (MDF).
    * set_obj_E_value: The object function is the minimum enzyme cost of a pathway (True or False).
    * set_obj_V_value: The object function is the pFBA of a pathway (True or False)
    * set_obj_TM_value: The object function is the thermodynamic driving force of a reaction (True or False).
    * set_obj_Met_value: The object function is the concentration of a metabolite (True or False).
    * set_obj_single_E_value: The object function is the enzyme cost of a reaction (True or False).     
    * E_total: Total amount constraint of effective enzymes pool (0.13).
    * Bottleneck_reaction_list: A list extracted from the result file automatically.
    * set_Bottleneck_reaction: Adding integer variable constraints for specific reaction (True or False).
    """
    Concretemodel = ConcreteModel() 
    
    reaction_list=Concretemodel_Need_Data['reaction_list']
    metabolite_list=Concretemodel_Need_Data['metabolite_list']
    coef_matrix=Concretemodel_Need_Data['coef_matrix']
    lb_list=Concretemodel_Need_Data['lb_list']
    ub_list=Concretemodel_Need_Data['ub_list']
    if   mode=='S':
        if set_obj_X_value:
            set_obj_Two_value=True
    if   mode=='ST':
        if set_obj_B_value:
            set_thermodynamics0=True
        if not set_obj_B_value:
            set_thermodynamics=True
        set_Df_value=True
        set_metabolite=True
        set_Df=True
        set_integer=True
        set_metabolite_ratio=True
        
        metabolites_lnC=Concretemodel_Need_Data['metabolites_lnC']
        reaction_g0=Concretemodel_Need_Data['reaction_g0']
        B_value=Concretemodel_Need_Data['B_value']
        K_value=Concretemodel_Need_Data['K_value']
        
        
    if  mode=='SE': 
        set_enzyme_constraint=True
        set_enzyme_value=True
        set_max_sum_e=True
        # if set_obj_sum_e or set_fix_sum_e:
        #     set_max_sum_e=False
        # if not set_obj_sum_e:
        #     set_max_sum_e=True
            
        kcat_dict=Concretemodel_Need_Data['kcat_dict']
        mw_dict=Concretemodel_Need_Data['mw_dict']
        enzyme_rxns_dict=Concretemodel_Need_Data['enzyme_rxns_dict']
        e0=Concretemodel_Need_Data['e0']
        totle_E=Concretemodel_Need_Data['totle_E']
        E_total=Concretemodel_Need_Data['E_total']
           
    if  mode=='SET':   
        if set_obj_X_value:
            set_obj_Two_value=True 
        if set_obj_B_value:
            set_thermodynamics0=True
        if not set_obj_B_value:
            set_thermodynamics=True
        if set_obj_sum_e or set_fix_sum_e:
            set_max_sum_e=False
        if not set_obj_sum_e:
            set_max_sum_e=True
        
        set_stoi_matrix=True
        set_bound=True   
        set_Df_value=True
        set_enzyme_value=True
        set_metabolite=True
        set_Df=True
        set_metabolite_ratio=True
        set_integer=True 
        set_enzyme_constraint=True    
        kcat_dict=Concretemodel_Need_Data['kcat_dict']
        mw_dict=Concretemodel_Need_Data['mw_dict']
        enzyme_rxns_dict=Concretemodel_Need_Data['enzyme_rxns_dict']
        e0=Concretemodel_Need_Data['e0']
        totle_E=Concretemodel_Need_Data['totle_E']
        E_total=Concretemodel_Need_Data['E_total']   
        metabolites_lnC=Concretemodel_Need_Data['metabolites_lnC']
        reaction_g0=Concretemodel_Need_Data['reaction_g0']
        B_value=Concretemodel_Need_Data['B_value']
        K_value=Concretemodel_Need_Data['K_value']        
                    
    if constr_coeff !=None:
        if 'substrate_constrain' in constr_coeff:
            set_substrate_ini=True
            substrate_value=constr_coeff['substrate_constrain'][1]
            substrate_name=constr_coeff['substrate_constrain'][0]
        if 'biomass_constrain' in constr_coeff:
            set_biomass_ini=True
            biomass_value=constr_coeff['biomass_constrain'][1]
            biomass_id=constr_coeff['biomass_constrain'][0]
        if 'product_constrain' in constr_coeff:
            set_product_ini=True
            product_value=constr_coeff['product_constrain'][1]
            product_name=constr_coeff['product_constrain'][0]
        if 'fix_reactions' in constr_coeff and constr_coeff['fix_reactions'] !={}:
            set_reactions_ini=True
        if 'fix_E_total' in constr_coeff:
            set_fix_sum_e=True
            set_max_sum_e=False
        # if 'fix_E_total' in constr_coeff:
        #     set_fix_sum_e=True
        #     fix_sum_e=constr_coeff['fix_E_total']
        #     set_max_sum_e=False          
    if mode =='S0':
        Concretemodel.reaction = pyo.Var(reaction_list)
    elif mode =='S' or mode =='SE' or mode =='SET' or mode =='ST':
           # Concretemodel.metabolite = pyo.Var(metabolite_list,  within=Reals)  
        Concretemodel.reaction = pyo.Var(reaction_list,  within=NonNegativeReals,bounds=lambda m, i: (lb_list[i], ub_list[i]))
    # Concretemodel.reaction = pyo.Var(reaction_list,  within=NonNegativeReals,bounds=lambda m, i: (lb_list[i], ub_list[i])) # reaction flux variable  

    # Concretemodel.X = pyo.Var() 
    # SV=0
    # Concretemodel.set_stoi_matrix = ConstraintList()
    # for i in metabolite_list:
    #     Concretemodel.set_stoi_matrix.add(sum(coef_matrix[i,j]*Concretemodel.reaction[j]  for j in reaction_list if (i,j) in coef_matrix.keys() )==0)
    def set_stoi_matrix(m,i):
        return sum(coef_matrix[i,j]*m.reaction[j]  for j in reaction_list if (i,j) in coef_matrix.keys() )==0
    Concretemodel.set_stoi_matrix = Constraint( metabolite_list,rule=set_stoi_matrix)
    # LB¡¢UB

    def set_bound(m,j):
        return inequality(lb_list[j],m.reaction[j],ub_list[j])
    Concretemodel.set_bound = Constraint(reaction_list,rule=set_bound)   
    if set_obj_Two_value:
        # def set_obj_Two_value(m,obj_name,biomass_name):
        #     return m.X==(m.reaction[obj_name]*m.reaction[biomass_name])
        Concretemodel.set_obj_Two_value = Constraint(expr=Concretemodel.reaction[obj_name]*Concretemodel.reaction[biomass_name]==Concretemodel.X)

    if reaction_kcat_MW is not None:
    # check if kcat_MW is a column in reaction_kcat_MW
        if 'kcat_MW' in reaction_kcat_MW.columns:
            kcatmw='kcat_MW'
        elif 'kcat/mw' in reaction_kcat_MW.columns:
            kcatmw='kcat/mw'
            
    #           ****ÉèÖÃ±äÁ¿****    
    # ÈÈÁ¦Ñ§±äÁ¿ÉèÖÃ    
    if set_Df_value:
        reaction_listDF=[j for j in reaction_g0.index if j in reaction_list]
        ConMet_list= list(set(metabolite_list).intersection(set(metabolites_lnC.index)))
        Concretemodel.metabolite = pyo.Var(ConMet_list,  within=Reals) # metabolite concentration variable
        Concretemodel.Df = pyo.Var(reaction_listDF,  within=Reals) # thermodynamic driving force variable--reactions
        Concretemodel.z = pyo.Var(reaction_listDF,  within=pyo.Binary)    # binary variable
        Concretemodel.B = pyo.Var()     # thermodynamic driving force variable--mdfz
 
    # µ°°×ÎªÖÐÐÄÃ¸±äÁ¿ÉèÖÃ
    if set_enzyme_value:
        Concretemodel.e = pyo.Var(e0, within=NonNegativeReals)
        Concretemodel.e1 = pyo.Var(totle_E, within=NonNegativeReals)              

    #           ****ÉèÖÃÄ¿±êº¯Êý****
    
    # ÒÔB_valueÎªÄ¿±êº¯Êý
    if set_obj_X_value:
        def set_obj_X_value(m):
            return m.X
        Concretemodel.obj = Objective(rule=set_obj_X_value, sense=maximize) 
    if set_obj_B_value:             
        def set_obj_B_value(m):
            return m.B  
        Concretemodel.obj = Objective(rule=set_obj_B_value, sense=maximize)  
    # ÒÔ×Üµ°°×Á¿¼ÓºÍÎªÄ¿±êº¯Êý
    if set_obj_sum_e:             
        def set_obj_sum_e(m):
            return sum(Concretemodel.e1[i] * mw_dict[i] for i in totle_E if i in mw_dict)
        if obj_target=='maximize':
            Concretemodel.obj = Objective(rule=set_obj_sum_e, sense=maximize)
        elif obj_target=='minimize':
            Concretemodel.obj = Objective(rule=set_obj_sum_e, sense=minimize)
    # ÒÔ×îÐ¡»¯Ã¸Á¿¼ÓºÍÎªÄ¿±êº¯Êý£¨²»ÊÇÒÔµ°°×ÎªÖÐÐÄ£©   
    if set_obj_E_value:             
        def set_obj_E_value(m):
            return sum(m.reaction[j]/(reaction_kcat_MW.loc[j,kcatmw]) for j in reaction_kcat_MW.index if j in reaction_list)
        Concretemodel.obj = Objective(rule=set_obj_E_value, sense=minimize)      
    # ÒÔµ¥Ò»µ°°×ÎªÄ¿±êº¯Êý
    if set_obj_value_e:   
        def set_obj_value_e(m):
            return m.e1[obj_name]*mw_dict[obj_name]
        if obj_target=='maximize':
            Concretemodel.obj = Objective(rule=set_obj_value_e, sense=maximize)
        elif obj_target=='minimize':
            Concretemodel.obj = Objective(rule=set_obj_value_e, sense=minimize)    
    # ÒÔµ¥Ò»µ°°×ÎªÄ¿±êº¯Êý£¨²»ÊÇÒÔµ°°×ÎªÖÐÐÄ£©
    if set_obj_single_E_value:             
        def set_obj_single_E_value(m):
            return m.reaction[obj_name]/(reaction_kcat_MW.loc[obj_name,kcatmw])
        if obj_target=='maximize':
            Concretemodel.obj = Objective(rule=set_obj_single_E_value, sense=maximize) 
        elif obj_target=='minimize':
            Concretemodel.obj = Objective(rule=set_obj_single_E_value, sense=minimize)
    # ÒÔÄ³¸ö·´Ó¦ÎªÄ¿±êº¯Êý
    if set_obj_value:   
        def set_obj_value(m):
            return m.reaction[obj_name]
        if obj_target=='maximize':
            Concretemodel.obj = Objective(rule=set_obj_value, sense=maximize)
        elif obj_target=='minimize':
            Concretemodel.obj = Objective(rule=set_obj_value, sense=minimize)
    # 2个变量
    # if set_obj_Two_value:   
    #     def set_obj_Two_value(m):
    #         return m.reaction[obj_name]*m.reaction[biomass_name]
    #     if obj_target=='maximize':
    #         Concretemodel.obj = Objective(rule=set_obj_Two_value, sense=maximize)
    #     elif obj_target=='minimize':
    #         Concretemodel.obj = Objective(rule=set_obj_Two_value, sense=minimize)    
    #×îÐ¡»¯Í¨Á¿ºÍÎªÄ¿±êº¯Êý (pFBA)
    if set_obj_V_value:             
        def set_obj_V_value(m):
            return sum(m.reaction[j] for j in reaction_list)
        Concretemodel.obj = Objective(rule=set_obj_V_value, sense=minimize)  
    #To calculate the variability of thermodynamic driving force (Dfi) of single reaction.
    if set_obj_TM_value:   
        def set_obj_TM_value(m):
            return m.Df[obj_name]
        if obj_target=='maximize':
            Concretemodel.obj = Objective(rule=set_obj_TM_value, sense=maximize)
        elif obj_target=='minimize':
            Concretemodel.obj = Objective(rule=set_obj_TM_value, sense=minimize)
    #ÒÔ´úÐ»ÎïÅ¨¶ÈÎªÄ¿±êº¯Êý£¨ÈÈÁ¦Ñ§£©
    if set_obj_Met_value:   
        def set_obj_Met_value(m):
            return m.metabolite[obj_name]
        if obj_target=='maximize':
            Concretemodel.obj = Objective(rule=set_obj_Met_value, sense=maximize)
        elif obj_target=='minimize':
            Concretemodel.obj = Objective(rule=set_obj_Met_value, sense=minimize)
 
 
 
 
    #           ****Ìí¼ÓÔ¼Êø****
    
    #  ÈÈÁ¦Ñ§Ô¼Êøset_metabolite¡¢set_Df¡¢set_obj_B_value£¨B_valueÎª±äÁ¿£©¡¢set_thermodynamics£¨B_valueÎª¶¨Öµ£©¡¢set_integer
    if set_metabolite:
        def set_metabolite(m,i):
            return  inequality(metabolites_lnC.loc[i,'lnClb'], m.metabolite[i], metabolites_lnC.loc[i,'lnCub'])
        Concretemodel.set_metabolite= Constraint(ConMet_list,rule=set_metabolite)        
    #thermodynamic driving force expression for reactions
    if set_Df:
        def set_Df(m,j):
            return  m.Df[j]==-reaction_g0.loc[j,'g0']-2.579*sum(coef_matrix[i,j]*m.metabolite[i]  for i in ConMet_list if (i,j) in coef_matrix.keys())
        Concretemodel.set_Df = Constraint(reaction_listDF,rule=set_Df)
    #Adding thermodynamic MDF(B) object function
    if set_thermodynamics0:
        def set_thermodynamics0(m,j):
            return m.B<=(m.Df[j]+(1-m.z[j])*K_value)
        Concretemodel.set_obj_B_value = Constraint(reaction_listDF, rule=set_thermodynamics0)
    #Adding thermodynamic constraints
    if set_thermodynamics:
        def set_thermodynamics(m,j):
            return (m.Df[j]+(1-m.z[j])*K_value)>= B_value
        Concretemodel.set_thermodynamics = Constraint(reaction_listDF, rule=set_thermodynamics)
    #Adding binary variables constraints
    if set_integer:
        def set_integer(m,j):
            return m.reaction[j]<=m.z[j]*ub_list[j] 
        Concretemodel.set_integer = Constraint(reaction_listDF,rule=set_integer)     
        
        
    # Ã¸Ô¼Êø
    if set_enzyme_constraint:
        # flux- enzyme relationship
        Concretemodel.set_constr2 = ConstraintList()
        for i in e0:
            Concretemodel.set_constr2.add(Concretemodel.e[i] *kcat_dict[i] >= Concretemodel.reaction[i])
        # protein  concentration
        Concretemodel.set_constr3 = ConstraintList()  
        for i in totle_E:   
            Concretemodel.set_constr3.add(Concretemodel.e1[i] == sum(Concretemodel.e[j] for j in enzyme_rxns_dict[i] if j in e0))
        # total protein concentration
    if set_max_sum_e:
        Concretemodel.set_max_sum_e = Constraint(expr=sum(Concretemodel.e1[i] * mw_dict[i] for i in totle_E if i in mw_dict) <= E_total)
    if set_fix_sum_e:
        Concretemodel.set_max_sum_e = Constraint(expr=sum(Concretemodel.e1[i] * mw_dict[i] for i in totle_E if i in mw_dict) == E_total) 
        
    # ×î´óµ×ÎïÉãÈëÔ¼Êø
    if set_substrate_ini:
        def set_substrate_ini(m): 
            return m.reaction[substrate_name] <= substrate_value
        Concretemodel.set_substrate_ini = Constraint(rule=set_substrate_ini)       
    # ¹Ì¶¨biomass·¶Î§
    if set_biomass_ini:
        def set_biomass_ini(m): 
            return m.reaction[biomass_id] >=biomass_value
        Concretemodel.set_biomass_ini = Constraint(rule=set_biomass_ini)
    # ¹Ì¶¨product·¶Î§
    if set_product_ini:
        Concretemodel.set_product_ini = Constraint(expr=Concretemodel.reaction[product_name] >= product_value)
    # ÅúÁ¿¹Ì¶¨·´Ó¦Í¨Á¿
    if set_reactions_ini:
        Concretemodel.set_reactions_ini = ConstraintList()
        for i in constr_coeff['fix_reactions']:
            Concretemodel.set_reactions_ini.add(Concretemodel.reaction[i] >= constr_coeff['fix_reactions'][i][0])
            Concretemodel.set_reactions_ini.add(Concretemodel.reaction[i] <= constr_coeff['fix_reactions'][i][1])

    #Adding concentration ratio constraints for metabolites
    if set_metabolite_ratio:
        def set_atp_adp(m):
            return m.metabolite['atp_c']-m.metabolite['adp_c']==np.log(10)
        def set_adp_amp(m):
            return m.metabolite['adp_c']-m.metabolite['amp_c']==np.log(1)
        def set_nad_nadh(m):
            return m.metabolite['nad_c']-m.metabolite['nadh_c']==np.log(10)
        def set_nadph_nadp(m):
            return m.metabolite['nadph_c']-m.metabolite['nadp_c']==np.log(10)
        def set_hco3_co2(m):
            return m.metabolite['hco3_c']-m.metabolite['co2_c']==np.log(2)

        Concretemodel.set_atp_adp = Constraint(rule=set_atp_adp) 
        Concretemodel.set_adp_amp = Constraint(rule=set_adp_amp) 
        Concretemodel.set_nad_nadh = Constraint(rule=set_nad_nadh) 
        Concretemodel.set_nadph_nadp = Constraint(rule=set_nadph_nadp) 
        Concretemodel.set_hco3_co2 = Constraint(rule=set_hco3_co2)

    #Adding Bottleneck reaction constraints
    if set_Bottleneck_reaction:
        def set_Bottleneck_reaction(m,j):
            return m.z[j]==1 
        Concretemodel.set_Bottleneck_reaction = Constraint(Bottleneck_reaction_list,rule=set_Bottleneck_reaction) 

    return Concretemodel


def Model_Solve(model, solver=None):
    """Solve a Pyomo model with the user-selected solver and parameters.

    solver defaults to OPTME_PYOMO_SOLVER (then OPTME_COBRA_SOLVER).
    OPTME_SOLVER_OPTIONS is a JSON object of solver option names and values.
    """
    if not solver:
        solver = os.environ.get("OPTME_PYOMO_SOLVER") or os.environ.get("OPTME_COBRA_SOLVER") or "cplex"
    opt = pyo.SolverFactory(solver)
    raw = os.environ.get("OPTME_SOLVER_OPTIONS", "").strip()
    options = json.loads(raw) if raw else {}
    if options:
        opt.solve(model, options=options)
    else:
        opt.solve(model)
    return model

#Maximum growth rate calculation


def conbine_flux(solved_Concretemodel,modelid):   
    #get flux of stoichiometry from E_protain_FBA

    flux_efba={ i: value(solved_Concretemodel.reaction[i])   for i in solved_Concretemodel.reaction}
    reaction_list=[]
    for rea in modelid.reactions:
        reaction_list.append(rea.id)
    fluxnew={}
    for i in reaction_list:
        fluxnew[i]=0
        for j in flux_efba.keys():
            if i == j:
                fluxnew[i]=fluxnew[i]+flux_efba[i]
            if i+'_num' ==j[:-1]:
                fluxnew[i]=fluxnew[i]+flux_efba[j]
            if i +'_reverse'  == j or  i+'_reverse_num' == j[:-1]:
                fluxnew[i]=fluxnew[i]-flux_efba[j]
    return fluxnew


def Get_Concretemodel_Need_Data_g0(Concretemodel_Need_Data,reaction_g0_file,metabolites_lnC_file,reaction_kcat_MW_file):
    reaction_g0=pd.read_csv(reaction_g0_file,index_col=0,sep='\t')
    Concretemodel_Need_Data['reaction_g0']=reaction_g0
    metabolites_lnC = pd.read_csv(metabolites_lnC_file, index_col=0,sep='\t')
    Concretemodel_Need_Data['metabolites_lnC']=metabolites_lnC
    reaction_kcat_MW=pd.read_csv(reaction_kcat_MW_file,index_col=0)
    Concretemodel_Need_Data['reaction_kcat_MW']=reaction_kcat_MW


def showflux(ECM_EGFBA):
    flux_positive={i:"{:.3e}".format(value(ECM_EGFBA.reaction[i])) for i in ECM_EGFBA.reaction if value(ECM_EGFBA.reaction[i]) >0 and ('EX_' in i or 'DM_' in i)}
    return flux_positive


