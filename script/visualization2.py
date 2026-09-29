# coding:utf-8
"""
Author: Jingyi Cai (2024-2026)
Function: Draw central-metabolism flux maps from wild-type, overproduction, and fold-change fluxes.
Input: Called by target.py with flux dictionaries and INPUT/map/iJO1366.Central_metabolism.json.
Output: d3flux HTML and TSV files written into the method result folder by the caller.
"""
import sys
import os
import cobra
from cobra import Model, Reaction, Metabolite
import json
import re
import d3flux
import pandas as pd
from cobra.medium import find_boundary_types
# from bigg_function import *


def get_model_met_from_reaction_string(model,reaction_str,fwd_arrow,rev_arrow,reversible_arrow,term_split):
    # set the arrows
    forward_arrow_finder = (
            _forward_arrow_finder
            if fwd_arrow is None
            else re.compile(re.escape(fwd_arrow))
        )
    reverse_arrow_finder = (
            _reverse_arrow_finder
            if rev_arrow is None
            else re.compile(re.escape(rev_arrow))
        )
    reversible_arrow_finder = (
            _reversible_arrow_finder
            if reversible_arrow is None
            else re.compile(re.escape(reversible_arrow))
        )
    compartment_finder = re.compile("^\s*(\[[A-Za-z]\])\s*:*")
    found_compartments = compartment_finder.findall(reaction_str)
    if len(found_compartments) == 1:
        compartment = found_compartments[0]
        reaction_str = compartment_finder.sub("", reaction_str)
    else:
        compartment = ""
    # reversible case
    arrow_match = reversible_arrow_finder.search(reaction_str)
    if arrow_match is not None:
        pass
    else:  # irreversible
        # try forward
        arrow_match = forward_arrow_finder.search(reaction_str)
        if arrow_match is not None:
            pass
        else:
            # must be reverse
            arrow_match = reverse_arrow_finder.search(reaction_str)
            if arrow_match is None:
                raise ValueError("no suitable arrow found in '%s'" % reaction_str)
            else:
                pass
    reactant_str = reaction_str[: arrow_match.start()].strip()
    product_str = reaction_str[arrow_match.end() :].strip()
    for substr, factor in ((reactant_str, -1), (product_str, 1)): # (('asn-L_c + h2o_c', -1), ('NH4+ + asp-L_c', 1)) -1和1代表方程左边还是右边
        if len(substr) == 0:
            continue
        # substr： ’NH4+ + asp-L_c‘， term_split==’+‘，此时以+号分割就会把NH4+分割成NH4和空字符串，导致后面的met_id为空，从而导致错误
        # asn-L_c + h2o_c --> NH4+ + asp-L_c，方程里就是给的错误ID:NH4+,所以改掉cobra的内置函数————reaction.build_reaction_from_string(new_equ)，将其中分割的+号前加一个空格，即" +",nh4+就是一个整体了
        for term in substr.split(term_split):
            term = term.strip()
            if term.lower() == "nothing":
                continue
            if " " in term:
                num_str, met_id = term.split()
                num = float(num_str.lstrip("(").rstrip(")")) * factor
            else:
                met_id = term
                num = factor
            met_id += compartment  
            try:
                met = model.metabolites.get_by_id(met_id)
            except KeyError:
                met = Metabolite(met_id)
                model.add_metabolites(met)
                
                #reaction.add_metabolites({met: num})
    #model.add_reaction(reaction)
    return model  

def get_reaction_uni_met(equation):
    uni_str_line=[]
    if re.search('-->',equation):
        str_line=equation.replace('_e','').replace('_p','').replace('_c','').split(' --> ')
        str_line = [i for i in str_line if i != '']
        uni_str_line=list(set(str_line))
    elif re.search('<--',equation):
        str_line=equation.replace('_e','').replace('_p','').replace('_c','').split(' <-- ')
        str_line = [i for i in str_line if i != '']
        uni_str_line=list(set(str_line))
    elif re.search('<=>',equation):
        str_line=equation.replace('_e','').replace('_p','').replace('_c','').split(' <=> ')
        str_line = [i for i in str_line if i != '']
        uni_str_line=list(set(str_line))

    return uni_str_line


def exchanges(model):
    """Exchange reactions in model.
    Reactions that exchange mass with the exterior. Uses annotations
    and heuristics to exclude non-exchanges such as sink reactions.
    """
    exchange_list=[]
    try: 
        find_boundary_types(model,'exchange')
    except:
        pass
    else:
        for reachr in find_boundary_types(model,'exchange'):
            exchange_list.append(reachr.id)
    return exchange_list


def create_model_change_json(model):
    # 生成model_change
    model_change_json={}
    model_change_json['reaction_infor']={}
    model_change_json['gene_infor']={}
    model_change_json['exchange_infor']={}  #add_reaction后 model.exchanges没了
    model_change_json['bigg_metabolite_ids']={}  
    model_change_json['add_reaction_infor']={}    # 与modify生成的model_change_json_final.keys()保持一致
    model_change_json['del_reaction_infor']={}
    model_change_json['del_gene_infor']={}
    for eachr in model.reactions:
        model_change_json['reaction_infor'][eachr.id]={}
        model_change_json['reaction_infor'][eachr.id]['reaction_name']=eachr.name
        model_change_json['reaction_infor'][eachr.id]['equation']=eachr.reaction
        model_change_json['reaction_infor'][eachr.id]['equation by name']=eachr.build_reaction_string(True).replace('O2 O2','O2').replace('CO2 CO2','CO2').replace('H2O H2O','H2O')
        model_change_json['reaction_infor'][eachr.id]['gene_reaction_rule']=eachr.gene_reaction_rule
        model_change_json['reaction_infor'][eachr.id]['lower_bound']=model.reactions.get_by_id(eachr.id).lower_bound
        model_change_json['reaction_infor'][eachr.id]['upper_bound']=model.reactions.get_by_id(eachr.id).upper_bound
    for eachg in model.genes:
        if eachg.name:
            model_change_json['gene_infor'][eachg.id]={}
            model_change_json['gene_infor'][eachg.id]['gene_name']=eachg.name
            reaction_id_list=[]
            reaction_name_list=[]
            reaction_equ_list=[]
            for eachr in model.genes.get_by_id(eachg.id).reactions:
                reaction_name_list.append(eachr.name)
                reaction_equ_list.append(eachr.reaction)
            model_change_json['gene_infor'][eachg.id]['reaction_name_list']='; '.join(reaction_name_list) 
            model_change_json['gene_infor'][eachg.id]['reaction_equ_list']='; '.join(reaction_equ_list)    
    if len(exchanges(model))>0:
        for eachr in exchanges(model):
            model_change_json['exchange_infor'][eachr]={}
            model_change_json['exchange_infor'][eachr]['reaction_name']=model.reactions.get_by_id(eachr).name
            model_change_json['exchange_infor'][eachr]['lower_bound']=model.reactions.get_by_id(eachr).lower_bound
            model_change_json['exchange_infor'][eachr]['upper_bound']=model.reactions.get_by_id(eachr).upper_bound
    else:
        model_change_json['exchange_infor']={}
    for eachm in model.metabolites:
        model_change_json['bigg_metabolite_ids'][eachm.id]={}
        model_change_json['bigg_metabolite_ids'][eachm.id]['bigg_id']=eachm.id
        model_change_json['bigg_metabolite_ids'][eachm.id]['name']=eachm.name
    return model_change_json


def build_model_from_FBA_result_plus(model_pfba_solution,rxn_switch,coordinates, inimodel, model_change_json, substrate_id,product_id):
    """由FBA结果构建模型
    Arguments
    ----------
    * model_pfba_solution: 模型计算结果.
    :return: cobra model.
    """
    #print(substrate_id,product_id)
    metlink_model = Model('metlink_model')
    #print(substrate_id,product_id)
    #Demand反应 
    if product_id.startswith('DM_'):
        product_id=product_id.replace('DM_','')
        product_add = inimodel.metabolites.get_by_id(product_id)
        demand =inimodel.add_boundary(product_add, type='demand')
    product_id =product_id.replace('__D_e','_e').replace('__L_e','_e')    
        
    #Sink反应
    if substrate_id.startswith('SK_'):
        substrate_id=substrate_id.replace('SK_','')
        substract_add=inimodel.metabolites.get_by_id(substrate_id)
        sink_reaction=inimodel.add_boundary(substract_add, type="sink")
        #substrate_id=sink_reaction.id
    elif re.search('_reverse',substrate_id):
        substrate_id =substrate_id.replace('_reverse','')
    substrate_id =substrate_id.replace('__D_e','_e').replace('__L_e','_e')    
    #    
    substrate_sub=substrate_id.replace('_e','').replace('_p','').replace('_c','')
    product_sub=product_id.replace('_e','').replace('_p','').replace('_c','')
    #print(substrate_id,product_id)
    for index, row in model_pfba_solution.iterrows():
        new_equ=row['equ'].replace('__D_e','_e').replace('__L_e','_e')
        uni_str_line=get_reaction_uni_met(new_equ)    
        if len(uni_str_line)==1 and substrate_sub== uni_str_line[0]:  
            pass
        elif len(uni_str_line)==1 and product_sub== uni_str_line[0]: 
            pass         
        else:
            metlink_model=get_model_met_from_reaction_string(metlink_model,new_equ,fwd_arrow='-->', rev_arrow='<--', reversible_arrow='<=>', term_split='+')
            reaction = Reaction(index) 
            try:
                reaction_infor=inimodel.reactions.get_by_id(index)
            except:
                reaction.name=index
                reaction_ano='Reaction id: <br/>'+str(index)+'<hr/>'+'Reaction name: <br/>'+str(index)+'<hr/>'+\
                        'Reaction subsystem: <br/>'+'none'+'<hr/>'+'Reaction equation: <br/>'+\
                        str(new_equ)+'<hr/>'+'Reaction gene_reaction_rule: <br/>'+'none'+'<hr/>'
            else:
                reaction.name=reaction_infor.name
                reaction_ano='Reaction id: <br/>'+str(index)+'<hr/>'+'Reaction name: <br/>'+str(reaction_infor.name)+'<hr/>'+\
                            'Reaction subsystem: <br/>'+str(reaction_infor.subsystem)+'<hr/>'+'Reaction equation: <br/>'+\
                            str(reaction_infor.reaction)+'<hr/>'+'Reaction gene_reaction_rule: <br/>'+reaction_infor.gene_reaction_rule+'<hr/>'
            reaction.notes['map_info'] = {}
            reaction.notes['map_info']['annotation']=reaction_ano
            if rxn_switch == 'T':
                #reaction.notes['map_info']['display_name'] = reaction.id + " "+ str("%.2f" %(abs(row['fluxes'])))
                reaction.notes['map_info']['display_name'] = str("%.2f" %(abs(row['fluxes'])))
            elif rxn_switch == 'F':
                #reaction.notes['map_info']['display_name'] = reaction.name + " "+str("%.2f"%(abs(row['fluxes'])))
                reaction.notes['map_info']['display_name'] = str("%.2f"%(abs(row['fluxes'])))
            #d.notes.map_info.group   
            #[undefined, 'ko', 1, 2, 3, 4, 5, 6, 7, 8]
            if coordinates:#坐标轴赋值
                if 'reactions' in coordinates.keys():
                    if index in coordinates['reactions']:
                        reaction.notes['map_info']['x'] = float(coordinates['reactions'][index]['x'])
                        reaction.notes['map_info']['y'] = float(coordinates['reactions'][index]['y'])

            # if index in model_change_json['exchange_infor'].keys():
            #     reaction.notes['map_info']['group']=4
            if index in model_change_json['add_reaction_infor'].keys():
                reaction.notes['map_info']['group']=5       
            elif index in model_change_json['del_reaction_infor'].keys():
                reaction.notes['map_info']['group']=1    
            elif index in model_change_json['del_gene_infor'].keys():
                reaction.notes['map_info']['group']=1 
            else:
                reaction.notes['map_info']['group']=3
            metlink_model.add_reactions([reaction])
            reaction.build_reaction_from_string(new_equ) 

    #去除H
    H_list=['h_c','h_p','h_e','C00080','C00080[c]', 'C00080[p]','C00080[e]','PROTON', 'MNXM1', 'cpd00067', 'cpd00067_c0','cpd00067_p0','cpd00067_e0', 'HMDB59597', 'GPRLSGONYQIRFK-UHFFFAOYSA-N']
    for each_h in H_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_h)
        except:
            pass
        else:
            hc_id.remove_from_model()

    #去除H2O
    H2O_list=['h2o_c','h2o_p','h2o_e','C00001', 'C01328','C00001[c]','C00001[p]','C00001[e]', 'C01328[c]','C01328[p]','C01328[e]', 'CPD-15815', 'HYDROXYL-GROUP', 'WATER', 'MNXM2', 'cpd00001', 'cpd00001_c0', 'cpd00001_p0','cpd00001_e0','cpd15275', 'cpd27222', 'HMDB01039', 'HMDB02111', 'XLYOFNOQVPJJNP-UHFFFAOYSA-N']
    for each_h2o in H2O_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_h2o)
        except:
            pass
        else:
            hc_id.remove_from_model()
        #去除Pi
    Pi_list=['pi_c','pi_p','pi_e','C00009', 'C00009[c]','C00009[p]','C00009[e]','CPD-16459', 'PHOSPHATE-GROUP', 'Pi', 'MNXM9', 'cpd00009', 'cpd00009_c0', 'cpd00009_p0','cpd00009_e0','cpd27787', 'HMDB00973', 'HMDB01429', 'HMDB02105', 'HMDB02142', 'HMDB05947', 'NBIIXXVUZAFLBC-UHF']
    for each_pi in Pi_list:
        try:
            hc_id=metlink_model.metabolites.get_by_id(each_pi)
        except:
            pass
        else:
            hc_id.remove_from_model()

    #print(len(metlink_model.reactions.get_by_id('DM_pyr_c').products))  
    display_name_format = (lambda met: re.sub('__[D,L]', '', met.id[:-2].upper()))
    model_met_list=[]
    for eachm in metlink_model.metabolites:
        model_met_list.append(eachm.id)
    print(substrate_id,product_id)
    if substrate_id.endswith('_e'):
        tmp_p="_e".join(substrate_id.split('_e')[:-1])+'_p'
        tmp_c="_e".join(substrate_id.split('_e')[:-1])+'_c'
        if substrate_id=='glc__D_e' or substrate_id=='glc_D_e':
            if tmp_c in model_met_list:
                substrate_id = tmp_c    
            elif tmp_p in model_met_list:
                substrate_id = tmp_p   
        else:
            if tmp_c in model_met_list:
                substrate_id = tmp_c    
            elif tmp_p in model_met_list:
                substrate_id = tmp_p 
 
    if product_id.endswith('_e'):
        tmp_p="_e".join(product_id.split('_e')[:-1])+'_p'
        tmp_c="_e".join(product_id.split('_e')[:-1])+'_c'
        if tmp_c in model_met_list:
            product_id = tmp_c    
        elif tmp_p in model_met_list:
            product_id = tmp_p 
  
    #print(substrate_id,product_id)
    for eachmet in metlink_model.metabolites:
        if re.search('_c',eachmet.id):
            eachmet.compartment='c'
        elif re.search('_p',eachmet.id):
            eachmet.compartment='p'
        elif re.search('_e',eachmet.id):
            eachmet.compartment='e'  
            
        try:
            inimodel.metabolites.get_by_id(eachmet.id)
        except:
            eachmet_id= "_e".join(eachmet.id.split('_e')[:-1])+'__D_e' 
            try:
                inimodel.metabolites.get_by_id(eachmet_id)
            except:
                eachmet_id= "_e".join(eachmet.id.split('_e')[:-1])+'__L_e' 
                try:
                    inimodel.metabolites.get_by_id(eachmet_id)
                except:
                    pass
                else:
                    metabolite_infor=inimodel.metabolites.get_by_id(eachmet_id)
            else:
                metabolite_infor=inimodel.metabolites.get_by_id(eachmet_id)
        else:
            metabolite_infor=inimodel.metabolites.get_by_id(eachmet.id)
        eachmet.name='Metabolite id: <br/>'+str(eachmet.id)+'<hr/>'+'Metabolite name: <br/>'+str(metabolite_infor.name)+'<hr/>'+\
                        'Metabolite formula: <br/>'+str(metabolite_infor.formula)+'<hr/>'
        eachmet.notes['map_info'] = {}
        #eachmet.notes['map_info']['annotation']=metabolite_ano
        new_list=['nadh_c','nad_c','atp_c','atp_p', 'atp_e','C00002', 'C00002[c]','C00002[p]','C00002[e]','ATP', 'MNXM3', 'cpd00002', 'cpd00002_c0','cpd00002_p0','cpd00002_e0', 'HMDB00538', 'ZKHQWZAMYRWXGA-KQYNXXCUSA-J','adp_c','adp_p','adp_e','C00008','C00008[c]', 'ADP', 'MNXM7','cpd00008','cpd00008_c0','cpd00008_p0','cpd00008_e0', 'HMDB01341', 'XTWYTFMLZFPYCI-KQYNXXCUSA-K','nadph_c','NADPH','C00005','C00005[c]','MNXM6','HMDB00221', 'HMDB00799', 'HMDB06341','cpd00005','ACFIXJIJDZMPPO-NNYOXOHSSA-J','nadph','META:NADPH','cpd00005_c0','nadp_c','NADP','MNXM5','C00006','C00006[c]','HMDB00217','cpd00006','XJLXINKUBYWONI-NNYOXOHSSA-K','META:NADP','nadp','cpd00006_c0']
        if eachmet.id in new_list:    # [ATP,'cpd00002',ADP,'cpd00008',鲜肉(鲑鱼)色]
            eachmet.notes['map_info']['color'] = "#FA8072"
        # elif met.id[:-3] in ['cpd00003', 'cpd00004', 'cpd00005', 'cpd00006']:# [nad、nadh、nadp、nadph]绿宝石
        #     met.notes['map_info']['color'] = "#008000"
        else:   #  其他为适中的碧绿色
            eachmet.notes['map_info']['color'] = "#00FA9A" 
        #print(eachmet.id,substrate_id) 
#         if eachmet.id=='glc__D_e':
#              #print(substrate_id)
#             eachmet.notes['map_info']['substrate'] ='glc__D_e' 
        if eachmet.id==substrate_id:
            #print(substrate_id)
            eachmet.notes['map_info']['substrate'] =substrate_id          
            #print(eachmet.id,substrate_id)    
        if eachmet.id==product_id:
            eachmet.notes['map_info']['product'] =product_id
        if coordinates:#坐标轴赋值
            cormet=display_name_format(eachmet)
            if 'metabolites' in coordinates.keys():
                #print(cormet,cormet.lower())
                if cormet.lower() in coordinates['metabolites'].keys():
                    eachmet.notes['map_info']['x'] = float(coordinates['metabolites'][cormet.lower()]['x'])
                    eachmet.notes['map_info']['y'] = float(coordinates['metabolites'][cormet.lower()]['y'])
        if eachmet.id.endswith('_c') or eachmet.id.endswith('_C') or eachmet.id.endswith('_p') or eachmet.id.endswith('_P') or eachmet.id.endswith('_e') or eachmet.id.endswith('_E'):
            pass
        else:
            eachmet.id=eachmet.id+'_c'
    return metlink_model


def get_model_information(model):
    All_met = []
    All_rxn = []
    for met in model.metabolites:
        All_met.append(met.id)
        # print(met.id,met.name,'............................')
    for rxn in model.reactions:
        All_rxn.append(rxn.id)  # 所有 rxn
    uni_All_met=list(set(All_met))
    uni_All_rxn=list(set(All_rxn))
    return uni_All_met, uni_All_rxn

def pathway_tsv_visualization(pathway_tsv_file,substrate_id,product_id,cofactor_switch,flow,outputdir, rxnId, model):
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
                            'didp_e', 'mn2_e', 'dcdp_e', 'nh4_e', 'amp_e', 'fe3_e', 'nadp_e', 'so2_e', 'h2o2_e', 'mg2_e', 'flxso_c','flxr_c','WATER_c',\
                            'PROTON_c','NADH_c','NAD_c','NADPH_c','NADP_c','ATP_c','ADP_c','Pi_c','pqq_c','pqqh2_c']
    #由通量数据构建model 便于后续获得metlink
    rxn_switch = 'F'
    cofactor_switch = "F"
    id_or_name = True
    # flow = True 
    coordinates= False

    if flow == "True":
        flow = True
    else:
        flow = False
    print('*'*20)
    print('xxxxxxxxxxxxx')

    # outfolder = pathway_outputpath # 读取途径生成文件
    # for file_name in os.listdir(outfolder): #得到文件夹下所有文件名
        # if 'summary' not in file_name and 'tsv' in file_name: 
    pathway_df = pd.read_table(pathway_tsv_file,header=None)
    pathway_df.columns = ['id','fluxes','equ','equName','bounds','check']
    pathway_df.set_index('id',inplace=True)

    model_change_json = create_model_change_json(model)

    metlink_model = build_model_from_FBA_result_plus(pathway_df,rxn_switch,coordinates, model, model_change_json, substrate_id,product_id)
    
    # 可视化
    All_met, All_rxn= get_model_information(metlink_model)
    #print(All_rxn)
    common_factor_list=list(set(total_common_factor_list).intersection(set(All_met))) #取交集
    print('common_factor_list:',common_factor_list)
    if cofactor_switch == "T":
        d3flux.update_cofactors(metlink_model, common_factor_list)

    

    html = d3flux.flux_map(metlink_model, overwrite_reversibility=True,
                        excluded_metabolites = common_factor_list,excluded_reactions =exclude_rxn,
                        figsize=(1500, 1500), display_name_format=id_or_name,
                        hide_unused="true", hide_unused_cofactors="true",
                        flux_dict=pathway_df['fluxes'],
                        flowLayout = flow)  # 是否要隐藏flux为0的

    # show_factor_file = './show_factor.html'  # 前端d3flux可视化文件
    html_file = os.path.join(outputdir,f"{rxnId}.html")
    with open(html_file, "w")as wf:
        wf.write(html.data)
    return html_file



