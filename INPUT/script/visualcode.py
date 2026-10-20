# coding:utf-8
"""
Author: Jingyi Cai (2024-2026)
Function: Identify the biomass reaction and build the metabolite-link model used for pathway figures.
Input: Called by target.py with a COBRA model and FBA result tables. Uses the identifier lists in defaults.py.
Output: Returns biomass reaction information and map tables to the caller.
"""
import json
import pandas as pd
import cobra
import itertools
import re
import os
from cobra import Model, Reaction, Metabolite
from cobra.util.solver import linear_reaction_coefficients
from defaults import *
from cobra.medium import find_boundary_types
from task_config import TEMP_DIR

def download_reference_from_s3(s3,bucket,obj):
    """
        download reference data /tmp/
    """
    object_name = obj.split('/')[-1]
    local_reference = os.path.join(TEMP_DIR, object_name)
    s3.Object(bucket, obj).download_file(local_reference)
    return local_reference

def download_file_from_s3(s3,bucket,obj,workdir):
    """
        download file data workdir
    """
    object_name = obj.split('/')[-1]
    local_reference = os.path.join(workdir,object_name)
    s3.Object(bucket, obj).download_file(local_reference)
    return local_reference

def change_file_acl(s3,bucket,object_key):
    """
    修改文件的acl，便于客户下载
    """
    try:
        s3.ObjectAcl(bucket,object_key).put(ACL="public-read")
        url = "https://"+bucket+".s3.amazonaws.com/%s" % object_key
        return url
    except Exception as e:
        print(str(e))
        print('%s is not existed in s3' %object_key )
        return False

def json_load(path):
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary

def get_fluxes_detail_in_model(model, call_method, fluxes_outfile):
    """Get the detailed information of each reaction

    Arguments
    ----------
    * model: cobra.Model.
    * fluxes_outfile: reaction flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.

    :return: fluxes, kcat, MW and kcat_MW in dataframe.
    """
    try:
        if call_method=='pFBA':
            model_pfba_solution = cobra.flux_analysis.pfba(model)
        elif call_method=='FBA':   
            model_pfba_solution = model.optimize()

    except:
        print('Can not solve the model!')
    else:
        model_pfba_solution = model_pfba_solution.to_frame()
        model_pfba_solution_detail = {}
        for index, row in model_pfba_solution.iterrows():
            reaction_detail = model.reactions.get_by_id(index)
            model_pfba_solution_detail[index]={}
            if abs(row['fluxes'])>1e-9:
                useflux=round(row['fluxes'],6)
            else:
                useflux=0
                
            model_pfba_solution_detail[index]['fluxes'] = useflux    
            #for eachkey in reaction_detail.annotation.keys():
            #    model_pfba_solution_detail[index][eachkey] = reaction_detail.annotation[eachkey]   
            if 'ec-code' in reaction_detail.annotation.keys():
                model_pfba_solution_detail[index]['ec-code'] = reaction_detail.annotation['ec-code']  
            else:
                model_pfba_solution_detail[index]['ec-code'] = 'none'
            if reaction_detail.gene_reaction_rule:
                model_pfba_solution_detail[index]['gene id'] =reaction_detail.gene_reaction_rule
            else:
                model_pfba_solution_detail[index]['gene id'] = 'none'
            model_pfba_solution_detail[index]['equ'] = reaction_detail.reaction
            
        model_pfba_solution_detail_df = pd.DataFrame(model_pfba_solution_detail).T
        model_pfba_solution_detail_df=model_pfba_solution_detail_df[abs(model_pfba_solution_detail_df['fluxes'])>0]
        model_pfba_solution_detail_df.to_csv(fluxes_outfile,sep='\t')
        return model_pfba_solution_detail_df
    
def json_write(path, dictionary):
    """Writes a JSON file at the given path with the given dictionary as content.

    Arguments
    ----------
    * path:   The path of the JSON file that shall be written
    * dictionary: The dictionary which shalll be the content of
      the created JSON file
    """
    json_output = json.dumps(dictionary, indent=4)
    with open(path, "w", encoding="utf-8") as f:
        f.write(json_output)


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
    for substr, factor in ((reactant_str, -1), (product_str, 1)):
        if len(substr) == 0:
            continue
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

def build_model_from_FBA_result_plus(inimodel,model_pfba_solution,rxn_switch,model_change_json,coordinates,substrate_id,product_id):
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
        
def get_model_information(model):#Mr.Mao改
    All_rxn = []
    All_met = []
    for met in model.metabolites:
        All_met.append(met.id)
    for rxn in model.reactions:
        All_rxn.append(rxn.id)  # 所有 rxn
    uni_All_met=list(set(All_met))
    uni_All_rxn=list(set(All_rxn))
    return uni_All_met, uni_All_rxn

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


def getfluxdata(norm_model,product_id, fluxA, fluxB,fluxC,model_name="e_coli_core"):
    
    flux_solution = cobra.flux_analysis.pfba(norm_model)  
    flux_solution_df=flux_solution.to_frame()
    Optimal_rate=flux_solution_df.loc[product_id,'fluxes']
    flux_solution_df_select=flux_solution_df[abs(flux_solution_df['fluxes'])>0.0000001]#只储存大于0的值e-6

    # for fluxA
    model_pfba_solutionA=pd.DataFrame() 
    if flux_solution_df_select.shape[0] > 0 and Optimal_rate>0:
        for index, row in flux_solution_df_select.iterrows():
                    model_pfba_solution.loc[index,'reaction_id']=index
                    model_pfba_solution.loc[index,'reaction_name']=norm_model.reactions.get_by_id(index).name
                    # model_pfba_solution.loc[index,'fluxes']=abs(round(row['fluxes'],3))
                    model_pfba_solution.loc[index,'fluxes']=round(row['fluxes'],7)
                    model_pfba_solution.loc[index,'gene id']=norm_model.reactions.get_by_id(index).gene_reaction_rule
                    #model_pfba_solution.loc[index,'ec-code']=norm_model.reactions.get_by_id(index).annotation['ec-code'] 
                    model_pfba_solution.loc[index,'equ']=norm_model.reactions.get_by_id(index).reaction
                    model_pfba_solution.loc[index,'equ_name']=norm_model.reactions.get_by_id(index).build_reaction_string(True).replace('O2 O2','O2').replace('CO2 CO2','CO2').replace('H2O H2O','H2O')
                    model_pfba_solution.loc[index,'abs_fluxes']=abs(round(row['fluxes'],7))

    flux_outfile='./%s-d3flux.tsv'%model_name
    model_pfba_solution=model_pfba_solution.sort_values(by=['abs_fluxes'],ascending=False)
    model_pfba_solution.to_csv(flux_outfile,index=False,sep='\t')
    return model_pfba_solution
def colse_model_default_carbon(model,met_C):
    #关闭碳源
    default_C_list=[]
    exclude_C_reaction_list=['EX_BIOTIN_reverse','EX_BIOTIN','EX_biotin_e', 'EX_btn', 'EX_btn_e', 'EX_btn_e_', 'EX_btn(e)', 'EX_btn[e]', 'EX_btn_LPAREN_e_RPAREN_','EX_biotin_e_reverse', 'EX_btn_reverse', 'EX_btn_e_reverse', 'EX_btn_e__reverse', 'EX_btn(e)_reverse', 'EX_btn[e]_reverse', 'EX_btn_LPAREN_e_RPAREN__reverse','META:TRANS-RXN0-240','MNXR96333']
    for eachr in medium(model):
        #print(eachr)
        #r_bounds=norm_model.reactions.get_by_id(eachr).bounds
        product_equ=model.reactions.get_by_id(eachr)
        if len(product_equ.reactants)>1:
            pass
        else:
            for eachm in product_equ.reactants:
                #print(eachr,eachm.id,eachm.formula)
                #met_coef=abs(product_equ.get_coefficient(eachm.id))
                try:
                    met_formula=get_composition_from_formula(eachm.formula)
                except:
                    try:
                        met_formula=get_composition_from_formula(met_C.loc[eachm.id,'formula'])
                        #print(eachr,eachm.id,met_formula)
                    except:
                        print('No formula!')
                    else:
                        if 'C' in met_formula.keys():#是碳源才关闭
#                             if eachr not in exclude_C_reaction_list:
#                                 model.reactions.get_by_id(eachr).bounds = (0, 1000) 
#                             default_C_list.append(eachr)
                            if model.reactions.get_by_id(eachr).lower_bound<0:
                                if eachr not in exclude_C_reaction_list:
                                    model.reactions.get_by_id(eachr).bounds = (0, 1000) 
                                default_C_list.append(eachr)
                            else:
                                if eachr not in exclude_C_reaction_list:
                                    model.reactions.get_by_id(eachr).bounds = (0, 0) 
                                default_C_list.append(eachr)  
                else:
                    if 'C' in met_formula.keys():#是碳源才关闭
#                         if eachr not in exclude_C_reaction_list:
#                             model.reactions.get_by_id(eachr).bounds = (0, 1000) 
#                         default_C_list.append(eachr)
                        if model.reactions.get_by_id(eachr).lower_bound<0:
                            if eachr not in exclude_C_reaction_list:
                                model.reactions.get_by_id(eachr).bounds = (0, 1000) 
                            default_C_list.append(eachr)
                        else:
                            if eachr not in exclude_C_reaction_list:
                                model.reactions.get_by_id(eachr).bounds = (0, 0) 
                            default_C_list.append(eachr)
        if len(product_equ.products)>1:
            pass
        else:
            for eachm in product_equ.products:
                #print(eachr,eachm.id,eachm.formula)
                #met_coef=abs(product_equ.get_coefficient(eachm.id))
                try:
                    met_formula=get_composition_from_formula(eachm.formula)
                except:
                    try:
                        met_formula=get_composition_from_formula(met_C.loc[eachm.id,'formula'])
                        #print(eachr,eachm.id,met_formula)
                    except:
                        print('No formula!')
                    else:
                        if 'C' in met_formula.keys():#是碳源才关闭
                            if model.reactions.get_by_id(eachr).lower_bound<0:
                                if eachr not in exclude_C_reaction_list:
                                    model.reactions.get_by_id(eachr).bounds = (0, 1000) 
                                default_C_list.append(eachr)
                            else:
                                if eachr not in exclude_C_reaction_list:
                                    model.reactions.get_by_id(eachr).bounds = (0, 0) 
                                default_C_list.append(eachr)                                
                else:
                    if 'C' in met_formula.keys():#是碳源才关闭
                        if model.reactions.get_by_id(eachr).lower_bound<0:
                            if eachr not in exclude_C_reaction_list:
                                model.reactions.get_by_id(eachr).bounds = (0, 1000) 
                            default_C_list.append(eachr)
                        else:
                            if eachr not in exclude_C_reaction_list:
                                model.reactions.get_by_id(eachr).bounds = (0, 0) 
                            default_C_list.append(eachr)                              
                    
    return [model,default_C_list]

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

def deal_dict(model_input,inputdic):
    # 加载任务字典


    # 提取产物反应名称
    productDM_name = inputdic['product'].replace("user_add_", "")
    productEX_name = None  # 初始化变量

    if 'user_add_DM_' in inputdic['product']:
        # 替换反应名称，并调整为外环境（`_c` 改为 `_e`）
        productEX_name = inputdic['product'].replace('user_add_DM_', 'EX_')
        productEX_name = productEX_name[::-1].replace('c', 'e', 1)[::-1]

    # 检查是否存在 EX_产品
    product_used = None  # 初始化使用的产品为 None
    if productEX_name and productEX_name in model_input.reactions:
        product_used = productEX_name
        model_input.reactions.get_by_id(productDM_name).bounds = (0, 0)
        model_input.reactions.get_by_id(product_used).bounds = (0, 1000)
    elif productDM_name in model_input.reactions:
        product_used = productDM_name

    # 更新字典中的产品信息
    inputdic['product'] = product_used

    return inputdic, model_input, product_used, productEX_name, productDM_name



def Error_Code(inputdic,model_input,productEX_name,productDM_name):
    class ModelError(Exception):
        def __init__(self, error_code, message):
            self.error_code = error_code
            self.message = message
            super().__init__(self.message)

        # 错误代码001：biomass不在表里
    class BiomassNotFoundError(ModelError):
        def __init__(self):
            error_code = '001'
            message = "biomass 不在表里"
            super().__init__(error_code, message)

        # 错误代码002：底物不在表里
    class SubstrateNotFoundError(ModelError):
        def __init__(self, substrate_name):
            error_code = '002'
            message = f"底物 '{substrate_name}' 不在表里"
            super().__init__(error_code, message)

    # 错误代码003：EX_产品、DM_产品都不在表里
    # 错误代码003：EX_产品、DM_产品都不在表里
    class ProductNotFoundError(ModelError):
        def __init__(self, productEX_name, productDM_name):
            error_code = '003'
            message = ""
            product_used = None  # 初始化使用的产品为 None
            if productEX_name in model_input.reactions:
                # 使用EX_产品进行后续操作
                product_used = productEX_name
                error_code = None  # 设置错误代码为 None
            elif productDM_name in model_input.reactions:
                # 使用DM_产品进行后续操作
                product_used = productDM_name
                error_code = None  # 设置错误代码为 None
            else:
                if productEX_name and productDM_name:
                    message = f"产品 '{productEX_name}' 和 '{productDM_name}' 都不在表里"
                elif productEX_name:
                    message = f"产品 '{productEX_name}' 不在表里"
                elif productDM_name:
                    message = f"产品 '{productDM_name}' 不在表里"
            super().__init__(error_code, message)
            self.product_used = product_used 
    



    # 错误代码004：biomass为0，不生长
    class GrowthError(ModelError):
        def __init__(self):
            error_code = '004'
            message = "不生长"
            super().__init__(error_code, message)


    # 错误代码005：product为0，无产品
    class ProductError(ModelError):
        def __init__(self):
            error_code = '005'
            message = "无产品"
            super().__init__(error_code, message)

        
    try:

        # 检查是否存在 EX_产品
        # productEX_name = inputdic['product']

        product_used = None  # 初始化使用的产品为 None

        if productEX_name in model_input.reactions:
            product_used = productEX_name
        elif productDM_name in model_input.reactions:
            product_used = productDM_name
        else:
            raise ProductNotFoundError(productEX_name, productDM_name)
        
        # 检查底物是否存在
        substrate_name = inputdic['substrate']
        if substrate_name not in model_input.reactions:
            raise SubstrateNotFoundError(substrate_name)
        
        # 检查biomass是否存在
        biomass_name = inputdic['biomass']
        if biomass_name not in model_input.reactions:
            raise BiomassNotFoundError()
        

        # 检查是否生长
        model_input.objective = inputdic['biomass']
        objvalue2 = model_input.optimize().objective_value
        
        if objvalue2 is None or objvalue2 <=1e-5:
            raise GrowthError()
        
        # 检查是否无产品
        model_input.reactions.get_by_id(product_used).bounds = (0, 1000)
        model_input.objective = product_used
        objvalue1 = model_input.optimize().objective_value

        if objvalue1 is None or objvalue1 <=1e-5:
            raise ProductError()
        

    except ModelError as e:
        if e.error_code is not None:   # 捕获自定义异常并处理
            print("错误代码:", e.error_code)
            print("错误描述:", e.message)
            error_code=e.error_code
    else:
        # 没有发生异常的情况下执行的代码
        error_code = '100'
    return error_code


def get_initial_obj_info( model):
    """"""
    biomassIds = []
    try:
        initial_rxn_id = list(linear_reaction_coefficients(model).keys())[0].id
    except:
        initial_rxn_id = ''
    first_biomass_rxn_id = first_get_biomass_rxn(model, initial_rxn_id)
    if first_biomass_rxn_id:
        biomass_rxn_id = first_biomass_rxn_id
    else:
        all_biomassIds = get_biomass_rxn(model)
        if len(all_biomassIds) > 1: # 如果找到了多个biomass方程，需要进行判断
            with model as model:
                for rxnId in all_biomassIds:
                    model.objective = rxnId
                    if model.slim_optimize() > 1e-6:
                        biomassIds.append(rxnId)
        else:
            biomassIds = all_biomassIds
        # biomassIds=['BIOMASS_glyc']
        if biomassIds:
            if initial_rxn_id in biomassIds:
                biomass_rxn_id = initial_rxn_id
            else:
                biomass_rxn_id = biomassIds[0]
        else:
            biomass_rxn_id = initial_rxn_id  #以前是在最后仍然找不到biomass的时候把目标反应作为biomass，这样可能会使一些非biomass反应作为biomass，所以现在改为直接告诉用户找不到
            # biomass_rxn_id = ''
    return biomass_rxn_id

def first_get_biomass_rxn(model, initial_rxn_id):
    """"""
    if initial_rxn_id:
        initial_rxn = model.reactions.get_by_id(initial_rxn_id)
        mets = list(initial_rxn.metabolites.keys())[0]
        if len(initial_rxn.metabolites) == 1 and 'biomass' in mets.id.lower() or mets.id.lower() in 'biomass' or mets.name in BIOMASS or 'biomass' in mets.name.lower():
            ids = mets.id
            for rxn in model.reactions:
                if any(ids in k.id for k in rxn.products if len(rxn.products) > 1):
                    return rxn.id
    for rxn in model.reactions:
        if len(rxn.metabolites) == 1:
            mets = list(rxn.metabolites.keys())[0]
            if 'biomass' in mets.id.lower() or mets.id.lower() in 'biomass' or mets.name in BIOMASS or 'biomass' in mets.name.lower():
                ids = mets.id
                for rxn in model.reactions:
                    if any(ids in k.id for k in rxn.products if len(rxn.products) > 1):
                        return rxn.id
    return ""

def get_biomass_rxn(model):
    """"""
    biomassId, num = [], 0
    for rxn in model.reactions:
        all_met_id = [met.id for met in rxn.metabolites]
        all_met_name = [met.name.lower() for met in rxn.metabolites]
        if len(set(PROTEIN) & set(all_met_id)) != 0 and len(set(DNA_RNA) & set(all_met_id)) >= 2:
            biomassId.append(rxn.id)
        if len(set(PROTEIN_COMPOSITION + DNA_COMPOSITION + RNA_COMPOSITION) & set(all_met_id)) >= 16 and len(set(DNA_COMPOSITION + RNA_COMPOSITION) & set(all_met_id)) >= 4:
            biomassId.append(rxn.id)
        if any(['protein' in k for k in all_met_name]) : num += 1
        if any(['dna' in k for k in all_met_name]) : num += 1
        if any(['rna' in k for k in all_met_name]) : num += 1
        if any(['lipid' in k for k in all_met_name]) : num += 1
        if num == 4:
            biomassId.append(rxn.id)
        num = 0
    biomassId = list(set(biomassId))
    return biomassId

def optmeprocess(norm_model,substrate_id,product_id,model_substrate,model_product,met_C,biomass):
    
    substrate_name = model_substrate.loc[substrate_id,'substrate_name']
    substrate_uptake_rate = 10


    input_json={
        "model": "wangry@tib.cas.cn_20241203-092527_iML1515.json",
        "substrate": "EX_glc__D_e",
        "biomass": "CG_biomass_cgl_ATCC13032",
        "substrate_uptake_rate": 10,
        "ATPM": 0,
        "product": "EX_ac_e",
        "oxygenstate": "aerobic",
        "email": "wangry@tib.cas.cn",
        "taskname": ["FSEOF","loopless_optforce_MUST","iBridge"],
        "excluded_rxns": [

        ],
        "CO2_flag": "True",
        "min_growth": 10,
        "O2": "EX_o2_e",
        "species": "Escherichia_coli_str._K-12_substr._MG1655",
        "ID": "wangry@tib.cas.cn_20241203-092527_iML1515"
    }

    
    input_json["substrate"]=substrate_id
    input_json["product"]=product_id
    input_json["substrate_uptake_rate"]=substrate_uptake_rate
    input_json["biomass"]=biomass

    
    substrate_name = model_substrate.loc[substrate_id,'substrate_name']
    if product_id:
        product_name = model_product.loc[product_id,'obj_name']
        product_m_id=model_product.loc[product_id,'obj_id']
        product_type=model_product.loc[product_id,'type']

    substrate_m_id=model_substrate.loc[substrate_id,'substrate_id']

    product_add = norm_model.metabolites.get_by_id(product_m_id)
    if product_id:
        if re.search('user_add_DM_',product_id):
            if 'DM_'+product_m_id not in norm_model.reactions:
                try:
                    product_add = norm_model.metabolites.get_by_id(product_m_id)
                    demand =norm_model.add_boundary(product_add, type='demand') #add demand reaction as the objective

                except:
                    print('The model cannot solve for this product!')
                    raise ValueError("The model cannot solve for this product!")
                else:
                    norm_model.objective = demand
                    product_id=demand.id
            else:
                norm_model.objective ='DM_'+product_m_id
        else:
            norm_model.objective = product_id

    ##关闭培养基其他底物
    if substrate_id !='substrate_use_in_model':
        try:
            find_boundary_types(norm_model,'exchange')
        except:
            print('No medium information!')
        else:  
            if product_type =='default':
                pass
            else:
                [norm_model,default_C_list]=colse_model_default_carbon(norm_model,met_C)
            #print(norm_model.optimize(),default_C_list,substrate_id,product_id)

        if substrate_id.startswith('SK_'):
            if 'SK_'+substrate_id not in norm_model.reactions:
                substract_add=norm_model.metabolites.get_by_id(substrate_m_id)
                sink_reaction=norm_model.add_boundary(substract_add, type="sink")
                substrate_id=sink_reaction.id
                
        if substrate_uptake_rate:
            if re.search('_reverse',substrate_id):
                try:
                    norm_model.reactions.get_by_id(substrate_id).upper_bound =substrate_uptake_rate
                except:
                    pass
            else:
                norm_model.reactions.get_by_id(substrate_id).lower_bound =-substrate_uptake_rate
    return input_json,norm_model
def medium(model):
    def is_active(reaction):
        """Determine if a boundary reaction permits flux towards creating
        metabolites
        """

        return (bool(reaction.products) and (reaction.upper_bound > 0)) or (
            bool(reaction.reactants) and (reaction.lower_bound < 0)
        )

    def get_active_bound(reaction):
        """For an active boundary reaction, return the relevant bound"""
        if reaction.reactants:
            return -reaction.lower_bound
        elif reaction.products:
            return reaction.upper_bound

    return {
        rxn.id: get_active_bound(rxn) for rxn in find_boundary_types(model,'exchange') if is_active(rxn)
    }
def build_model_from_FBA_result_cb_plus2(model_comparison_results,typelist,rxn_switch,bigg_models_metabolites,coordinates,substrate_id_list,product_id_list):
    """由FBA结果构建模型
    Arguments
    ----------
    * model_pfba_solution: 模型计算结果.
    :return: cobra model.
    """
    #print(substrate_id_list,product_id_list)
    metlink_model = Model('metlink_model')
    new_substrate_id_list=[]
    new_product_id_list=[]
    substrate_sub_list=[]
    product_sub_list=[]
    for eachproduct in product_id_list:
        #Demand反应
        eachproduct =eachproduct.replace('__D_e','_e').replace('__L_e','_e')
        if eachproduct.startswith('DM_'):
            eachproduct=eachproduct.replace('DM_','')
        product_sub=eachproduct.replace('_e','').replace('_p','').replace('_c','')
        product_sub_list.append(product_sub)
        new_product_id_list.append(eachproduct)
        
    for eachsubstrate in substrate_id_list:   
        eachsubstrate =eachsubstrate.replace('__D_e','_e').replace('__L_e','_e')
        #Sink
        if eachsubstrate.startswith('SK_'):
            eachsubstrate=eachsubstrate.replace('SK_','')
        elif re.search('_reverse',eachsubstrate):
            eachsubstrate =eachsubstrate.replace('_reverse','')
        substrate_sub=eachsubstrate.replace('_e','').replace('_p','').replace('_c','')
        new_substrate_id_list.append(eachsubstrate)
        substrate_sub_list.append(substrate_sub)

    model_comparison_results=model_comparison_results.fillna('')
    #print(model_comparison_results)
    for index, row in model_comparison_results.iterrows():
        new_equ=row['equ'].replace('__D_e','_e').replace('__L_e','_e')
        #print(substrate_id,product_id)
        uni_str_line=get_reaction_uni_met(new_equ)    
        if len(uni_str_line)==1 and uni_str_line[0] in substrate_sub_list:  
            pass
        elif len(uni_str_line)==1 and uni_str_line[0] in product_sub_list:  
            pass        
        else:   
            metlink_model=get_model_met_from_reaction_string(metlink_model,new_equ,fwd_arrow='-->', rev_arrow='<--', reversible_arrow='<=>', term_split='+')
            reaction = Reaction(index) 
            reaction.name=row['reaction_name']
            reaction_ano='Reaction id: <br/>'+str(index)+'<hr/>'+'Reaction name: <br/>'+str(row['reaction_name'])+'<hr/>'+\
                        'Reaction fluxes: <br/>'+str(row['fc'])+'<hr/>'+'Reaction equation: <br/>'+\
                        str(row['equ'])+'<hr/>'+'Reaction gene_reaction_rule: <br/>'+str(row['gene'])+'<hr/>'
            reaction.notes['map_info'] = {}
            reaction.notes['map_info']['annotation']=reaction_ano
            if rxn_switch == 'T':
                reaction.notes['map_info']['display_name'] = index
            elif rxn_switch == 'F':
                reaction.notes['map_info']['display_name'] = row['reaction_name']

            if coordinates:#坐标轴赋值
                if 'reactions' in coordinates.keys():
                    if index in coordinates['reactions']:
                        reaction.notes['map_info']['x'] = float(coordinates['reactions'][index]['x'])
                        reaction.notes['map_info']['y'] = float(coordinates['reactions'][index]['y'])
            #d.notes.map_info.group   
            #[undefined, 'ko', 1, 2, 3, 4, 5, 6, 7, 8]
            for i in range(len(typelist)):           
                if row['reaction type'] ==typelist[i]:
                    reaction.notes['map_info']['group']=i+1
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

    model_met_list=[]
    for eachm in metlink_model.metabolites:
        model_met_list.append(eachm.id)           
    #print(substrate_id_list,product_id_list) 

    for i in range(len(substrate_id_list)):
        if substrate_id_list[i].endswith('_e'):
            tmp_p="_e".join(substrate_id_list[i].split('_e')[:-1])+'_p'
            tmp_c="_e".join(substrate_id_list[i].split('_e')[:-1])+'_c'
            if substrate_id_list[i]=='glc__D_e' or substrate_id_list[i]=='glc_D_e':
                if tmp_c in model_met_list:
                    substrate_id_list[i] = tmp_c    
                elif tmp_p in model_met_list:
                    substrate_id_list[i] = tmp_p   
            else:
                if tmp_c in model_met_list:
                    substrate_id_list[i] = tmp_c    
                elif tmp_p in model_met_list:
                    substrate_id_list[i] = tmp_p 

    for j in range(len(product_id_list)):                
        if product_id_list[j].endswith('_e'):
            tmp_p="_e".join(product_id_list[j].split('_e')[:-1])+'_p'
            tmp_c="_e".join(product_id_list[j].split('_e')[:-1])+'_c'
            if tmp_c in model_met_list:
                product_id_list[j] = tmp_c    
            elif tmp_p in model_met_list:
                product_id_list[j] = tmp_p 

    display_name_format = (lambda met: re.sub('__[D,L]', '', met.id[:-2].upper()))    
    #bigg_models_metabolites=bigg_models_metabolites.fillna('')
    for eachmet in metlink_model.metabolites:
        if re.search('_c',eachmet.id):
            eachmet.compartment='c'
        elif re.search('_p',eachmet.id):
            eachmet.compartment='p'
        elif re.search('_e',eachmet.id):
            eachmet.compartment='e'  
        eachmet.notes['map_info'] = {}
        show_met=''
        if eachmet.id in list(bigg_models_metabolites.index):
            show_met=bigg_models_metabolites.loc[eachmet.id,'name']
        else:
            new_met="_e".join(eachmet.id.split('_e')[:-1])+'__D_e'  
            if new_met in list(bigg_models_metabolites.index):
                show_met=bigg_models_metabolites.loc[new_met,'name']
            else:
                new_met="_e".join(eachmet.id.split('_e')[:-1])+'__L_e' 
                if new_met in list(bigg_models_metabolites.index):
                    show_met=bigg_models_metabolites.loc[new_met,'name']
                else:
                    print(eachmet.id)
        eachmet.name='Metabolite id: <br/>'+str(eachmet.id)+'<hr/>'+'Metabolite name: <br/>'+str(show_met)+'<hr/>'
        #eachmet.notes['map_info']['annotation']=metabolite_ano
        new_list=['nadh_c','nad_c','atp_c','atp_p', 'atp_e','C00002', 'C00002[c]','C00002[p]','C00002[e]','ATP', 'MNXM3', 'cpd00002', 'cpd00002_c0','cpd00002_p0','cpd00002_e0', 'HMDB00538', 'ZKHQWZAMYRWXGA-KQYNXXCUSA-J','adp_c','adp_p','adp_e','C00008','C00008[c]', 'ADP', 'MNXM7','cpd00008','cpd00008_c0','cpd00008_p0','cpd00008_e0', 'HMDB01341', 'XTWYTFMLZFPYCI-KQYNXXCUSA-K','nadph_c','NADPH','C00005','C00005[c]','MNXM6','HMDB00221', 'HMDB00799', 'HMDB06341','cpd00005','ACFIXJIJDZMPPO-NNYOXOHSSA-J','nadph','META:NADPH','cpd00005_c0','nadp_c','NADP','MNXM5','C00006','C00006[c]','HMDB00217','cpd00006','XJLXINKUBYWONI-NNYOXOHSSA-K','META:NADP','nadp','cpd00006_c0']
        if eachmet.id in new_list:    # [ATP,'cpd00002',ADP,'cpd00008',鲜肉(鲑鱼)色]
            eachmet.notes['map_info']['color'] = "#FA8072"
        # elif met.id[:-3] in ['cpd00003', 'cpd00004', 'cpd00005', 'cpd00006']:# [nad、nadh、nadp、nadph]绿宝石
        #     met.notes['map_info']['color'] = "#008000"
        else:   #  其他为适中的碧绿色
            eachmet.notes['map_info']['color'] = "#00FA9A"  
        #print(eachmet.id) 
        if eachmet.id in new_substrate_id_list:
            eachmet.notes['map_info']['substrate'] =eachmet.id
        if eachmet.id in new_product_id_list:
            eachmet.notes['map_info']['product'] =eachmet.id
 
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

def getinputtable(fc_data,norm_model2):
    # make a new dataframe named "inputdownup", first column is "reaction_id",the same as the first column of df_sample, the second column is reaction_name, can got by norm_model.reactions.get_by_id(reaction_id).name
    inputdownup = pd.DataFrame(columns=['reaction_id','reaction_name'])
    inputdownup['reaction_id'] = fc_data['reaction_id']
    inputdownup['reaction_name'] = inputdownup['reaction_id'].apply(lambda x: norm_model2.reactions.get_by_id(x).name)
    # the thrid column is "equ" which is inputdownup['reaction_id'].apply(lambda x: norm_model.reactions.get_by_id(x).reaction)
    inputdownup['equ'] = inputdownup['reaction_id'].apply(lambda x: norm_model2.reactions.get_by_id(x).reaction)

    # the forth column is "down" which is the df_sample['FC'] where df_sample['FC']<1 and no value if  otherwise
    inputdownup['down'] = fc_data.apply(
    lambda row: row['FC'] if row['manipulation'] == 'down' else None, axis=1)
    # the fifth column is "up" which is the df_sample['FC'] where df_sample['FC']>1 and no value if  otherwise
    inputdownup['up'] = fc_data.apply(
    lambda row: row['FC'] if row['manipulation'] == 'up' else None, axis=1)

    inputdownup['none'] = fc_data.apply(
    lambda row: row['FC'] if row['manipulation'] not in ['up', 'down'] else None, axis=1)


    inputdownup['gene'] = ['']*inputdownup.shape[0]
    inputdownup['fc'] = ['']*inputdownup.shape[0]
    inputdownup['reaction type'] = ['']*inputdownup.shape[0]
    inputdownup['fc_down']= ['']*inputdownup.shape[0]
    inputdownup['fc_up']= ['']*inputdownup.shape[0]
    inputdownup['fc_none']= ['']*inputdownup.shape[0]
    # using for loop: for each line in inputdownup['gene'], if inputdownup['down'] is not None, then inputdownup['gene'] is 'down':model.reaction.get_by_id(inputdownup['reaction_id']).gpr, if inputdownup['up'] is not None, then inputdownup['gene'] is 'up':model.reaction.get_by_id(inputdownup['reaction_id']).gpr, if both not none then use 'and' to connect them, if .gpr is None

    for i in range(len(inputdownup)):
        thestr=''
        fcstr=''
        reaction_type=''
        if not pd.isna(inputdownup['down'][i]):
            downstr = norm_model2.reactions.get_by_id(inputdownup['reaction_id'][i]).gene_reaction_rule
            # if downstr is '' then downstr is 'nan'
            if downstr == '':
                downstr = 'nan'
            thestr = thestr + 'Down: '+downstr
            fcstr=fcstr+'Down: '+str(-inputdownup['down'][i])
            reaction_type=reaction_type+'Down'
        else:
            downstr = ''

        if not pd.isna(inputdownup['up'][i]):
            upstr = norm_model2.reactions.get_by_id(inputdownup['reaction_id'][i]).gene_reaction_rule
            if upstr == '':
                upstr = 'nan'
            if thestr != '':
                thestr = thestr + '; '
                fcstr=fcstr+'; '
                reaction_type=reaction_type+' and '
            thestr = thestr + 'Up: '+upstr
            fcstr=fcstr+'Up: '+str(inputdownup['up'][i])
            reaction_type=reaction_type+'Up'
        else:
            upstr = ''

        if not pd.isna(inputdownup['none'][i]):
            nonestr = norm_model2.reactions.get_by_id(inputdownup['reaction_id'][i]).gene_reaction_rule
            if nonestr == '':
                nonestr = 'nan'
            if thestr != '':
                thestr = thestr + '; '
                fcstr=fcstr+'; '
                reaction_type=reaction_type+' and '
            thestr = thestr + 'none: '+nonestr
            fcstr=fcstr+'none: '+str(abs(inputdownup['none'][i]))
            reaction_type=reaction_type+'none'
        else:
            upstr = ''


        # using ';' to connect "Down:"+downstr and "Up:"+upstr and save it to inputdownup['up'][i]
        inputdownup['gene'][i] =thestr
        inputdownup['fc'][i]=fcstr
        inputdownup['reaction type'][i]=reaction_type
        inputdownup['fc_down'][i]=inputdownup['down'][i]
        inputdownup['fc_up'][i]=inputdownup['up'][i]
        inputdownup['fc_none'][i]=inputdownup['none'][i]
    # if inputdownup['down'][i] is not 'nan'

    # remove two columns "down" and "up"
    inputdownup.drop(['down','up','none'],axis=1,inplace=True)  
    # set the index of the dataframe to the reaction_id
    inputdownup.set_index('reaction_id',inplace=True)

    # export csv
    # inputdownup.to_csv('/hpcfs/fhome/xuwenqi/project/OptMetarget/OUTPUT/inputdownup.csv',index=False)
    return inputdownup
