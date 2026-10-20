# coding:utf-8
"""
Author: Jingyi Cai, Wenqi Xu (2024-2026)
Function: Check the task model before target prediction. Confirms the product reaction and metabolite formulas, then writes the checked model back.
Input: python Error.py <model_dir> <task_dir> <task_id>. Reads <task_dir>/<task_id>.json and the model file named in that task.
Output: Prints a status code (100 means pass). Overwrites the task JSON and the model file in place when the check finishes.
"""
import os
import json
import os
import cobra
import sys

import difflib


def reaction_has_complete_formula(model_input, reaction_id):
    if reaction_id not in model_input.reactions:
        return False
    reaction = model_input.reactions.get_by_id(reaction_id)
    for metabolite in reaction.metabolites.keys():
        formula = getattr(metabolite, 'formula', None)
        if formula is None or str(formula).strip() == '':
            return False
    return True

def deal_dict(path_model, path_task, taskname):
    chosen = os.environ.get("OPTME_COBRA_SOLVER", "").strip()
    if chosen:
        cobra.Configuration().solver = chosen
    else:
        try:
            cobra.Configuration().solver = 'glpk'
        except Exception:
            try:
                cobra.Configuration().solver = 'cplex'
            except Exception:
                print("Warning: Could not set a solver. Using the COBRA default.")
    # 加载任务字典
    with open(os.path.join(path_task, taskname) + '.json', encoding='utf-8') as fp:
        inputdic = json.load(fp)

    # 加载模型
    if '.mat' in inputdic['model']:
        model_input = cobra.io.load_matlab_model(os.path.join(path_model, inputdic['model']))
    elif '.xml' in inputdic['model'] or '.sbml' in inputdic['model']:
        model_input = cobra.io.read_sbml_model(os.path.join(path_model, inputdic['model']))
    elif '.json' in inputdic['model']:
        model_input = cobra.io.load_json_model(os.path.join(path_model, inputdic['model']))
    else:
        model_input = cobra.io.load_model(inputdic['model'])

    # 提取产物反应名称
    productDM_name = inputdic['product'].replace("user_add_", "")
    productEX_name = None  # 初始化变量

    if 'user_add_DM_' in inputdic['product']:
        # 修改反应名称：第一个 "_" 前修改为 EX，最后一个 "_" 前修改为 e
        productEX_name = inputdic['product'].replace('user_add_DM_', 'EX_')
        # 找到最后一个 "_"，并将其前部分替换为 'e'
        productEX_name = '_'.join(productEX_name.rsplit('_', 2)[:-1]) + '_e' + productEX_name.rsplit('_', 1)[-1]

    # 获取模型中最相似的反应ID
    most_similar_reaction_id = None
    highest_similarity = 0

    if productEX_name is None:
        productEX_name = ""

    for reaction in model_input.reactions:
        # 计算反应 ID 和 productEX_name 的相似度
        similarity = difflib.SequenceMatcher(None, productEX_name, reaction.id).ratio()

        # 更新最相似的反应 ID
        if similarity > highest_similarity:
            highest_similarity = similarity
            most_similar_reaction_id = reaction.id

    # 输出最相似的反应ID和相似度
    print(f"Most similar reaction ID: {most_similar_reaction_id} with similarity: {highest_similarity}")

    # 检查是否存在 EX_产品
    product_used = None  # 初始化使用的产品为 None
    if most_similar_reaction_id and most_similar_reaction_id in model_input.reactions:
        product_used = most_similar_reaction_id
        model_input.reactions.get_by_id(productDM_name).bounds = (0, 0)
        model_input.reactions.get_by_id(product_used).bounds = (0, 1000)
    elif productDM_name in model_input.reactions:
        product_used = productDM_name

    # 更新字典中的产品信息
    inputdic['product'] = product_used
    inputdic['yield_na'] = False
    inputdic['yield_na_reason'] = ''

    return inputdic, model_input, product_used, productEX_name, productDM_name




def Error_Code(inputdic,model_input):
    class ModelError(Exception):
        def __init__(self, error_code, message):
            self.error_code = error_code
            self.message = message
            super().__init__(self.message)

        # 错误代码001：biomass不在表里
    class BiomassNotFoundError(ModelError):
        def __init__(self):
            error_code = '001'
            message = "biomass not found"
            super().__init__(error_code, message)

        # 错误代码002：底物不在表里
    class SubstrateNotFoundError(ModelError):
        def __init__(self, substrate_name):
            error_code = '002'
            message = f"substrate '{substrate_name}' not found"
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
                    message = f"product '{productEX_name}' and '{productDM_name}' not found"
                elif productEX_name:
                    message = f"product '{productEX_name}' not found"
                elif productDM_name:
                    message = f"product '{productDM_name}' not found"
            super().__init__(error_code, message)
            self.product_used = product_used 




    # 错误代码004：biomass为0，不生长
    class GrowthError(ModelError):
        def __init__(self):
            error_code = '004'
            message = "no growth"
            super().__init__(error_code, message)


    # 错误代码005：product为0，无产品
    class ProductError(ModelError):
        def __init__(self):
            error_code = '005'
            message = "no production"
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


        # Check formula completeness for substrate/product reactions and skip yield when missing
        missing_formula_msgs = []
        if not reaction_has_complete_formula(model_input, substrate_name):
            missing_formula_msgs.append(f"substrate '{substrate_name}' has metabolites without formula")
        if not reaction_has_complete_formula(model_input, product_used):
            missing_formula_msgs.append(f"product '{product_used}' has metabolites without formula")

        if missing_formula_msgs:
            inputdic['yield_na'] = True
            inputdic['yield_na_reason'] = '; '.join(missing_formula_msgs)
            print('Warning:', inputdic['yield_na_reason'])
        else:
            inputdic['yield_na'] = False
            inputdic['yield_na_reason'] = ''

        # 检查是否生长
        model_input.objective = inputdic['biomass']
        objvalue2 = model_input.optimize().objective_value
        
        if objvalue2 == 0:
            raise GrowthError()
        
        # 检查是否无产品
        model_input.reactions.get_by_id(product_used).bounds = (0, 1000)
        model_input.objective = product_used
        objvalue1 = model_input.optimize().objective_value

        if objvalue1 == 0:
            raise ProductError()
        

    except ModelError as e:
        if e.error_code is not None:   # 捕获自定义异常并处理
            error_code = e.error_code
            message = e.message
            print("Error Code:", e.error_code)
            print("Error description:", e.message)


    else:
        # 没有发生异常的情况下执行的代码
        error_code = '100'
        message = "pass"
    
    finally:
        # 无论是否发生异常，最终都会执行的代码
        print("Code:", error_code)  # 打印错误代码
        print("Description:", message)  # 打印错误描述


def write_new_json(new_inputdic,path_task,taskname):
    with open(os.path.join(path_task,taskname)+'.json',"w") as f:
        json.dump(new_inputdic,f)
        

def write_new_model(inputdic,model_input,path_model):
    if  '.mat' in inputdic['model']: 
        cobra.io.save_matlab_model(model_input,os.path.join(path_model,inputdic['model']))
    elif '.xml' in inputdic['model']:
        cobra.io.write_sbml_model(model_input,os.path.join(path_model,inputdic['model']))
    elif '.sbml' in inputdic['model']:
        cobra.io.write_sbml_model(model_input,os.path.join(path_model,inputdic['model']))   #  need add judge
    elif '.json' in inputdic['model']:
        cobra.io.save_json_model(model_input,os.path.join(path_model,inputdic['model']))
    else:
        pass
    print('write new model file success.')




if __name__=="__main__":
    path_model=sys.argv[1]
    path_task=sys.argv[2]
    taskname=sys.argv[3]

    inputdic,model_input,product_used,productEX_name,productDM_name=deal_dict(path_model,path_task,taskname)
    Error_Code(inputdic,model_input)
    write_new_json(inputdic,path_task,taskname)
    write_new_model(inputdic,model_input,path_model)
        #write_new_model(inputdic,model_input,path_model)
# python Error.py '/home/xuwenqi/sun/input1/MCS/model' '/home/xuwenqi/sun/input1/MCS/task' 'maoyf@tib.cas.cn_20230628-163545_mcs'
# python Error.py '/home/sun/optme/input1/target/model' '/home/sun/optme/input1/target/task' 'wangry@tib.cas.cn_20240305-150420_optforce'
# python Error.py '/hpcfs/fhome/xuwenqi/project/optme-master/OptMe_target/input1/target/model' '/hpcfs/fhome/xuwenqi/project/optme-master/OptMe_target/input1/target/task' 'wangry@tib.cas.cn_20241209-163548_iML1515'
