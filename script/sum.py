# coding:utf-8
"""
Author: Jingyi Cai, Wenqi Xu (2024-2026)
Function: Merge upregulation and downregulation targets from the finished methods, attach pathway and reaction information, and write the combined tables.
Input: python sum.py <model_dir> <task_dir> <map_dir> <results_dir> <task_id>. Reads each method's results.xlsx or output.json, INPUT/task/reaction_groups_mapping.csv, and the model named in the task.
Output: UP_modification.tsv, down_modification.tsv, up.json, down.json, flux.json, and map.json in <results_dir>/<task_id>/. Enzyme methods also write UP_modification_E.tsv and down_modification_E.tsv when those results exist.
"""
import cobra
import json
import pandas as pd
import numpy as np
import sys
import os
import re
import shutil

_METHOD_FOLDER_MAP = {
    'loopless_optforce_MUST': 'OptForce',
    'FSEOF': 'FSEOF',
    'iBridge': 'iBridge',
}
_GEM_METHODS_WITH_VIZ = frozenset(_METHOD_FOLDER_MAP.keys())
def read_method_results(path, method_name, sheet_name=None, filter_col=None, filter_val=None):
    """Helper function to read and filter method results"""
    if not os.path.exists(path):
        return None, None
    try:
        kwargs = {'sheet_name': sheet_name} if sheet_name else {}
        df = pd.read_excel(path, **kwargs).iloc[:, 1:]
        up = df[df[filter_col] == filter_val[0]]
        down = df[df[filter_col] == filter_val[1]]
        return up, down
    except Exception as e:
        print(f"Error reading {method_name} results: {e}")
        return None, None

def read_file(path_results, path_task, path_model, path_map, taskname):
    """Read and process all method results with improved error handling"""
    try:
        with open(os.path.join(path_task, taskname) + '.json', encoding='utf-8') as fp:
            inputdic = json.load(fp)
        method = inputdic['taskname']
        if isinstance(method, str):
            method = [method]
    except Exception as e:
        print(f"Error reading task configuration: {e}")
        return None, None, None, None, None, None, None, None, None, None, []

    # Define method configurations
    method_config = {
        'loopless_optforce_MUST': {
            'path': os.path.join(path_results, taskname, 'OptForce', 'results.xlsx'),
            'filter_col': 'manipulations',
            'filter_val': ('Up', 'down')
        },
        'FSEOF': {
            'path': os.path.join(path_results, taskname, 'FSEOF', 'results.xlsx'),
            'filter_col': 'manipulation',
            'filter_val': ('up', 'down')
        },
        'iBridge': {
            'path': {
                'up': os.path.join(path_results, taskname, 'iBridge', 'Application_final_up_summary.txt'),
                'down': os.path.join(path_results, taskname, 'iBridge', 'Application_final_down_summary.txt')
            },
            'delimiter': '\t'
        },
        'E_FSEOF': {
            'path': os.path.join(path_results, taskname, 'E_FSEOF', 'results.xlsx'),
            'sheet_name': 'test',
            'filter_col': 'manipulations',
            'filter_val': ('up', 'down')
        },
        'E_OptForce': {
            'path': os.path.join(path_results, taskname, 'E_OptForce', 'results.xlsx'),
            'sheet_name': 'MUST',
            'filter_col': 'manipulations',
            'filter_val': ('Up', 'down')
        },
        'llm': {
            'path': os.path.join(path_results, taskname, 'llm', 'output.json'),
            'type': 'json_genes'
        }
    }
    
    results = {}
    valid_methods = []
    
    for m in method:
        if m not in method_config:
            continue
        config = method_config[m]
        try:
            if m == 'iBridge':
                if os.path.exists(config['path']['up']) and os.path.exists(config['path']['down']):
                    results['ibridge_up'] = pd.read_csv(config['path']['up'], delimiter=config['delimiter'])
                    results['ibridge_down'] = pd.read_csv(config['path']['down'], delimiter=config['delimiter'])
                    valid_methods.append(m)
            elif m == 'llm':
                llm_path = config['path']
                if os.path.exists(llm_path):
                    with open(llm_path, encoding='utf-8') as fp:
                        llm_data = json.load(fp)
                    genes = llm_data.get('genes', {})
                    rows = []
                    for gene_key, info in genes.items():
                        rows.append({
                            'gene': info.get('Gene_Name', gene_key),
                            'Enzyme_Name': info.get('Enzyme_Name', 'N/A'),
                            'Reaction_Formula': info.get('Reaction Formula', 'N/A'),
                            'Modification': info.get('Modification', ''),
                            'Fold_change': float(info.get('Fold_change', 0)),
                            'Rationale': info.get('Rationale', '')
                        })
                    if rows:
                        df_llm = pd.DataFrame(rows)
                        mod_upper = df_llm['Modification'].str.upper()
                        results['llm_up'] = df_llm[mod_upper == 'UP'].copy()
                        results['llm_down'] = df_llm[mod_upper.isin(['DOWN', 'KO'])].copy()
                        valid_methods.append('llm')
            else:
                up, down = read_method_results(
                    config['path'],
                    m,
                    config.get('sheet_name', None),
                    config['filter_col'],
                    config['filter_val']
                )
                if up is not None and down is not None:
                    results[f"{m.lower()}_up"] = up
                    results[f"{m.lower()}_down"] = down
                    valid_methods.append(m)
        except Exception as e:
            print(f"Error processing {m} results: {e}")
            continue
    
    return (
        results.get('loopless_optforce_must_up', pd.DataFrame()),
        results.get('loopless_optforce_must_down', pd.DataFrame()),
        results.get('fseof_up', pd.DataFrame()),
        results.get('fseof_down', pd.DataFrame()),
        results.get('ibridge_up', pd.DataFrame()),
        results.get('ibridge_down', pd.DataFrame()),
        results.get('e_optforce_up', pd.DataFrame()),
        results.get('e_optforce_down', pd.DataFrame()),
        results.get('e_fseof_up', pd.DataFrame()),
        results.get('e_fseof_down', pd.DataFrame()),
        results.get('llm_up', pd.DataFrame()),
        results.get('llm_down', pd.DataFrame()),
        valid_methods
    )


def normalized_format(optf_up, optf_down, fseof_up, fseof_down, ibridge_up, ibridge_down, etoptf_up, etoptf_down, etfseof_up, etfseof_down, llm_up, llm_down, method, use_mean=False):
    # Initialize all variables to None
    optf_up_sorted = optf_down_sorted = fseof_up_sorted = fseof_down_sorted = ibridge_up_sorted = ibridge_down_sorted = etoptf_up_sorted = etoptf_down_sorted = etfseof_up_sorted = etfseof_down_sorted = llm_up_sorted = llm_down_sorted = None
    print(method)
    
    if 'loopless_optforce_MUST' in method:
        # Step 1: Sort optf_up by 'range_change' in descending order
        optf_up_sorted = optf_up.sort_values(by='range_change', ascending=False)
        max_range_change = optf_up_sorted['range_change'].mean() if use_mean else optf_up_sorted['range_change'].max()
        optf_up_sorted['normalized_range_change'] = optf_up_sorted['range_change'] / max_range_change
        optf_up_sorted['FC'] = abs(optf_up_sorted['normalized_range_change']).round(3)

        # DOWN
        optf_down_sorted = optf_down.sort_values(by='range_change', ascending=False)
        max_range_change = optf_down_sorted['range_change'].mean() if use_mean else optf_down_sorted['range_change'].max()
        optf_down_sorted['normalized_range_change'] = optf_down_sorted['range_change'] / max_range_change
        optf_down_sorted['FC'] = abs(optf_down_sorted['normalized_range_change']).round(3)

    if 'FSEOF' in method:
        fseof_up_sorted = fseof_up.sort_values(by='result', ascending=False)
        max_range_change = fseof_up_sorted['result'].mean() if use_mean else fseof_up_sorted['result'].max()
        fseof_up_sorted['normalized_range_change'] = fseof_up_sorted['result'] / max_range_change
        fseof_up_sorted['FC'] = abs(fseof_up_sorted['normalized_range_change']).round(3)

        # DOWN
        fseof_down_sorted = fseof_down.sort_values(by='result', ascending=False)
        max_range_change = fseof_down_sorted['result'].mean() if use_mean else fseof_down_sorted['result'].max()
        fseof_down_sorted['normalized_range_change'] = fseof_down_sorted['result'] / max_range_change
        fseof_down_sorted['FC'] = abs(fseof_down_sorted['normalized_range_change']).round(3)
        
        # Filter out rows where FC > 50
        fseof_down_sorted = fseof_down_sorted[fseof_down_sorted['FC'] <= 50]

    if 'iBridge' in method:
        ibridge_up['difference'] = abs(ibridge_up['Positive score']) - abs(ibridge_up['Negative score'])
        ibridge_up_sorted = ibridge_up.sort_values(by='difference', ascending=False)
        max_range_change = ibridge_up_sorted['difference'].mean() if use_mean else ibridge_up_sorted['difference'].max()
        ibridge_up_sorted['normalized_range_change'] = ibridge_up_sorted['difference'] / max_range_change
        ibridge_up_sorted['FC'] = abs(ibridge_up_sorted['normalized_range_change']).round(3)
        ibridge_up_sorted = ibridge_up_sorted.rename(columns={'Reaction': 'reaction'})

        # DOWN
        ibridge_down['difference'] = ibridge_down['Positive score'] - abs(ibridge_down['Negative score'])
        ibridge_down_sorted = ibridge_down.sort_values(by='difference', ascending=False)
        max_range_change = ibridge_down_sorted['difference'].mean() if use_mean else ibridge_down_sorted['difference'].max()
        ibridge_down_sorted['normalized_range_change'] = ibridge_down_sorted['difference'] / max_range_change
        ibridge_down_sorted['FC'] = abs(ibridge_down_sorted['normalized_range_change']).round(3)
        ibridge_down_sorted = ibridge_down_sorted.rename(columns={'Reaction': 'reaction'})

    if 'E_OptForce' in method:
        etoptf_up_sorted = etoptf_up.sort_values(by='Fold_change', ascending=False)
        max_range_change = etoptf_up_sorted['Fold_change'].mean() if use_mean else etoptf_up_sorted['Fold_change'].max()
        etoptf_up_sorted['normalized_range_change'] = etoptf_up_sorted['Fold_change'] / max_range_change
        etoptf_up_sorted['FC'] = abs(etoptf_up_sorted['normalized_range_change']).apply(lambda x: f"{x:.3e}")

        # DOWN
        etoptf_down_sorted = etoptf_down.sort_values(by='Fold_change', ascending=False)
        max_range_change = etoptf_down_sorted['Fold_change'].mean() if use_mean else etoptf_down_sorted['Fold_change'].max()
        etoptf_down_sorted['normalized_range_change'] = etoptf_down_sorted['Fold_change'] / max_range_change
        etoptf_down_sorted['FC'] = abs(etoptf_down_sorted['normalized_range_change']).apply(lambda x: f"{x:.3e}")

    # ET_FSEOF processing
    if 'E_FSEOF' in method:
        selected_columns = etfseof_up.iloc[:, 1:11]
        etfseof_up['flux'] = selected_columns.mean(axis=1)

        data_slice = etfseof_up.iloc[:, 1:11]
        sum_first_three = data_slice.iloc[:, :3].sum(axis=1)
        sum_last_three = data_slice.iloc[:, -3:].sum(axis=1)
        result = ((sum_last_three - sum_first_three) / sum_first_three).round(3)
        etfseof_up['range_fc'] = result
        etfseof_up['range_fc'] = etfseof_up['range_fc'].apply(
            lambda x: 1/x if 0 < x < 1 else x
        )
        # 仅保留 range_fc 值大于 2 的行
        etfseof_up = etfseof_up[etfseof_up['range_fc'] > 2].reset_index(drop=True)
        etfseof_up_sorted = etfseof_up.sort_values(by='range_fc', ascending=False)
        selected_columns_down = etfseof_down.iloc[:, 1:11]
        etfseof_down['flux'] = selected_columns_down.mean(axis=1)
        data_slice = etfseof_down.iloc[:, 1:11]
        sum_first_three = data_slice.iloc[:, :3].sum(axis=1)
        sum_last_three = data_slice.iloc[:, -3:].sum(axis=1)
        result = abs((sum_last_three - sum_first_three) / sum_first_three).round(3)
        etfseof_down['range_fc'] = result
        etfseof_down['range_fc'] = etfseof_down['range_fc'].apply(
            lambda x: 1/x if 0 < x < 1 else x
        )
        # 仅保留 range_fc 值大于 2 的行
        etfseof_down = etfseof_down[etfseof_down['range_fc'] > 2].reset_index(drop=True)
        etfseof_down_sorted = etfseof_down.sort_values(by='range_fc', ascending=False)

        max_range_change = etfseof_up_sorted['range_fc'].mean() if use_mean else etfseof_up_sorted['range_fc'].max()
        etfseof_up_sorted['normalized_range_change'] = etfseof_up_sorted['range_fc'] / max_range_change
        etfseof_up_sorted['FC'] = (etfseof_up_sorted['normalized_range_change']).apply(lambda x: f"{x:.3e}")

        max_range_change = etfseof_down_sorted['range_fc'].mean() if use_mean else etfseof_down_sorted['range_fc'].max()
        etfseof_down_sorted['normalized_range_change'] = etfseof_down_sorted['range_fc'] / max_range_change
        etfseof_down_sorted['FC'] = (etfseof_down_sorted['normalized_range_change']).apply(lambda x: f"{x:.3e}")

    # LLM: assign FC=1.0 for all predictions (Fold_change is often 0 from LLM)
    if 'llm' in method:
        if llm_up is not None and not llm_up.empty:
            llm_up_sorted = llm_up.copy()
            fc = llm_up_sorted['Fold_change'].clip(lower=1.0)
            norm = fc.max() if not use_mean else fc.mean()
            llm_up_sorted['FC'] = (fc / norm).round(3) if norm > 0 else 1.0
        if llm_down is not None and not llm_down.empty:
            llm_down_sorted = llm_down.copy()
            fc = llm_down_sorted['Fold_change'].clip(lower=1.0)
            norm = fc.max() if not use_mean else fc.mean()
            llm_down_sorted['FC'] = (fc / norm).round(3) if norm > 0 else 1.0

    return optf_up_sorted, optf_down_sorted, fseof_up_sorted, fseof_down_sorted, ibridge_up_sorted, ibridge_down_sorted, etoptf_up_sorted, etoptf_down_sorted, etfseof_up_sorted, etfseof_down_sorted, llm_up_sorted, llm_down_sorted


def merged(method, optf_up_sorted=None, optf_down_sorted=None, 
           fseof_up_sorted=None, fseof_down_sorted=None, 
           ibridge_up_sorted=None, ibridge_down_sorted=None):
    """
    动态合并上调和下调DataFrame，根据method动态决定合并逻辑。

    Args:
        method (list): 包含所需方法的列表，如 ['optforce', 'fseof'] 或 ['optforce', 'ibridge', 'fseof']。
        optf_up_sorted (pd.DataFrame): optforce上调表。
        optf_down_sorted (pd.DataFrame): optforce下调表。
        fseof_up_sorted (pd.DataFrame): fseof上调表。
        fseof_down_sorted (pd.DataFrame): fseof下调表。
        ibridge_up_sorted (pd.DataFrame): ibridge上调表。
        ibridge_down_sorted (pd.DataFrame): ibridge下调表。

    Returns:
        pd.DataFrame, pd.DataFrame: 合并后的上调和下调DataFrame。
    """
    # 映射method到对应的DataFrame
    up_dfs = {
        'OptForce': optf_up_sorted,  # 直接使用OptForce，不区分loopless_optforce_MUST
        'FSEOF': fseof_up_sorted,
        'iBridge': ibridge_up_sorted
    }
    down_dfs = {
        'OptForce': optf_down_sorted,  # 直接使用OptForce，不区分loopless_optforce_MUST
        'FSEOF': fseof_down_sorted,
        'iBridge': ibridge_down_sorted
    }

    # 合并的通用函数
    def merge_dfs(dfs_dict):
        merged_df = None
        suffix_map = {key: f"FC_{key}" for key in dfs_dict.keys()}  # 定义后缀
        
        for key, df in dfs_dict.items():
            if df is not None:  # 确保表格存在
                if merged_df is None:
                    merged_df = df[['reaction', 'FC']].rename(columns={'FC': suffix_map[key]})
                else:
                    merged_df = merged_df.merge(
                        df[['reaction', 'FC']].rename(columns={'FC': suffix_map[key]}),
                        on='reaction',
                        how='outer'
                    )
        
        # 如果没有表格可合并，返回空DataFrame
        if merged_df is None:
            return pd.DataFrame(columns=['reaction', 'Score_sum', 'target_number'])
        
        # 填充空值为0
        merged_df.fillna({col: 0 for col in suffix_map.values()}, inplace=True)
        
        # 计算FC总和
        merged_df['Score_sum'] = merged_df[list(suffix_map.values())].sum(axis=1)
        
        # 生成target_number列
        merged_df['target_number'] = merged_df.apply(
            lambda row: ', '.join(
                method_name for method_name, col in suffix_map.items()
                if row[col] > 0
            ),
            axis=1
        )
        
        return merged_df

    # 替换 method 中的 'loopless_optforce_MUST' 为 'OptForce'
    method = ['OptForce' if m == 'loopless_optforce_MUST' else m for m in method]
    
    # 获取对应的上调和下调DataFrame子集
    up_dfs_subset = {key: up_dfs[key] for key in method if key in up_dfs}
    down_dfs_subset = {key: down_dfs[key] for key in method if key in down_dfs}
    
    # 合并上调和下调表格
    merged_df_up = merge_dfs(up_dfs_subset)
    merged_df_down = merge_dfs(down_dfs_subset)
    merged_df_up = merged_df_up.sort_values(by='Score_sum', ascending=False)
    merged_df_down = merged_df_down.sort_values(by='Score_sum', ascending=False)
    
    return merged_df_up, merged_df_down

def merge_and_process(etfseof_sorted, etoptf_sorted):
    """
    合并 FSEOF 和 OptForce 表格，处理上调或下调数据。
    
    参数:
    - etfseof_sorted: 已排序的 FSEOF 数据表（包含 'gene' 和 'FC' 列）
    - etoptf_sorted: 已排序的 OptForce 数据表（包含 'gene' 和 'FC' 列）
    
    返回:
    - 处理后的 DataFrame，包含计算后的 FC 总和和 target_number 列
    """
    
    # 给每个表添加方法标记
    etfseof_sorted['method'] = 'E_FSEOF'
    etoptf_sorted['method'] = 'E_OptForce'

    # 重命名 FC 列为方法特定的列名
    etfseof_sorted.rename(columns={'FC': 'FC_E_FSEOF'}, inplace=True)
    etoptf_sorted.rename(columns={'FC': 'FC_E_OptForce'}, inplace=True)

    # 按 gene 列排序
    etfseof_sorted.sort_values(by='gene', inplace=True)
    etoptf_sorted.sort_values(by='gene', inplace=True)

    # 合并两个表格，按 gene 匹配
    merged_df = pd.merge(
        etfseof_sorted[['gene', 'FC_E_FSEOF']], 
        etoptf_sorted[['gene', 'FC_E_OptForce']],
        on='gene',
        how='outer'
    )

    # 填充 NaN 值为 0
    merged_df.fillna({'FC_E_FSEOF': 0, 'FC_E_OptForce': 0}, inplace=True)

    # 确保 FC 列是数值类型
    merged_df['FC_E_FSEOF'] = pd.to_numeric(merged_df['FC_E_FSEOF'], errors='coerce')
    merged_df['FC_E_OptForce'] = pd.to_numeric(merged_df['FC_E_OptForce'], errors='coerce')

    # 填充 NaN 值为 0
    merged_df.fillna({'FC_E_FSEOF': 0, 'FC_E_OptForce': 0}, inplace=True)

    # 计算 FC 总和
    merged_df['Score_sum'] = merged_df['FC_E_FSEOF'] + merged_df['FC_E_OptForce']

    # 生成 target_number 列，记录方法名称
    merged_df['target_number'] = merged_df.apply(
        lambda row: ', '.join(
            algo for algo, fc in {
                'E_FSEOF': row['FC_E_FSEOF'],
                'E_OptForce': row['FC_E_OptForce']
            }.items() if fc > 0
        ),
        axis=1
    )

    # 排序
    merged_df = merged_df.sort_values(by='Score_sum', ascending=False)
    
    return merged_df

def add_reactions_to_dataframe(df, model, expression_col='gene', output_col='reactions'):
    """
    根据 DataFrame 中的基因表达式列，找到相关反应并添加到新列中，结果以逗号分隔的字符串格式存储。
    
    参数:
    - df: 包含基因表达式的 DataFrame
    - model: COBRA 模型对象
    - expression_col: 包含基因表达式的列名
    - output_col: 输出的反应列名
    
    返回:
    - 包含反应列的 DataFrame
    """
    def get_reactions_from_expression(expression, model):
        """
        根据基因表达式获取相关反应的 ID 列表，支持解析 'and' 和 'or'。
        """
        def get_gene_reactions(gene_id):
            """
            获取单个基因 ID 对应的反应集合。
            """
            try:
                gene = model.genes.get_by_id(gene_id)
                return {reaction.id for reaction in gene.reactions}
            except KeyError:
                print(f"Gene ID '{gene_id}' not found in the model.")
                return set()
        
        # 将表达式分割为 'and' 和 'or' 逻辑组
        or_groups = [group.strip() for group in expression.split(' or ')]
        all_reactions = set()
        
        for group in or_groups:
            # 对于每个 'or' 组，解析其中的 'and' 关系
            and_genes = [gene.strip() for gene in group.split(' and ')]
            # 获取所有 'and' 相关基因的反应集合的交集
            and_reactions = get_gene_reactions(and_genes[0])
            for gene in and_genes[1:]:
                and_reactions &= get_gene_reactions(gene)
            # 将结果添加到总反应集合中
            all_reactions |= and_reactions
        
        return all_reactions
    
    # 获取反应并将其以逗号分隔的字符串格式存储
    df[output_col] = df[expression_col].apply(
        lambda i: ', '.join(get_reactions_from_expression(i, model))
    )
    
    return df



def expand_reactions_to_rows(df, reaction_col='reaction'):
    """
    根据 reactions 列的内容，将反应拆分为多行。
    
    参数:
    - df: 包含反应列的 DataFrame
    - reaction_col: 需要拆分的列名
    
    返回:
    - 拆分后的 DataFrame
    """
    # 拆分 reactions 列，并扩展为多行
    expanded_df = df.explode(reaction_col).reset_index(drop=True)
    return expanded_df




def pathway(merged_df_down,merged_df_up,path_task):
    # 读取csv文件
    tujing = pd.read_csv(os.path.join(path_task, 'reaction_groups_mapping.csv'))
    # 假设 tujing 和 merged_df_down 已加载为 DataFrame

    # 确保 tujing 数据有 'Reaction ID' 和 'Group Name' 两列
    if 'Reaction ID' in tujing.columns and 'Group Name' in tujing.columns:
        # 将 'Reaction ID' 和 'Group Name' 转换为字典，方便快速查找
        reaction_to_group = dict(zip(tujing['Reaction ID'], tujing['Group Name']))
        
        # 根据 merged_df_down 中的 reaction 匹配，添加新列 '途径'
        merged_df_up['pathway'] = merged_df_up['reaction'].map(reaction_to_group)
    else:
        print("Error: 'Reaction ID' or 'Group Name' column not found in tujing.")
    if 'Reaction ID' in tujing.columns and 'Group Name' in tujing.columns:
        # 将 'Reaction ID' 和 'Group Name' 转换为字典，方便快速查找
        reaction_to_group = dict(zip(tujing['Reaction ID'], tujing['Group Name']))
        
        # 根据 merged_df_down 中的 reaction 匹配，添加新列 '途径'
        merged_df_down['pathway'] = merged_df_down['reaction'].map(reaction_to_group)
    else:
        print("Error: 'Reaction ID' or 'Group Name' column not found in tujing.")
    # 2. 移除 pathway 列中为 'Transport, Inner Membrane' 的行
    merged_df_up = merged_df_up[merged_df_up['pathway'] != 'Transport, Inner Membrane']
    
        # 1. 确保列名 'pathway' 和 'Score_sum' 存在
    if 'pathway' in merged_df_up.columns and 'Score_sum' in merged_df_up.columns:
        # 2. 按照 pathway 分组并计算每组 Score_sum 的总和
        pathway_sums = merged_df_up.groupby('pathway')['Score_sum'].sum()
        
        # 3. 将分组后的结果按总和降序排列
        sorted_pathways = pathway_sums.sort_values(ascending=False).index
        
        # 4. 根据排序的路径顺序重新排列 merged_df_up
        merged_df_up['pathway'] = pd.Categorical(merged_df_up['pathway'], categories=sorted_pathways, ordered=True)
        # merged_df_up_sorted = merged_df_up.sort_values(by='pathway')
        merged_df_up_sorted = merged_df_up.sort_values(by=['pathway', 'Score_sum'], ascending=[True, False])
    else:
        print("Error: 'pathway' or 'Score_sum' column not found in merged_df_up.")
    # 1. 确保列名 'pathway' 和 'Score_sum' 存在
    if 'pathway' in merged_df_down.columns and 'Score_sum' in merged_df_down.columns:
        # 2. 按照 pathway 分组并计算每组 Score_sum 的总和
        pathway_sums = merged_df_down.groupby('pathway')['Score_sum'].sum()
        
        # 3. 将分组后的结果按总和降序排列
        sorted_pathways = pathway_sums.sort_values(ascending=False).index
        
        # 4. 根据排序的路径顺序重新排列 merged_df_down
        merged_df_down['pathway'] = pd.Categorical(merged_df_down['pathway'], categories=sorted_pathways, ordered=True)
        # merged_df_down_sorted = merged_df_down.sort_values(by='pathway')
        merged_df_down_sorted = merged_df_down.sort_values(by=['pathway', 'Score_sum'], ascending=[True, False])
    else:
        print("Error: 'pathway' or 'Score_sum' column not found in merged_df_down.")

        # 添加 'others' 到类别中
    merged_df_up_sorted['pathway'] = merged_df_up_sorted['pathway'].cat.add_categories('others')

    # 再替换 NaN 值为 'others'
    merged_df_up_sorted['pathway'] = merged_df_up_sorted['pathway'].fillna('others')
    merged_df_down = merged_df_down[merged_df_down['pathway'] != 'Transport, Inner Membrane']
        # 添加 'others' 到类别中
    merged_df_down_sorted['pathway'] = merged_df_down_sorted['pathway'].cat.add_categories('others')

    # 再替换 NaN 值为 'others'
    merged_df_down_sorted['pathway'] = merged_df_down_sorted['pathway'].fillna('others')

    return merged_df_up_sorted,merged_df_down_sorted


def add_uniprot_ids_to_dataframe(df, model, gene_col='gene', uniprot_col='uniprot_id'):
    """
    根据基因列为 DataFrame 添加 Uniprot ID 列。
    
    参数:
    - df: 输入的 DataFrame，包含基因列
    - model: COBRA 模型对象
    - gene_col: 基因列的列名
    - uniprot_col: 新增 Uniprot ID 列的列名
    
    返回:
    - 添加了 Uniprot ID 列的 DataFrame
    """
    uniprot_list = []
    
    for gene_expression in df[gene_col]:
        # 提取基因 ID
        gene_ids = re.findall(r'\b\w+\b', gene_expression)  # 提取基因 ID
        uniprot_ids = []
        for gene_id in gene_ids:
            try:
                gene = model.genes.get_by_id(gene_id)
                uniprot_id = gene.annotation.get('uniprot', None)
                if uniprot_id:
                    uniprot_ids.append(uniprot_id)
            except KeyError:
                print(f"Gene ID '{gene_id}' not found in the model.")
        
        # 合并多个 Uniprot ID，用逗号分隔
        uniprot_list.append(", ".join(uniprot_ids))
    
    # 添加新列到 DataFrame
    df[uniprot_col] = uniprot_list
    
    return df


def calculate_pathway_scores_with_details(df, gene_col='gene', score_col='Score_sum', pathway_col='pathway'):
    """
    按照 pathway 分类，计算每个 pathway 下的 Score_sum 总和，按总和排序，保留其他列。
    同一个基因的 Score_sum 只计算一次。

    参数:
    - df: 包含基因、分值和途径信息的 DataFrame
    - gene_col: 基因列的名称
    - score_col: 分值列的名称
    - pathway_col: 途径列的名称

    返回:
    - result_df: 计算后的 DataFrame，包含所有列并按照 Score_sum 排序
    """
    # 按 pathway 和 gene 去重
    unique_gene_df = df.drop_duplicates(subset=[pathway_col, gene_col])

    # 按 pathway 分组并计算 Score_sum 总和
    pathway_scores = unique_gene_df.groupby(pathway_col).agg({
        score_col: 'sum'
    }).rename(columns={score_col: 'Total_Score'}).reset_index()

    # 合并原始数据和总分值
    merged_df = df.merge(pathway_scores, on=pathway_col, how='left')

    # 按 Total_Score 排序
    result_df = merged_df.sort_values(by='Total_Score', ascending=False).reset_index(drop=True)

    return result_df



# 定义一个函数，用于获取 reaction 的相关信息
def get_reaction_info(model, reaction_id):
    """Get reaction info with robust error handling and validation"""
    # Initialize default values
    defaults = ('N/A', 'N/A', 'N/A')
    
    # Validate inputs
    if not model or not reaction_id or not isinstance(reaction_id, str):
        return defaults
        
    try:
        # Get reaction object
        reaction = model.reactions.get_by_id(reaction_id)
        if not reaction:
            return defaults
            
        # Extract info with validation
        gene_id = str(reaction.gene_reaction_rule) if hasattr(reaction, 'gene_reaction_rule') and reaction.gene_reaction_rule else 'N/A'
        equation = str(reaction.reaction) if hasattr(reaction, 'reaction') and reaction.reaction else 'N/A'
        rxn_name = str(reaction.name) if hasattr(reaction, 'name') and reaction.name else 'N/A'
        
        # Ensure non-empty strings
        gene_id = gene_id if gene_id.strip() else 'N/A'
        equation = equation if equation.strip() else 'N/A'
        rxn_name = rxn_name if rxn_name.strip() else 'N/A'
        
        # Validate tuple length
        result = (gene_id, equation, rxn_name)
        if len(result) != 3:
            return defaults
            
        return result
        
    except KeyError:
        print(f"Reaction ID '{reaction_id}' not found in model")
    except AttributeError as e:
        print(f"Error accessing reaction attributes: {e}")
    except Exception as e:
        print(f"Unexpected error getting reaction info: {e}")
    
    # Return defaults if any error occurs
    return defaults


def _read_json_file(path):
    if not os.path.isfile(path):
        return {}
    with open(path, encoding='utf-8') as fp:
        data = json.load(fp)
    return data if isinstance(data, dict) else {}


def _merge_json_files(paths):
    merged = {}
    for path in paths:
        merged.update(_read_json_file(path))
    return merged


def prepare_target_root_json(path_results2, method):
    """
    Aggregate per-method map/flux/up/down into task-root JSON files
    (expected by upload and frontend via output.json keys).
    """
    method_dirs = []
    for m in method:
        if m not in _GEM_METHODS_WITH_VIZ:
            continue
        folder = _METHOD_FOLDER_MAP[m]
        method_dir = os.path.join(path_results2, folder)
        if os.path.isdir(method_dir):
            method_dirs.append(method_dir)
    if not method_dirs:
        return

    map_dst = os.path.join(path_results2, 'map.json')
    if not os.path.isfile(map_dst):
        for method_dir in method_dirs:
            map_src = os.path.join(method_dir, 'map.json')
            if os.path.isfile(map_src):
                shutil.copy2(map_src, map_dst)
                break

    flux_merged = {}
    for method_dir in method_dirs:
        for flux_name in ('fcflux_map.json', 'flux.json'):
            flux_src = os.path.join(method_dir, flux_name)
            data = _read_json_file(flux_src)
            for key, value in data.items():
                if key not in ('range_lb', 'range_ub'):
                    flux_merged[key] = value
    if flux_merged:
        numeric_values = [v for v in flux_merged.values() if isinstance(v, (int, float, np.floating))]
        if numeric_values:
            mean = float(np.mean(numeric_values))
            std = float(np.std(numeric_values))
            flux_merged['range_lb'] = round(mean - 2 * std, 3)
            flux_merged['range_ub'] = round(mean + 2 * std, 3)
        with open(os.path.join(path_results2, 'flux.json'), 'w', encoding='utf-8') as fp:
            json.dump(flux_merged, fp)

    up_merged = _merge_json_files(os.path.join(d, 'up.json') for d in method_dirs)
    down_merged = _merge_json_files(os.path.join(d, 'down.json') for d in method_dirs)
    with open(os.path.join(path_results2, 'up.json'), 'w', encoding='utf-8') as fp:
        json.dump(up_merged, fp, indent=2, ensure_ascii=False)
    with open(os.path.join(path_results2, 'down.json'), 'w', encoding='utf-8') as fp:
        json.dump(down_merged, fp, indent=2, ensure_ascii=False)


def update_output_json_target_paths(path_results2):
    """Add map/flux/up/down paths to output.json (local_app.py upload step)."""
    output_json = os.path.join(path_results2, 'output.json')
    if not os.path.isfile(output_json):
        return
    dictall = _read_json_file(output_json)
    dictall['map'] = os.path.join(path_results2, 'map.json')
    dictall['flux'] = os.path.join(path_results2, 'flux.json')
    dictall['up'] = os.path.join(path_results2, 'up.json')
    dictall['down'] = os.path.join(path_results2, 'down.json')
    with open(output_json, 'w', encoding='utf-8') as fp:
        json.dump(dictall, fp, indent=4, ensure_ascii=False)


if __name__=="__main__":
    cobra.Configuration().solver = os.environ.get("OPTME_COBRA_SOLVER", "cplex")
    path_model=sys.argv[1]
    path_task=sys.argv[2]
    path_map=sys.argv[3]
    path_results=sys.argv[4]
    taskname=sys.argv[5]
    path_results2=os.path.join(path_results,taskname)

    



    if not os.path.exists(path_results2):
        os.makedirs(path_results2)

    optf_up,optf_down,fseof_up,fseof_down,ibridge_up,ibridge_down,etoptf_up,etoptf_down,etfseof_up,etfseof_down,llm_up,llm_down,method = read_file(path_results, path_task, path_model, path_map, taskname)




    optf_up_sorted, optf_down_sorted, fseof_up_sorted, fseof_down_sorted, ibridge_up_sorted, ibridge_down_sorted, etoptf_up_sorted, etoptf_down_sorted, etfseof_up_sorted, etfseof_down_sorted, llm_up_sorted, llm_down_sorted = normalized_format(optf_up, optf_down, fseof_up, fseof_down, ibridge_up, ibridge_down, etoptf_up, etoptf_down, etfseof_up, etfseof_down, llm_up, llm_down, method)
    # method = [m.replace("loopless_optforce_MUST", "OptForce") for m in method]
    print(method)
    # 对 method 进行排序以支持任意顺序
    method_gem = [m for m in method if m not in ['E_OptForce', 'E_FSEOF', 'llm']]

    method_gem = sorted(method_gem)
    # Check the size of method_gem and adjust the flow
    if len(method_gem) == 1:
        if method_gem == ['FSEOF']:
            fseof_down_sorted.rename(columns={"FC": "Score_sum"}, inplace=True)
            fseof_up_sorted.rename(columns={"FC": "Score_sum"}, inplace=True)
            fseof_up_sorted['target_number'] = 'FSEOF'
            fseof_down_sorted['target_number'] = 'FSEOF'
            merged_df_up_sorted, merged_df_down_sorted = pathway(fseof_down_sorted, fseof_up_sorted, path_task)
        if method_gem == ['loopless_optforce_MUST']:
            optf_down_sorted.rename(columns={"FC": "Score_sum"}, inplace=True)
            optf_up_sorted.rename(columns={"FC": "Score_sum"}, inplace=True)
            optf_up_sorted['target_number'] = 'OptForce'
            optf_down_sorted['target_number'] = 'OptForce'
            merged_df_up_sorted, merged_df_down_sorted = pathway(optf_down_sorted, optf_up_sorted, path_task)
        if method_gem == ['iBridge']:
            ibridge_down_sorted.rename(columns={"FC": "Score_sum"}, inplace=True)
            ibridge_up_sorted.rename(columns={"FC": "Score_sum"}, inplace=True)
            ibridge_up_sorted['target_number'] = 'iBridge'
            ibridge_down_sorted['target_number'] = 'iBridge'
            merged_df_up_sorted, merged_df_down_sorted = pathway(ibridge_down_sorted, ibridge_up_sorted, path_task)
    elif len(method_gem) > 1:
        # Existing logic for when method_gem has more than 2 methods
        if set(method_gem) == {'loopless_optforce_MUST', 'FSEOF'}:
            merged_up, merged_down = merged(
                method=method_gem,
                optf_up_sorted=optf_up_sorted,
                optf_down_sorted=optf_down_sorted,
                fseof_up_sorted=fseof_up_sorted,
                fseof_down_sorted=fseof_down_sorted,
            )
        if set(method_gem) == {'loopless_optforce_MUST', 'iBridge'}:
            merged_up, merged_down = merged(
                method=method_gem,
                optf_up_sorted=optf_up_sorted,
                optf_down_sorted=optf_down_sorted,
                ibridge_up_sorted=ibridge_up_sorted,
                ibridge_down_sorted=ibridge_down_sorted
            )
        if set(method_gem) == {'FSEOF', 'iBridge'}:
            merged_up, merged_down = merged(
                method=method_gem,
                fseof_up_sorted=fseof_up_sorted,
                fseof_down_sorted=fseof_down_sorted,
                ibridge_up_sorted=ibridge_up_sorted,
                ibridge_down_sorted=ibridge_down_sorted
            )
        if set(method_gem) == {'loopless_optforce_MUST', 'FSEOF', 'iBridge'}:
            merged_up, merged_down = merged(
                method=method_gem,
                optf_up_sorted=optf_up_sorted,
                optf_down_sorted=optf_down_sorted,
                fseof_up_sorted=fseof_up_sorted,
                fseof_down_sorted=fseof_down_sorted,
                ibridge_up_sorted=ibridge_up_sorted,
                ibridge_down_sorted=ibridge_down_sorted
            )
        
        # Check if any of the methods is in ["FSEOF", "loopless_optforce_MUST", "iBridge"]
        if any(m in ["FSEOF", "loopless_optforce_MUST", "iBridge"] for m in method):
            merged_df_up_sorted, merged_df_down_sorted = pathway(merged_down, merged_up, path_task)
    if any(m in method for m in method_gem):
        # Load model with error handling
        try:
            with open(os.path.join(path_task,taskname)+'.json',encoding='utf-8') as fp:
                inputdic=json.load(fp)
                
            model_path = os.path.join(path_model,inputdic['model'])
            if not os.path.exists(model_path):
                raise FileNotFoundError(f"Model file not found: {model_path}")
                
            model = cobra.io.load_json_model(model_path)
            
            # Get reaction info with validation
            reaction_info = []
            for reaction_id in merged_df_up_sorted['reaction']:
                info = get_reaction_info(model, reaction_id)
                if len(info) != 3:
                    print(f"Invalid reaction info for {reaction_id}: {info}")
                    info = ('N/A', 'N/A', 'N/A')
                reaction_info.append(info)
                
            # Add columns with validation
            if len(reaction_info) == len(merged_df_up_sorted):
                merged_df_up_sorted['gene_id'], merged_df_up_sorted['equation'], merged_df_up_sorted['rxn_name'] = zip(*reaction_info)
            else:
                print("Warning: Reaction info length mismatch")
                merged_df_up_sorted['gene_id'] = 'N/A'
                merged_df_up_sorted['equation'] = 'N/A'
                merged_df_up_sorted['rxn_name'] = 'N/A'
                
        except FileNotFoundError as e:
            print(f"Error: {e}")
            merged_df_up_sorted['gene_id'] = 'N/A'
            merged_df_up_sorted['equation'] = 'N/A'
            merged_df_up_sorted['rxn_name'] = 'N/A'
        except json.JSONDecodeError as e:
            print(f"Error parsing JSON file: {e}")
            merged_df_up_sorted['gene_id'] = 'N/A'
            merged_df_up_sorted['equation'] = 'N/A'
            merged_df_up_sorted['rxn_name'] = 'N/A'
        except Exception as e:
            print(f"Unexpected error: {e}")
            merged_df_up_sorted['gene_id'] = 'N/A'
            merged_df_up_sorted['equation'] = 'N/A'
            merged_df_up_sorted['rxn_name'] = 'N/A'
        columns_to_keep = ['pathway', 'reaction', 'rxn_name', 'gene_id', 'equation', 'Score_sum', 'target_number']
        merged_df_up_sorted = merged_df_up_sorted[columns_to_keep]

        # 对表格中的每个 reaction 获取相关信息
        reaction_info = merged_df_down_sorted['reaction'].apply(lambda r: get_reaction_info(model, r))
        
        # 确保每个 reaction_info 包含 3 个值
        reaction_info = [(gene, eq, name) if len((gene, eq, name)) == 3 else ('N/A', 'N/A', 'N/A')
                        for gene, eq, name in reaction_info]

        # 将结果解包并添加到新的列中
        merged_df_down_sorted['gene_id'], merged_df_down_sorted['equation'], merged_df_down_sorted['rxn_name'] = zip(*reaction_info)
        columns_to_keep = ['pathway', 'reaction', 'rxn_name', 'gene_id', 'equation', 'Score_sum', 'target_number']
        merged_df_down_sorted = merged_df_down_sorted[columns_to_keep]

        # merged_df_up_sorted.to_csv((os.path.join(path_results2, 'UP_modification.csv')))
        # merged_df_down_sorted.to_csv((os.path.join(path_results2, 'down_modification.csv')))
        # 保存为 TSV 文件
        # 保留三位小数
        merged_df_up_sorted = merged_df_up_sorted.round(3)
        merged_df_down_sorted = merged_df_down_sorted.round(3)
        merged_df_up_sorted.to_csv(os.path.join(path_results2, 'UP_modification.tsv'), sep='\t', index=False)
        merged_df_down_sorted.to_csv(os.path.join(path_results2, 'down_modification.tsv'), sep='\t', index=False)
    method_ecm = [m for m in method if m.startswith('E_')]
    method_llm = ['llm'] if 'llm' in method else []
    print(method_ecm)
    print(method)
    if any(m in method for m in method_ecm):
        print('ok')    
        with open(os.path.join(path_task,taskname)+'.json',encoding='utf-8') as fp:
            inputdic=json.load(fp)    
        model = cobra.io.load_json_model(os.path.join(path_model,inputdic['model']))
        if len(method_ecm) == 1:
            print('ok')
            if method_ecm == ['E_FSEOF']:
                updated_df_up = add_reactions_to_dataframe(etfseof_up_sorted, model, expression_col='gene', output_col='reaction')
                updated_df_down = add_reactions_to_dataframe(etfseof_down_sorted, model, expression_col='gene', output_col='reaction')
                updated_df_down.rename(columns={"FC": "Score_sum"}, inplace=True)
                updated_df_up.rename(columns={"FC": "Score_sum"}, inplace=True)
                updated_df_up['target_number'] = 'E_FSEOF'
                updated_df_down['target_number'] = 'E_FSEOF'
            if method_ecm == ['E_OptForce']:
                updated_df_up = add_reactions_to_dataframe(etoptf_up, model, expression_col='gene', output_col='reaction')
                updated_df_down = add_reactions_to_dataframe(etoptf_down, model, expression_col='gene', output_col='reaction')
                updated_df_down.rename(columns={"FC": "Score_sum"}, inplace=True)
                updated_df_up.rename(columns={"FC": "Score_sum"}, inplace=True)
                updated_df_up['target_number'] = 'E_OptForce'
                updated_df_down['target_number'] = 'E_OptForce'
            updated_df_up['reaction'] = updated_df_up['reaction'].str.split(', ')  # 确保 reactions 是列表
            updated_df_down['reaction'] = updated_df_down['reaction'].str.split(', ')  # 确保 reactions 是列表
            expanded_df_up = expand_reactions_to_rows(updated_df_up, reaction_col='reaction')
            expanded_df_down = expand_reactions_to_rows(updated_df_down, reaction_col='reaction')
            merged_df_up_sorted,merged_df_down_sorted = pathway(expanded_df_down,expanded_df_up,path_task)
            merged_df_up_sorted_with_uniprot = add_uniprot_ids_to_dataframe(merged_df_up_sorted, model, gene_col='gene', uniprot_col='uniprot_id')
        # 添加 Uniprot ID 列
            merged_df_down_sorted_with_uniprot = add_uniprot_ids_to_dataframe(merged_df_down_sorted, model, gene_col='gene', uniprot_col='uniprot_id')
            # 示例使用
            result_df_up = calculate_pathway_scores_with_details(
                merged_df_up_sorted_with_uniprot,
                gene_col='gene',
                score_col='Score_sum',
                pathway_col='pathway'
            )
            # 示例使用
            result_df = calculate_pathway_scores_with_details(
                merged_df_down_sorted_with_uniprot,
                gene_col='gene',
                score_col='Score_sum',
                pathway_col='pathway'
            )        
            result_df_up.rename(columns={"gene": "Enzyme"}, inplace=True)
            result_df_up.rename(columns={"target_number": "Methods"}, inplace=True)
            result_df.rename(columns={"gene": "Enzyme"}, inplace=True)
            result_df.rename(columns={"target_number": "Methods"}, inplace=True)
            columns_to_keep = ['pathway', 'uniprot_id', 'Enzyme', 'reaction','Score_sum', 'Methods']
            result_df_up = result_df_up[columns_to_keep]
            result_df = result_df[columns_to_keep]        
                # 保存为 TSV 文件
            result_df_up.to_csv(os.path.join(path_results2, 'UP_modification_E.tsv'), sep='\t', index=False)
            result_df.to_csv(os.path.join(path_results2, 'down_modification_E.tsv'), sep='\t', index=False)    

        elif len(method_ecm) == 2:
            if any(m in ['E_OptForce', 'E_FSEOF'] for m in method_ecm):     
                # 使用该函数处理上调数据
                merged_df_up = merge_and_process(etfseof_up_sorted, etoptf_up_sorted)
                # 使用该函数处理下调数据
                merged_df_down = merge_and_process(etfseof_down_sorted, etoptf_down_sorted) 
                updated_df_up = add_reactions_to_dataframe(merged_df_up, model, expression_col='gene', output_col='reaction')
                # 调用函数
                updated_df_down = add_reactions_to_dataframe(merged_df_down, model, expression_col='gene', output_col='reaction')
                updated_df_up['reaction'] = updated_df_up['reaction'].str.split(', ')  # 确保 reactions 是列表
                updated_df_down['reaction'] = updated_df_down['reaction'].str.split(', ')  # 确保 reactions 是列表
                expanded_df_up = expand_reactions_to_rows(updated_df_up, reaction_col='reaction')
                expanded_df_down = expand_reactions_to_rows(updated_df_down, reaction_col='reaction')
                merged_df_up_sorted,merged_df_down_sorted = pathway(expanded_df_down,expanded_df_up,path_task)
                # 添加 Uniprot ID 列
                merged_df_up_sorted_with_uniprot = add_uniprot_ids_to_dataframe(merged_df_up_sorted, model, gene_col='gene', uniprot_col='uniprot_id')
            # 添加 Uniprot ID 列
                merged_df_down_sorted_with_uniprot = add_uniprot_ids_to_dataframe(merged_df_down_sorted, model, gene_col='gene', uniprot_col='uniprot_id')
                # 示例使用
                result_df_up = calculate_pathway_scores_with_details(
                    merged_df_up_sorted_with_uniprot,
                    gene_col='gene',
                    score_col='Score_sum',
                    pathway_col='pathway'
                )
                # 示例使用
                result_df = calculate_pathway_scores_with_details(
                    merged_df_down_sorted_with_uniprot,
                    gene_col='gene',
                    score_col='Score_sum',
                    pathway_col='pathway'
                )
                result_df_up.rename(columns={"gene": "Enzyme"}, inplace=True)
                result_df_up.rename(columns={"Total_Score": "Score_Sum"}, inplace=True)
                result_df_up.rename(columns={"target_number": "Methods"}, inplace=True)
                result_df.rename(columns={"gene": "Enzyme"}, inplace=True)
                result_df.rename(columns={"Total_Score": "Score_Sum"}, inplace=True)
                result_df.rename(columns={"target_number": "Methods"}, inplace=True)
                columns_to_keep = ['pathway', 'uniprot_id', 'Enzyme', 'reaction','Score_sum', 'Methods']
                result_df_up = result_df_up[columns_to_keep]
                result_df = result_df[columns_to_keep]
                # 保留小数点后 3 位
                result_df_up['Score_sum'] = result_df_up['Score_sum'].apply(lambda x: f"{x:.3f}").astype(float)
                result_df['Score_sum'] = result_df_up['Score_sum'].apply(lambda x: f"{x:.3f}").astype(float)
                # 保存为 TSV 文件
                result_df_up.to_csv(os.path.join(path_results2, 'UP_modification_E.tsv'), sep='\t', index=False)
                result_df.to_csv(os.path.join(path_results2, 'down_modification_E.tsv'), sep='\t', index=False)    

    # ─── LLM predictions output ────────────────────────────────────────────────
    if method_llm and (llm_up_sorted is not None or llm_down_sorted is not None):
        try:
            with open(os.path.join(path_task, taskname) + '.json', encoding='utf-8') as fp:
                inputdic = json.load(fp)
            model = cobra.io.load_json_model(os.path.join(path_model, inputdic['model']))

            for df_sorted, direction, key in [
                (llm_up_sorted, 'UP', 'UP'),
                (llm_down_sorted, 'Down', 'Down')
            ]:
                if df_sorted is None or df_sorted.empty:
                    continue
                df_sorted = df_sorted.copy()
                df_sorted['Score_sum'] = df_sorted['FC']
                df_sorted['target_number'] = 'LLM'

                # Try reaction lookup via model gene
                try:
                    df_sorted = add_reactions_to_dataframe(df_sorted, model, expression_col='gene', output_col='reaction')
                    df_sorted['reaction'] = df_sorted['reaction'].str.split(', ')
                    df_sorted = expand_reactions_to_rows(df_sorted, reaction_col='reaction')
                    # Fill empty reactions with Reaction_Formula (when model gene lookup fails)
                    if 'Reaction_Formula' in df_sorted.columns:
                        empty_mask = df_sorted['reaction'].isna() | df_sorted['reaction'].astype(str).str.strip().eq('')
                        df_sorted.loc[empty_mask, 'reaction'] = df_sorted.loc[empty_mask, 'Reaction_Formula']
                    df_sorted, _ = pathway(df_sorted, df_sorted, path_task)
                    df_sorted = add_uniprot_ids_to_dataframe(df_sorted, model, gene_col='gene', uniprot_col='uniprot_id')
                    result_df = calculate_pathway_scores_with_details(df_sorted, gene_col='gene', score_col='Score_sum', pathway_col='pathway')
                    result_df.rename(columns={"gene": "Enzyme", "target_number": "Methods"}, inplace=True)
                    columns_to_keep = ['pathway', 'uniprot_id', 'Enzyme', 'reaction', 'Score_sum', 'Methods', 'Rationale']
                    result_df['Rationale'] = result_df.get('Rationale', 'N/A') if 'Rationale' in result_df.columns else 'N/A'
                    result_df = result_df[[c for c in columns_to_keep if c in result_df.columns]]
                except Exception as e:
                    print(f"LLM pathway lookup failed ({direction}): {e}, using raw output")
                    result_df = df_sorted.copy()
                    # Use Reaction_Formula as reaction content when model lookup unavailable
                    if 'Reaction_Formula' in result_df.columns:
                        result_df['reaction'] = result_df['Reaction_Formula']
                    elif 'reaction' not in result_df.columns:
                        result_df['reaction'] = 'N/A'
                    result_df.rename(columns={"gene": "Enzyme", "target_number": "Methods"}, inplace=True)
                    result_df['pathway'] = 'LLM prediction'
                    result_df['uniprot_id'] = 'N/A'
                    columns_to_keep = ['pathway', 'uniprot_id', 'Enzyme', 'Enzyme_Name', 'reaction', 'Score_sum', 'Methods', 'Rationale']
                    result_df = result_df[[c for c in columns_to_keep if c in result_df.columns]]

                fname = 'UP_modification_E.tsv' if direction == 'UP' else 'down_modification_E.tsv'
                result_df.to_csv(os.path.join(path_results2, fname), sep='\t', index=False)
                print(f"LLM {direction}: saved {len(result_df)} rows to {fname}")

        except Exception as e:
            print(f"Error processing LLM results: {e}")

    prepare_target_root_json(path_results2, method)
    update_output_json_target_paths(path_results2)

#  python sum.py '/hpcfs/fhome/xuwenqi/project/OptMetarget/INPUT/model' '/hpcfs/fhome/xuwenqi/project/OptMetarget/INPUT/task' '/hpcfs/fhome/xuwenqi/project/OptMetarget/INPUT/map' '/hpcfs/fhome/xuwenqi/project/OptMetarget/OUTPUT' 'wangry@tib.cas.cn_20241203-092527_iML1515'
