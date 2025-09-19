import json
import os
import sys
import importlib.util
from pathlib import Path

# 动态导入getrule和specpfam模块
def import_module_dynamically(module_name):
    """动态导入模块"""
    try:
        # 首先尝试直接导入
        return importlib.import_module(module_name)
    except ImportError:
        # 如果失败，尝试从当前目录导入
        current_dir = Path(__file__).parent
        module_path = current_dir / f"{module_name}.py"
        if module_path.exists():
            spec = importlib.util.spec_from_file_location(module_name, module_path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module
        else:
            raise ImportError(f"Cannot find module {module_name}")

# 导入所需模块
get_rule = import_module_dynamically('get_rule')
selfbuild_hmm = import_module_dynamically('selfbuild_hmm')

# 获取函数引用
parse_rule_file = get_rule.parse_rule_file
parse_pfam_spec = selfbuild_hmm.parse_pfam_spec

# 全局变量，将在运行时设置
filtered_result = None
result_dir = None

def set_result_dir(dir_path):
    """设置结果目录路径"""
    global result_dir
    result_dir = dir_path

def load_filtered_result(filtered_data=None):
    """加载过滤后的结果，优先使用内存中的数据，否则从文件读取
    
    Args:
        filtered_data (dict, optional): 内存中的过滤数据
    
    Returns:
        dict: 加载的数据，如果失败则返回默认结构
    """
    global result_dir
    
    # 如果提供了内存中的数据，直接使用
    if filtered_data is not None:
        print("使用内存中的过滤数据")
        return filtered_data
    
    # 否则从文件读取（保持向后兼容）
    if not result_dir:
        print("[错误] 结果目录未设置")
        return {"result": {"match": []}}
        
    processed_ipr_path = os.path.join(result_dir, "processed_ipr_domains.json")
    
    if os.path.exists(processed_ipr_path):
        print(f"从文件加载过滤数据: {processed_ipr_path}")
        with open(processed_ipr_path, "r") as f:
            return json.load(f)
    else:
        print(f"[错误] 找不到processed_ipr_domains.json文件: {processed_ipr_path}")
        return {"result": {"match": []}}

def initialize_module(result_directory, filtered_data=None):
    """初始化模块，设置结果目录并加载数据
    
    Args:
        result_directory (str): 结果目录路径
        filtered_data (dict, optional): 内存中的过滤数据
    """
    global filtered_result
    set_result_dir(result_directory)
    filtered_result = load_filtered_result(filtered_data)

# 获取当前文件所在目录的路径
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

# 全局变量，将在运行时设置
rules_dict = {}
result_spec = {}

def initialize_rules(rule_file_path):
    """初始化规则文件"""
    global rules_dict
    try:
        if not os.path.exists(rule_file_path):
            raise FileNotFoundError(f"规则文件不存在：{rule_file_path}")
            
        rules_dict = parse_rule_file(rule_file_path)
        print("[成功] 成功加载规则文件")
        return True
    except Exception as e:
        print(f"[错误] 错误：加载规则文件时发生异常 - {str(e)}")
        rules_dict = {}
        return False

def merge_results(dict1, dict2):
    """合并两个结果字典"""
    print(f"[信息] 开始合并结果...")
    print(f"dict1 中有 {len(dict1['result']['match'])} 个基因")
    print(f"dict2 中有 {len(dict2['result']['match'])} 个基因")
    
    # 创建新的match列表
    merged_matches = []
    
    # 将第一个字典的所有基因添加到merged_matches
    for gene_match in dict1["result"]["match"]:
        merged_matches.append(gene_match)
    
    # 将第二个字典的基因添加或合并到merged_matches
    merged_count = 0
    new_count = 0
    
    for gene_match in dict2["result"]["match"]:
        for gene_id, gene_data in gene_match.items():
            # 检查这个基因是否已经存在
            found = False
            for existing_match in merged_matches:
                if gene_id in existing_match:
                    print(f"[信息] 合并基因 {gene_id}")
                    # 找到了现有基因，合并域数据
                    existing_gene_data = existing_match[gene_id]
                    
                    # 两种格式现在都是 [sequence_obj, {domain: [matches]}]
                    if isinstance(existing_gene_data, list) and len(existing_gene_data) >= 2:
                        if isinstance(existing_gene_data[1], dict) and isinstance(gene_data, list) and len(gene_data) >= 2:
                            # 合并第二个元素（域字典）
                            for domain, hits in gene_data[1].items():
                                print(f"  添加域: {domain}")
                                if domain in existing_gene_data[1]:
                                    # 如果域已存在，扩展匹配列表
                                    if isinstance(existing_gene_data[1][domain], list):
                                        existing_gene_data[1][domain].extend(hits)
                                    else:
                                        existing_gene_data[1][domain] = [existing_gene_data[1][domain]] + hits
                                else:
                                    # 新域，直接添加
                                    existing_gene_data[1][domain] = hits
                    
                    merged_count += 1
                    found = True
                    break
            
            if not found:
                print(f"[信息] 添加新基因 {gene_id}")
                # 如果是新基因，创建新的字典并添加
                merged_matches.append({gene_id: gene_data})
                new_count += 1
    
    print(f"[成功] 合并完成: {merged_count} 个基因合并, {new_count} 个新基因添加")
    return {"result": {"match": merged_matches}}

# 其他函数保持不变
def get_ipr_set(gene_data):
    """从基因数据中提取IPR集合"""
    # gene_data是一个列表，第二个元素（索引1）包含IPR信息
    if len(gene_data) > 1 and isinstance(gene_data[1], dict):
        return set(ipr for ipr in gene_data[1].keys())
    return set()

def check_rule_match(ipr_set, rule_data):
    """检查IPR集合是否匹配规则"""
    # 获取规则的模式和条件
    mode = rule_data.get("mode", [])[0] if rule_data.get("mode") else None
    required = set(rule_data.get("required", []))
    forbidden = set(rule_data.get("forbidden", []))
    
    # 移除NA
    if "NA" in forbidden:
        forbidden.remove("NA")
    
    # 检查forbidden条件
    if forbidden and not forbidden.isdisjoint(ipr_set):
        return False
    
    # 根据mode检查required条件
    if mode == "a":
        # required必须是ipr_set的子集
        return required.issubset(ipr_set)
    elif mode == "b":
        # required必须与ipr_set有交集
        return not required.isdisjoint(ipr_set)
    
    return False

def classify_genes(input_dict):
    """对基因进行分类并生成结果"""
    result = {}
    
    # 处理输入字典中的每个基因
    for gene_match in input_dict["result"]["match"]:
        for gene_id, gene_data in gene_match.items():
            # 获取基因的IPR集合
            ipr_set = get_ipr_set(gene_data)
            
            # 存储匹配的所有家族
            matched_families = []
            final_rule = None
            
            # 检查每个规则
            for rule_id, rule_data in rules_dict.items():
                if check_rule_match(ipr_set, rule_data):
                    matched_families.append(rule_data)
            
            # 如果有匹配的规则
            if matched_families:
                # 选择最后一个匹配的规则作为主要规则
                final_rule = matched_families[-1]
                
                # 收集所有匹配家族的名称
                all_families = [rule["family"] for rule in matched_families]
                # 使用逗号和空格连接所有家族名称
                other_families = ", ".join(all_families)
                
                # 构建结果字典
                result[gene_id] = {
                    "name": final_rule["name"],
                    "family": final_rule["family"],
                    "type": final_rule["type"],
                    "desc": final_rule.get("desc", []),
                    "other_family": other_families
                }
    
    # 将结果写入文件（仅在直接调用时，不在process_with_data模式下）
    try:
        # 只有在result_dir未设置时才使用默认路径（向后兼容）
        if result_dir is None:
            # 构建JSON输出路径
            json_path = os.path.join(os.path.dirname(CURRENT_DIR), 'match.json')
            # 构建TBL输出路径
            tbl_path = os.path.join(os.path.dirname(CURRENT_DIR), 'match_tbl.txt')
            
            # 写入JSON文件
            with open(json_path, 'w', encoding='utf-8') as f:
                json.dump(result, f, indent=2, ensure_ascii=False)
                
            # 写入TBL文件
            with open(tbl_path, 'w', encoding='utf-8') as f:
                for gene_id, data in result.items():
                    desc_str = ';'.join(data['desc']) if data['desc'] else 'NA'
                    line = f"{gene_id}\t{data['name']}\t{data['family']}\t{data['type']}\t{desc_str}\t{data['other_family']}\n"
                    f.write(line)
                
        print("[成功] 结果已成功保存到文件")
    except Exception as e:
        print(f"[错误] 错误：保存结果时发生异常 - {str(e)}")
        
    return result

def process_with_data(result_directory, rule_file, filtered_data=None, spec_data=None, debug=False):
    """
    使用内存中的数据进行转录因子分类处理
    
    Args:
        result_directory (str): 结果目录路径
        rule_file (str): 规则文件路径
        filtered_data (dict, optional): 内存中的过滤数据
        spec_data (dict, optional): 内存中的spec数据
        debug (bool): 是否启用调试模式
    
    Returns:
        dict: 分类结果
    """
    # 初始化模块
    initialize_module(result_directory, filtered_data)
    
    # 初始化规则
    if not initialize_rules(rule_file):
        print("初始化规则文件失败")
        return None
    
    # 获取filtered_result
    if filtered_data is not None:
        current_filtered_result = filtered_data
    else:
        current_filtered_result = load_filtered_result()
    
    if current_filtered_result is None:
        print("无法获取过滤数据")
        return None
    
    # 处理spec数据
    if spec_data is not None:
        print("使用内存中的spec数据")
        current_spec_result = spec_data
    elif debug:
        # 在debug模式下尝试从文件加载spec数据
        pfamspec_file = os.path.join(result_directory, 'pfamspec.json')
        if os.path.exists(pfamspec_file):
            print(f"从文件加载spec数据: {pfamspec_file}")
            with open(pfamspec_file, 'r') as f:
                current_spec_result = json.load(f)
        else:
            print("debug模式下未找到pfamspec.json文件，使用空数据")
            current_spec_result = {"result": {"match": []}}
    else:
        # 非debug模式下使用空数据
        current_spec_result = {"result": {"match": []}}
    
    # 合并结果
    merged_dict = merge_results(current_filtered_result, current_spec_result)
    
    # 进行分类
    classification_result = classify_genes(merged_dict)
    
    return classification_result

if __name__ == "__main__":
    # 合并filtered_result和result_spec
    merged_dict = merge_results(filtered_result, result_spec)
    
    # 使用合并后的字典进行分类
    classify_genes(merged_dict)
    print("[成功] 分类完成，结果已保存为 match.json")