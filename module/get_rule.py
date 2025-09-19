import os
import json
import argparse

def parse_rule_file(file_path):
    rules_dict = {}
    # 检查输入文件是否存在
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The input file '{file_path}' was not found.")
    
    with open(file_path, 'r', encoding='utf-8') as file:
        # 将输入文件所有行合并为一个字符串，去除空白字符，按照//分割为多个块
        blocks = ''.join(file.readlines()).strip().split('//')
        
        for block in blocks:
            # 遍历每个块，跳过每个以#开头或为空的块
            if not block.strip() or block.startswith('#'):
                continue
            # 将块分割成多行，去除掉以#开头的注释行
            lines = [line.strip() for line in block.split('\n') if line.strip() and not line.startswith('#')]
            # 创建一个空字典
            rule = {}
            for line in lines:
                # 使用第一个出现的':'将行分割成键和值两部分
                key, value = line.split(':', 1)
                rule[key.strip()] = value.strip()
            
            # 从字典中获取必要字段
            id = rule.get('ID', None)
            if not id:
                raise KeyError(f"Missing 'ID' field in rule block: {block}")
            
            # 构建规则数据结构
            rule_data = {
                'name': rule['Name'],
                'family': rule['Family'],
                'type': rule['Type'],
                'desc': [] if rule['Desc'] == 'NA' else [rule['Desc']],
                'mode': ['a'] if ':' not in rule['Required'] else ['b'],
                'required': rule['Required'].replace(':', '#').split('#'),
                'forbidden': rule.get('Forbidden', '').split(':') if 'Forbidden' in rule else []
            }
            
            rules_dict[id] = rule_data
    
    return rules_dict

def main():
    # 添加命令行参数解析
    parser = argparse.ArgumentParser(description='解析规则文件并生成JSON格式输出')
    parser.add_argument('-i', '--input', required=True, help='输入规则文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出目录路径')
    
    args = parser.parse_args()
    
    input_file = args.input
    output_dir = args.output
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 输出文件路径
    output_file = os.path.join(output_dir, "getrule.json")
    
    try:
        # 解析规则文件
        rules = parse_rule_file(input_file)
        
        # 将结果写入JSON文件
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(rules, f, indent=2, ensure_ascii=False)
        
        print(f"✅ 规则已成功解析并保存到 {output_file}")
        
    except Exception as e:
        print(f"❌ 处理过程中出现错误: {str(e)}")

if __name__ == "__main__":
    main()