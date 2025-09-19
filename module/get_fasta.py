import os
import json
import argparse
from Bio import SeqIO

def get_output_dir(base_dir=None):
    """获取输出目录路径
    
    Args:
        base_dir (str, optional): 基础目录路径，如果不指定则使用项目根目录
    
    Returns:
        str: result目录的完整路径
    """
    if base_dir is None:
        # 获取当前脚本所在目录的父目录（项目根目录）
        current_dir = os.path.dirname(os.path.abspath(__file__))
        base_dir = os.path.dirname(current_dir)
    
    # 创建result目录路径
    output_dir = os.path.join(base_dir, "result")
    
    # 如果目录不存在，创建它
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    return output_dir

def read_fasta_to_dict(fasta_file):
    """
    读取FASTA文件并转换为字典格式
    
    Args:
        fasta_file (str): FASTA文件路径
    
    Returns:
        dict: 基因ID到序列的映射字典
    """
    fasta_dict = {}
    
    try:
        with open(fasta_file, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                # 使用record.id作为键，序列作为值
                fasta_dict[record.id] = str(record.seq)
        
        print(f"✅ 成功读取 {len(fasta_dict)} 个序列")
        return fasta_dict
        
    except Exception as e:
        print(f"❌ 读取FASTA文件时出错: {e}")
        return {}

def format_tf_fasta_with_classification(fasta_file, classification_result, output_file):
    """
    根据分类结果生成包含TF分类信息的标准格式FASTA文件
    
    Args:
        fasta_file (str): 输入的FASTA文件路径
        classification_result (dict): 分类结果字典，格式为 {gene_id: {name, family, type, desc, other_family}}
        output_file (str): 输出FASTA文件路径
    
    Returns:
        bool: 是否成功生成文件
    """
    try:
        # 读取FASTA文件
        fasta_dict = read_fasta_to_dict(fasta_file)
        if not fasta_dict:
            print("❌ 无法读取FASTA文件")
            return False
        
        # 生成格式化的FASTA文件，只包含有分类结果的序列
        with open(output_file, 'w') as f:
            written_count = 0
            for gene_id, sequence in fasta_dict.items():
                # 只处理在分类字典中的基因
                if gene_id in classification_result:
                    tf_info = classification_result[gene_id]
                    
                    # 构建FASTA头部信息
                    # 格式: >gene_id | family | type
                    header = f">{gene_id} | {tf_info['family']} | {tf_info['type']}"
                    
                    # 写入头部和序列
                    f.write(header + '\n')
                    f.write(sequence + '\n')
                    written_count += 1
            
            print(f"✅ 成功生成包含 {written_count} 个分类序列的FASTA文件")
            print(f"✅ 所有序列都包含TF分类信息")
            return True
            
    except Exception as e:
        print(f"❌ 生成FASTA文件时出错: {e}")
        return False

def generate_classified_fasta(fasta_file, classification_result, output_file=None, output_dir=None):
    """
    供其他模块调用的函数，生成包含TF分类信息的FASTA文件
    只包含有分类结果的序列
    
    Args:
        fasta_file (str): 输入的FASTA文件路径
        classification_result (dict): 从内存中获取的分类结果字典
        output_file (str, optional): 输出文件路径，如果不指定则自动生成
        output_dir (str, optional): 输出目录路径，主脚本可以指定项目子文件夹
    
    Returns:
        str: 输出文件路径，如果失败返回None
    """
    try:
        # 如果没有指定输出文件，自动生成路径
        if not output_file:
            input_name = os.path.splitext(os.path.basename(fasta_file))[0]
            result_dir = get_output_dir(output_dir)
            output_file = os.path.join(result_dir, f"{input_name}_tf_classified.fasta")
        
        # 调用格式化函数
        success = format_tf_fasta_with_classification(fasta_file, classification_result, output_file)
        
        if success:
            return output_file
        else:
            return None
            
    except Exception as e:
        print(f"❌ 生成分类FASTA文件时出错: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description='生成包含TF分类信息的FASTA文件')
    parser.add_argument('-i', '--input', required=True, help='输入FASTA文件路径')
    parser.add_argument('-o', '--output', help='输出FASTA文件路径（可选）')
    parser.add_argument('--classification', help='分类结果JSON文件路径（用于测试）')
    
    args = parser.parse_args()
    
    # 如果没有指定输出文件，使用默认路径
    if not args.output:
        input_name = os.path.splitext(os.path.basename(args.input))[0]
        output_dir = get_output_dir()
        args.output = os.path.join(output_dir, f"{input_name}_tf_classified.fasta")
    
    # 如果提供了分类结果文件（用于测试），从文件读取
    if args.classification:
        try:
            with open(args.classification, 'r', encoding='utf-8') as f:
                classification_result = json.load(f)
            print(f"✅ 从文件加载了 {len(classification_result)} 个分类结果")
        except Exception as e:
            print(f"❌ 读取分类结果文件时出错: {e}")
            classification_result = {}
    else:
        # 在实际使用中，这里应该从内存中获取分类结果
        # 目前为空字典，表示没有分类信息
        classification_result = {}
        print("⚠️  未提供分类结果，将生成不包含TF分类信息的FASTA文件")
    
    # 生成FASTA文件
    success = format_tf_fasta_with_classification(args.input, classification_result, args.output)
    
    if success:
        print(f"✅ FASTA文件已保存到: {args.output}")
    else:
        print("❌ FASTA文件生成失败")

if __name__ == "__main__":
    main()