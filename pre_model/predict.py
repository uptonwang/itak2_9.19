#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
One-Hot编码多尺度CNN模型蛋白质预测脚本
- 基于train_onehot_cnn.py的多尺度OneHotLightweightCNN模型
- 支持单个FASTA文件的转录因子预测
- 输出蛋白质预测概率表
- 支持自定义置信率阈值
- 适配train_data8.28的多尺度CNN架构
"""

import os
import argparse

# 解决MKL库冲突问题 - 必须在导入torch之前设置
os.environ['MKL_SERVICE_FORCE_INTEL'] = '1'
os.environ['MKL_THREADING_LAYER'] = 'GNU'
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import numpy as np
from Bio import SeqIO
import csv
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')



# 氨基酸到数字的映射
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
aa_to_idx = {aa: idx for idx, aa in enumerate(AMINO_ACIDS)}
aa_to_idx['X'] = len(AMINO_ACIDS)  # 未知氨基酸

class OneHotLightweightCNN(nn.Module):
    """使用One-Hot编码的多尺度CNN模型
    
    使用三种不同配置的多层卷积核组合，通过并行结构实现多尺度特征提取：
    配置1: (4×21, 4×1, 16×1)
    配置2: (12×21, 8×1, 4×1) 
    配置3: (16×21, 4×1, 4×1)
    """
    def __init__(self, vocab_size=21, num_classes=2, max_length=1000):
        super(OneHotLightweightCNN, self).__init__()
        
        # 输入: (batch_size, max_length, vocab_size)
        # 需要转置为: (batch_size, vocab_size, max_length) 用于Conv1d
        
        # 配置1: (4×21, 4×1, 16×1) - 小尺度精细特征
        self.branch1_conv1 = nn.Conv1d(vocab_size, 64, kernel_size=4, padding=2)  # 4×21等效
        self.branch1_conv2 = nn.Conv1d(64, 128, kernel_size=4, padding=2)        # 4×1
        self.branch1_conv3 = nn.Conv1d(128, 256, kernel_size=16, padding=8)      # 16×1
        
        # 配置2: (12×21, 8×1, 4×1) - 中等尺度特征
        self.branch2_conv1 = nn.Conv1d(vocab_size, 64, kernel_size=12, padding=6) # 12×21等效
        self.branch2_conv2 = nn.Conv1d(64, 128, kernel_size=8, padding=4)        # 8×1
        self.branch2_conv3 = nn.Conv1d(128, 256, kernel_size=4, padding=2)       # 4×1
        
        # 配置3: (16×21, 4×1, 4×1) - 大尺度全局特征
        self.branch3_conv1 = nn.Conv1d(vocab_size, 64, kernel_size=16, padding=8) # 16×21等效
        self.branch3_conv2 = nn.Conv1d(64, 128, kernel_size=4, padding=2)        # 4×1
        self.branch3_conv3 = nn.Conv1d(128, 256, kernel_size=4, padding=2)       # 4×1
        
        # 批归一化层 - 分支1
        self.bn1_1 = nn.BatchNorm1d(64)
        self.bn1_2 = nn.BatchNorm1d(128)
        self.bn1_3 = nn.BatchNorm1d(256)
        
        # 批归一化层 - 分支2
        self.bn2_1 = nn.BatchNorm1d(64)
        self.bn2_2 = nn.BatchNorm1d(128)
        self.bn2_3 = nn.BatchNorm1d(256)
        
        # 批归一化层 - 分支3
        self.bn3_1 = nn.BatchNorm1d(64)
        self.bn3_2 = nn.BatchNorm1d(128)
        self.bn3_3 = nn.BatchNorm1d(256)
        
        # 特征融合层
        self.fusion_conv = nn.Conv1d(768, 512, kernel_size=3, padding=1)  # 256*3=768
        self.fusion_bn = nn.BatchNorm1d(512)
        
        # 深层特征提取
        self.deep_conv1 = nn.Conv1d(512, 256, kernel_size=3, padding=1)
        self.deep_conv2 = nn.Conv1d(256, 128, kernel_size=3, padding=1)
        self.deep_bn1 = nn.BatchNorm1d(256)
        self.deep_bn2 = nn.BatchNorm1d(128)
        
        # 注意力机制
        self.attention = nn.MultiheadAttention(embed_dim=128, num_heads=8, batch_first=True)
        
        # 池化层
        self.global_avg_pool = nn.AdaptiveAvgPool1d(1)
        self.global_max_pool = nn.AdaptiveMaxPool1d(1)
        
        # Dropout层（增强正则化）
        self.dropout1 = nn.Dropout(0.2)
        self.dropout2 = nn.Dropout(0.3)
        self.dropout3 = nn.Dropout(0.4)
        self.dropout4 = nn.Dropout(0.5)
        
        # 分类层（多尺度特征分类）
        self.classifier = nn.Sequential(
            nn.Linear(256, 512),  # 128*2=256 (avg+max pooling)
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(256, 64),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(64, num_classes)
        )
    
    def forward(self, x):
        # x: (batch_size, max_length, vocab_size)
        x = x.transpose(1, 2)  # (batch_size, vocab_size, max_length)
        
        # 分支1: (4×21, 4×1, 16×1) - 小尺度精细特征
        branch1 = torch.relu(self.bn1_1(self.branch1_conv1(x)))
        branch1 = self.dropout1(branch1)
        branch1 = torch.relu(self.bn1_2(self.branch1_conv2(branch1)))
        branch1 = self.dropout2(branch1)
        branch1 = torch.relu(self.bn1_3(self.branch1_conv3(branch1)))
        branch1 = self.dropout1(branch1)
        
        # 分支2: (12×21, 8×1, 4×1) - 中等尺度特征
        branch2 = torch.relu(self.bn2_1(self.branch2_conv1(x)))
        branch2 = self.dropout1(branch2)
        branch2 = torch.relu(self.bn2_2(self.branch2_conv2(branch2)))
        branch2 = self.dropout2(branch2)
        branch2 = torch.relu(self.bn2_3(self.branch2_conv3(branch2)))
        branch2 = self.dropout1(branch2)
        
        # 分支3: (16×21, 4×1, 4×1) - 大尺度全局特征
        branch3 = torch.relu(self.bn3_1(self.branch3_conv1(x)))
        branch3 = self.dropout1(branch3)
        branch3 = torch.relu(self.bn3_2(self.branch3_conv2(branch3)))
        branch3 = self.dropout2(branch3)
        branch3 = torch.relu(self.bn3_3(self.branch3_conv3(branch3)))
        branch3 = self.dropout1(branch3)
        
        # 多尺度特征融合
        x = torch.cat([branch1, branch2, branch3], dim=1)  # (batch_size, 768, max_length)
        x = torch.relu(self.fusion_bn(self.fusion_conv(x)))
        x = self.dropout3(x)
        
        # 深层特征提取
        x = torch.relu(self.deep_bn1(self.deep_conv1(x)))
        x = self.dropout2(x)
        
        x = torch.relu(self.deep_bn2(self.deep_conv2(x)))
        x = self.dropout1(x)
        
        # 注意力机制
        x = x.transpose(1, 2)  # (batch_size, max_length, 128)
        x_att, _ = self.attention(x, x, x)
        x = x + x_att  # 残差连接
        x = x.transpose(1, 2)  # (batch_size, 128, max_length)
        
        # 双重池化
        x_avg = self.global_avg_pool(x)  # (batch_size, 128, 1)
        x_max = self.global_max_pool(x)  # (batch_size, 128, 1)
        x = torch.cat([x_avg, x_max], dim=1)  # (batch_size, 256, 1)
        x = x.squeeze(-1)  # (batch_size, 256)
        
        # 分类
        x = self.classifier(x)
        return x

def parse_fasta(fasta_file):
    """解析FASTA文件"""
    sequences = []
    current_header = None
    current_sequence = ""
    
    with open(fasta_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences.append({
                        'header': current_header,
                        'sequence': current_sequence
                    })
                current_header = line[1:]
                current_sequence = ""
            else:
                current_sequence += line
        
        if current_header is not None:
            sequences.append({
                'header': current_header,
                'sequence': current_sequence
            })
    
    return sequences



def sequence_to_onehot(sequence, max_length=1000, vocab_size=21):
    """将氨基酸序列转换为one-hot编码"""
    # 初始化one-hot矩阵
    onehot = np.zeros((max_length, vocab_size), dtype=np.float32)
    
    # 编码序列
    for i, aa in enumerate(sequence[:max_length]):
        aa_idx = aa_to_idx.get(aa, aa_to_idx['X'])
        onehot[i, aa_idx] = 1.0
    
    # 对于短于max_length的序列，剩余位置保持全零（padding）
    return onehot

def predict_sequences_batch(model, sequences, device, max_length=1000, batch_size=16):
    """批量预测序列"""
    model.eval()
    all_predictions = []
    
    with torch.no_grad():
        for i in range(0, len(sequences), batch_size):
            batch_sequences = sequences[i:i+batch_size]
            
            # 编码序列为One-Hot
            encoded_batch = []
            for seq_data in batch_sequences:
                onehot_encoded = sequence_to_onehot(seq_data['sequence'], max_length)
                encoded_batch.append(onehot_encoded)
            
            # 转换为tensor
            batch_tensor = torch.tensor(np.stack(encoded_batch), dtype=torch.float).to(device)
            
            # 预测
            outputs = model(batch_tensor)
            probabilities = torch.softmax(outputs, dim=1)
            
            # 提取概率
            for j, seq_data in enumerate(batch_sequences):
                non_tf_prob = probabilities[j][0].item()
                tf_prob = probabilities[j][1].item()
                
                all_predictions.append({
                    'header': seq_data['header'],
                    'tf_probability': tf_prob,
                    'non_tf_probability': non_tf_prob
                })
    
    return all_predictions



def predict_fasta(fasta_file, model, device, threshold=0.1, max_length=1000, batch_size=16):
    """预测FASTA文件中的蛋白质序列"""
    print(f"\n处理文件: {fasta_file}")
    
    # 解析FASTA文件
    sequences = parse_fasta(fasta_file)
    print(f"序列数量: {len(sequences)}")
    
    if not sequences:
        print(f"警告: 没有找到有效序列")
        return [], []
    
    # 预测
    predictions = predict_sequences_batch(model, sequences, device, max_length, batch_size)
    
    # 应用阈值分类
    results = []
    tf_headers = []  # 存储符合阈值的蛋白质名称
    tf_count = 0
    
    for pred in predictions:
        tf_prob = pred['tf_probability']
        predicted_class = 'TF' if tf_prob >= threshold else 'Non-TF'
        confidence = max(tf_prob, pred['non_tf_probability'])
        
        result = {
            'header': pred['header'],
            'predicted_class': predicted_class,
            'tf_probability': tf_prob,
            'non_tf_probability': pred['non_tf_probability'],
            'confidence': confidence
        }
        
        results.append(result)
        
        if predicted_class == 'TF':
            tf_count += 1
            tf_headers.append(pred['header'])  # 添加到TF名称列表
    
    print(f"预测结果: {tf_count} 个TF, {len(predictions) - tf_count} 个Non-TF")
    print(f"TF比例: {tf_count/len(predictions)*100:.2f}%")
    
    return results, tf_headers

def save_predictions(predictions, output_file):
    """保存预测结果到CSV文件"""
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Header', 'Predicted_Class', 'TF_Probability', 'Non_TF_Probability', 'Confidence'])
        for pred in predictions:
            writer.writerow([
                pred['header'],
                pred['predicted_class'],
                f"{pred['tf_probability']:.4f}",
                f"{pred['non_tf_probability']:.4f}",
                f"{pred['confidence']:.4f}"
            ])
    print(f"预测结果已保存到: {output_file}")

def save_tf_predictions(predictions, output_file):
    """保存TF预测结果到CSV文件"""
    tf_predictions = [p for p in predictions if p['predicted_class'] == 'TF']
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Header', 'TF_Probability', 'Confidence'])
        for pred in tf_predictions:
            writer.writerow([
                pred['header'],
                f"{pred['tf_probability']:.4f}",
                f"{pred['confidence']:.4f}"
            ])
    print(f"TF预测结果已保存到: {output_file}")

def load_model(model_path, device, max_length=1000):
    """加载训练好的多尺度CNN模型"""
    print(f"加载多尺度CNN模型: {model_path}")
    
    # 加载检查点
    checkpoint = torch.load(model_path, map_location=device)
    
    # 获取模型配置
    if 'model_config' in checkpoint:
        config = checkpoint['model_config']
        vocab_size = config.get('vocab_size', 21)
        max_length = config.get('max_length', 1000)
        num_classes = config.get('num_classes', 2)
    else:
        # 默认配置
        vocab_size = 21
        num_classes = 2
    
    # 创建多尺度CNN模型
    model = OneHotLightweightCNN(
        vocab_size=vocab_size,
        num_classes=num_classes,
        max_length=max_length
    ).to(device)
    
    # 加载权重
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    
    # 计算模型参数
    total_params = sum(p.numel() for p in model.parameters())
    
    print(f"多尺度CNN模型加载成功:")
    print(f"  - 模型架构: 三分支并行多尺度CNN")
    print(f"  - 总参数量: {total_params:,}")
    print(f"  - 词汇表大小: {vocab_size}")
    print(f"  - 最大序列长度: {max_length}")
    print(f"  - 类别数: {num_classes}")
    if 'val_acc' in checkpoint:
        print(f"  - 验证准确率: {checkpoint['val_acc']:.4f}")
    if 'epoch' in checkpoint:
        print(f"  - 训练轮数: {checkpoint['epoch'] + 1}")
    
    return model, max_length

def main():
    # 获取脚本所在目录的绝对路径
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_model_path = os.path.join(script_dir, 'model.pth')
    
    parser = argparse.ArgumentParser(description='多尺度CNN模型蛋白质预测')
    parser.add_argument('--fasta', type=str, required=True,
                       help='输入FASTA文件路径')
    parser.add_argument('--model', type=str, 
                        default=default_model_path,
                        help=f'模型文件路径 (默认: {default_model_path})')
    parser.add_argument('--threshold', type=float, default=0.1,
                       help='TF预测阈值 (默认: 0.1, 即10%%)')
    parser.add_argument('--output', type=str, default=None,
                       help='输出CSV文件路径 (默认: 基于输入文件名自动生成)')
    parser.add_argument('--batch_size', type=int, default=16,
                       help='批次大小 (默认: 16)')
    parser.add_argument('--debug', action='store_true',
                       help='调试模式，生成CSV文件 (默认: False)')
    parser.add_argument('--output-tf-list', action='store_true',
                       help='输出TF列表到标准输出，用于内存传递 (默认: False)')
    
    args = parser.parse_args()
    
    # 检查输入文件
    if not os.path.exists(args.fasta):
        print(f"错误: FASTA文件不存在: {args.fasta}")
        return
    
    # 检查模型文件
    if not os.path.exists(args.model):
        print(f"错误: 模型文件不存在: {args.model}")
        return
    
    # 生成输出文件名
    if args.output is None:
        base_name = os.path.splitext(os.path.basename(args.fasta))[0]
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.output = f"{base_name}_predictions_{timestamp}.csv"
    
    # 设备配置
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"使用设备: {device}")
    if torch.cuda.is_available():
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"GPU内存: {torch.cuda.get_device_properties(0).total_memory / 1024**3:.1f} GB")
    
    print(f"\n参数设置:")
    print(f"输入文件: {args.fasta}")
    print(f"模型文件: {args.model}")
    print(f"预测阈值: {args.threshold} ({args.threshold*100:.1f}%)")
    print(f"输出文件: {args.output}")
    print(f"批次大小: {args.batch_size}")
    
    # 加载模型
    model, max_length = load_model(args.model, device)
    
    # 预测
    predictions, tf_headers = predict_fasta(
        fasta_file=args.fasta,
        model=model,
        device=device,
        threshold=args.threshold,
        max_length=max_length,
        batch_size=args.batch_size
    )
    
    if predictions:
        # 仅在debug模式下保存CSV文件
        if args.debug:
            # 保存所有预测结果
            save_predictions(predictions, args.output)
            
            # 保存TF预测结果
            tf_output = args.output.replace('.csv', '_tf_only.csv')
            save_tf_predictions(predictions, tf_output)
        else:
            print("非debug模式，跳过CSV文件生成")
        
        # 输出TF列表到标准输出（用于内存传递）
        if args.output_tf_list:
            tf_headers = [p['header'] for p in predictions if p['predicted_class'] == 'TF']
            print("TF_LIST_START")
            for header in tf_headers:
                print(f"TF_HEADER:{header}")
            print("TF_LIST_END")
        
        print(f"\n=== 预测完成 ===")
        print(f"总序列数: {len(predictions)}")
        tf_count = sum(1 for p in predictions if p['predicted_class'] == 'TF')
        print(f"预测TF数: {tf_count}")
        print(f"TF比例: {tf_count/len(predictions)*100:.2f}%")
        print(f"使用阈值: {args.threshold} ({args.threshold*100:.1f}%)")
    else:
        print("错误: 预测失败或没有有效序列")

if __name__ == '__main__':
    main()