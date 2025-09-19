================================================================================
                              iTAK2 - 转录因子预测和分析工具
================================================================================

项目简介
--------
iTAK2 是一个用于转录因子预测和分析的生物信息学工具。该工具结合了深度学习模型预测
和传统的序列分析方法，能够准确识别和分类蛋白质序列中的转录因子。

主要功能
--------
1. 转录因子预测：使用深度学习模型预测蛋白质序列是否为潜在TF/TR序列
2. 序列分析：通过InterProScan和hmmscan进行结构域分析
3. 转录因子分类：基于规则引擎对转录因子进行家族分类
4. 结果输出：生成详细的分析报告和分类序列文件

系统要求
--------
- Python 3.7+
- Java 8+ (InterProScan需要)
- HMMER 3.0+ (hmmscan工具)

依赖包
------
必需的Python包：
- Biopython (Bio)
- pandas
- numpy
- 标准库：json, csv, argparse, subprocess, pathlib, os, sys, time, datetime, shutil

可选的Python包：
- PyTorch (仅预测功能需要)

外部工具：
- hmmscan (HMMER包)
- java (InterProScan需要)

项目结构
--------
itak2/
├── itak2.0.py              # 主程序脚本
├── module/                 # 功能模块目录
│   ├── check_dependencies.py    # 依赖检测模块
│   ├── validate_fasta.py        # FASTA文件验证模块
│   ├── classification.py        # 转录因子分类模块
│   ├── get_fasta.py             # 序列提取模块
│   ├── get_json.py              # InterProScan结果处理模块
│   ├── get_rule.py              # 规则文件解析模块
│   └── selfbuild_hmm.py         # hmmscan结果处理模块
├── pre_model/              # 预测模型目录
│   ├── model.pth               # 深度学习模型文件
│   └── predict.py              # 预测脚本
├── db/                     # 数据库目录
│   ├── interproscan/           # InterProScan程序和数据库
│   ├── self_build_hmm/         # 自建HMM数据库
│   └── rule.txt               # 转录因子分类规则文件
├── output/                 # 输出目录
├── temp/                   # 临时文件目录
└── test_protein.fasta      # 示例序列文件

安装说明
--------
1. 确保Python 3.7+已安装
2. 安装必需的Python包：
   pip install biopython pandas numpy

3. 安装可选的PyTorch（如需使用预测功能）：
   # CPU版本
   pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
   
   # GPU版本 (CUDA 11.8)
   pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

4. 安装HMMER：
   # Ubuntu/Debian
   sudo apt-get install hmmer
   
   # CentOS/RHEL
   sudo yum install hmmer
   
   # macOS
   brew install hmmer

5. 确保Java已安装（InterProScan需要）：
   java -version

使用方法
--------

基本语法：
python itak2.0.py [选项] -i <输入文件>

主要参数：
-i, --input          输入FASTA文件路径 (必需)
-o, --output         输出目录路径 (可选，不选将输出到output文件下的项目文件夹内)
-t, --threshold      预测阈值，范围0-1 (默认: 0.1)
--predict           启用预测功能
--debug             启用调试模式
--score             InterProScan结果的score阈值 (默认: 1.0)

功能选项：
--predict           使用深度学习模型进行转录因子预测
--extract-sequences 提取预测的TF序列 (仅在--predict时有效，无需选择，仅用作测试)
--test-mode         测试模式，使用指定文件进行分类验证
--appl              InterProScan应用程序列表 (默认: 全部)

依赖检查：
--check-deps        仅检查依赖并退出
--skip-deps-check   跳过依赖检查 (不推荐)

使用示例
--------

1. 基本分析（直接分析输入序列）：
   python itak2.0.py -i test_protein.fasta

2. 使用预测功能：
   python itak2.0.py --predict -i test_protein.fasta -t 0.1

3. 指定输出目录：
   python itak2.0.py --predict -i test_protein.fasta -o /path/to/output

4. 启用调试模式：
   python itak2.0.py --predict -i test_protein.fasta --debug

5. 自定义InterProScan数据库：
   python itak2.0.py -i test_protein.fasta --appl CDD,Pfam,SMART

6. 测试模式：
   python itak2.0.py --test-mode -i input.fasta -json ipr_result.json -spechmm hmmscan_result.tbl

7. 检查依赖：
   python itak2.0.py --check-deps

输出文件
--------
程序会在输出目录中创建以下文件和目录：

result/                     # 主要结果目录
├── match_tbl.txt              # 转录因子分类结果表格
├── match.json                 # 转录因子分类结果JSON (调试模式)
├── *_tf_classified.fasta      # 分类后的转录因子序列
├── processed_ipr_domains.json # 处理后的InterProScan结果 (调试模式)
└── pfamspec.json             # hmmscan结果 (调试模式)

fasta/                      # 序列文件目录 (预测模式)
├── *_tf_sequences.fasta       # 预测的转录因子序列
└── prediction_results.csv     # 预测结果详细表格

ipr/                        # InterProScan结果目录
└── *.json                     # InterProScan原始结果

hmmscan/                    # hmmscan结果目录
└── result.tbl                 # hmmscan原始结果

getrule.json               # 解析后的分类规则

工作流程
--------
1. 输入验证：检查FASTA文件格式和蛋白质序列
2. 依赖检查：验证所需工具和库是否可用
3. 预测阶段（可选）：使用深度学习模型预测转录因子
4. 序列分析：运行InterProScan和hmmscan进行结构域分析
5. 结果处理：解析分析结果并过滤
6. 转录因子分类：基于规则引擎进行家族分类
7. 结果输出：生成分类报告和序列文件

故障排除
--------

常见问题：

1. "依赖检查失败"
   - 检查Python包是否正确安装
   - 确认外部工具在PATH中可用
   - 使用 --check-deps 详细检查依赖状态

2. "FASTA文件验证失败"
   - 确认文件格式正确
   - 检查是否为蛋白质序列
   - 移除序列中的星号(*)字符

3. "InterProScan运行失败"
   - 确认Java已正确安装
   - 检查InterProScan数据库完整性
   - 确保有足够的磁盘空间

4. "预测功能不可用"
   - 安装PyTorch: pip install torch
   - 确认model.pth文件存在
   - 检查predict.py脚本完整性