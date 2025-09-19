#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
iTAK 2.0 依赖检测模块
检查程序运行所需的Python包和外部工具是否已安装
"""

import sys
import subprocess
import importlib
import shutil
from pathlib import Path
import os

class DependencyChecker:
    """依赖检测器类"""
    
    def __init__(self):
        self.missing_dependencies = []
        self.missing_optional_dependencies = []
        self.warnings = []
        
        # 定义必需的Python包（核心功能必须）
        self.required_python_packages = {
            'Bio': 'Biopython生物信息学库',
            'pandas': 'Pandas数据处理库',
            'numpy': 'NumPy数值计算库',
            'json': 'JSON处理库（标准库）',
            'csv': 'CSV处理库（标准库）',
            'argparse': '命令行参数解析库（标准库）',
            'subprocess': '子进程管理库（标准库）',
            'pathlib': '路径处理库（标准库）',
            'os': '操作系统接口库（标准库）',
            'sys': '系统相关库（标准库）',
            'time': '时间处理库（标准库）',
            'datetime': '日期时间库（标准库）',
            'shutil': '文件操作库（标准库）',
            'collections': '集合类库（标准库）',
            'warnings': '警告处理库（标准库）'
        }
        
        # 定义可选的Python包（特定功能需要）
        self.optional_python_packages = {
            'torch': 'PyTorch深度学习框架（仅预测功能需要）'
        }
        
        # 定义所需的外部工具
        self.required_external_tools = {
            'hmmscan': 'HMMER包中的hmmscan工具，用于HMM搜索',
            'java': 'Java运行环境，InterProScan需要',
            'perl': 'Perl解释器，InterProScan需要',
            'python3': 'Python3解释器'
        }
        
        # 定义关键文件路径（模块在module目录中，需要向上一级找到项目根目录）
        self.script_dir = Path(__file__).parent.parent.absolute()
        self.required_files = {
            'interproscan.sh': self.script_dir / "db" / "interproscan" / "interproscan.sh",
            'self_build.hmm': self.script_dir / "db" / "self_build_hmm" / "self_build.hmm",
            'predict.py': self.script_dir / "pre_model" / "predict.py",
            'model.pth': self.script_dir / "pre_model" / "model.pth"
        }
    
    def check_python_package(self, package_name):

        try:
            importlib.import_module(package_name)
            return True
        except ImportError:
            return False
    
    def check_external_tool(self, tool_name):

        return shutil.which(tool_name) is not None
    
    def check_file_exists(self, file_path):

        return file_path.exists()
    
    def check_interproscan_setup(self):

        interproscan_script = self.required_files['interproscan.sh']
        
        if not interproscan_script.exists():
            return False
        
        # 检查InterProScan脚本是否可执行
        if not os.access(interproscan_script, os.X_OK):
            self.warnings.append(f"InterProScan脚本不可执行: {interproscan_script}")
            return False
        
        # 检查InterProScan数据目录
        interproscan_dir = interproscan_script.parent
        data_dir = interproscan_dir / "data"
        
        if not data_dir.exists():
            self.warnings.append(f"InterProScan数据目录不存在: {data_dir}")
            return False
        
        return True
    
    def check_pytorch_gpu_support(self):

        try:
            import torch
            
            # 检查PyTorch版本
            torch_version = torch.__version__
            
            if torch.cuda.is_available():
                gpu_count = torch.cuda.device_count()
                gpu_name = torch.cuda.get_device_name(0) if gpu_count > 0 else "Unknown"
                cuda_version = torch.version.cuda
                return f"GPU支持可用 - PyTorch {torch_version} (CUDA {cuda_version}), 设备数量: {gpu_count}, 主设备: {gpu_name}"
            else:
                # 检查是否是CPU版本
                if '+cpu' in torch_version:
                    return f"PyTorch {torch_version} (CPU版本) - 如需GPU加速，请安装GPU版本"
                else:
                    return f"PyTorch {torch_version} - GPU支持不可用，将使用CPU进行计算"
        except ImportError:
            return "PyTorch未安装，无法检查GPU支持"
    
    def check_prediction_dependencies(self):

        # 检查PyTorch
        if not self.check_python_package('torch'):
            return False, "预测功能需要PyTorch，但PyTorch未安装"
        
        # 检查预测模型文件
        model_file = self.required_files.get('model.pth')
        if not self.check_file_exists(model_file):
            return False, f"预测功能需要模型文件，但文件不存在: {model_file}"
        
        predict_script = self.required_files.get('predict.py')
        if not self.check_file_exists(predict_script):
            return False, f"预测功能需要预测脚本，但文件不存在: {predict_script}"
        
        return True, "预测功能依赖检查通过"
    
    def run_full_check(self):

        print("开始检查iTAK 2.0依赖...")
        print("=" * 60)
        
        all_dependencies_met = True
        
        # 检查必需的Python包
        print("\n检查必需Python包依赖:")
        for package, description in self.required_python_packages.items():
            if self.check_python_package(package):
                print(f"  [成功] {package:<15} - {description}")
            else:
                print(f"  [错误] {package:<15} - {description}")
                self.missing_dependencies.append(f"Python包: {package}")
                all_dependencies_met = False
        
        # 检查可选的Python包
        print("\n检查可选Python包依赖:")
        for package, description in self.optional_python_packages.items():
            if self.check_python_package(package):
                print(f"  [成功] {package:<15} - {description}")
            else:
                print(f"  [警告] {package:<15} - {description}")
                self.missing_optional_dependencies.append(f"Python包: {package}")
                print(f"      注意: 缺少{package}不会影响基本功能，但预测功能将不可用")
        
        # 检查外部工具
        print("\n检查外部工具:")
        for tool, description in self.required_external_tools.items():
            if self.check_external_tool(tool):
                print(f"  [成功] {tool:<15} - {description}")
            else:
                print(f"  [错误] {tool:<15} - {description}")
                self.missing_dependencies.append(f"外部工具: {tool}")
                all_dependencies_met = False
        
        # 检查关键文件
        print("\n检查关键文件:")
        for file_name, file_path in self.required_files.items():
            if self.check_file_exists(file_path):
                print(f"  [成功] {file_name:<20} - {file_path}")
            else:
                print(f"  [错误] {file_name:<20} - {file_path}")
                self.missing_dependencies.append(f"文件: {file_path}")
                all_dependencies_met = False
        
        # 检查InterProScan配置
        print("\n检查InterProScan配置:")
        if self.check_interproscan_setup():
            print("  [成功] InterProScan配置正常")
        else:
            print("  [错误] InterProScan配置有问题")
            self.missing_dependencies.append("InterProScan配置")
            all_dependencies_met = False
        
        # 检查PyTorch GPU支持（仅在PyTorch可用时）
        print("\n检查GPU支持（预测功能相关）:")
        if self.check_python_package('torch'):
            gpu_status = self.check_pytorch_gpu_support()
            print(f"  [信息] {gpu_status}")
        else:
            print("  [警告] PyTorch未安装，无法检查GPU支持状态")
        
        # 显示警告信息
        if self.warnings:
            print("\n[警告] 警告信息:")
            for warning in self.warnings:
                print(f"  [警告] {warning}")
        
        # 总结
        print("\n" + "=" * 60)
        if all_dependencies_met:
            print("[成功] 所有必需依赖检查通过！iTAK 2.0可以正常运行。")
            if self.missing_optional_dependencies:
                print("[警告] 注意：以下可选依赖缺失，某些功能可能不可用:")
                for dep in self.missing_optional_dependencies:
                    print(f"  - {dep}")
                print("  如需使用预测功能，请安装PyTorch。")
        else:
            print("[错误] 发现缺失的必需依赖！请安装以下组件:")
            for dep in self.missing_dependencies:
                print(f"  - {dep}")
            print("\n安装建议:")
            self._print_installation_suggestions()
        
        return all_dependencies_met
    
    def _print_installation_suggestions(self):
        """打印安装建议"""
        print("\n安装建议:")
        
        # Python包安装建议
        python_packages_missing = [dep.split(': ')[1] for dep in self.missing_dependencies if dep.startswith('Python包')]
        if python_packages_missing:
            print("\n  Python包安装:")
            for pkg in python_packages_missing:
                if pkg == 'torch':
                    print("    PyTorch安装 (选择适合您系统的版本):")
                    print("      CPU版本: pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu")
                    print("      GPU版本 (CUDA): pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118")
                    print("      更多版本请访问: https://pytorch.org/get-started/locally/")
                elif pkg == 'Bio':
                    print(f"    pip install biopython")
                elif pkg in ['pandas', 'numpy']:
                    print(f"    pip install {pkg}")
        
        # 外部工具安装建议
        external_tools_missing = [dep.split(': ')[1] for dep in self.missing_dependencies if dep.startswith('外部工具')]
        if external_tools_missing:
            print("\n  外部工具安装:")
            for tool in external_tools_missing:
                if tool == 'hmmscan':
                    print("    安装HMMER: http://hmmer.org/download.html")
                    print("    Ubuntu/Debian: sudo apt-get install hmmer")
                    print("    CentOS/RHEL: sudo yum install hmmer")
                elif tool == 'java':
                    print("    安装Java: https://www.oracle.com/java/technologies/downloads/")
                    print("    Ubuntu/Debian: sudo apt-get install openjdk-11-jdk")
                    print("    CentOS/RHEL: sudo yum install java-11-openjdk")
                elif tool == 'perl':
                    print("    安装Perl: https://www.perl.org/get.html")
                    print("    Ubuntu/Debian: sudo apt-get install perl")
                    print("    CentOS/RHEL: sudo yum install perl")

def main():
    """主函数"""
    checker = DependencyChecker()
    success = checker.run_full_check()
    
    if not success:
        sys.exit(1)
    
    return success

if __name__ == "__main__":
    main()