from __future__ import annotations

# helper：动态加载上一级目录的 00_config.py。
# 这样主脚本既可以直接运行，也不需要把 scripts_format_python 安装成 Python package。
import importlib.util
from pathlib import Path


def load_config():
    # importlib 按文件路径加载配置模块，避免和系统中同名模块冲突。
    config_path = Path(__file__).resolve().parents[1] / "00_config.py"
    spec = importlib.util.spec_from_file_location("scripts_format_python_config", config_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load config from {config_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
