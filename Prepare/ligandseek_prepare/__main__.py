"""Prepare 独立命令. 将一个原始样本转换为 Infer 可以直接读取的 sample.json."""

import argparse
import os
from pathlib import Path

from .pipeline import prepare_sample


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path, help="单样本原始输入 JSON")
    parser.add_argument("--models", type=Path, default=os.environ.get("LIGANDSEEK_MODELS"), required="LIGANDSEEK_MODELS" not in os.environ)
    parser.add_argument("--output", type=Path, required=True, help="本次准备产物的新目录")
    args = parser.parse_args()
    print(prepare_sample(args.raw, args.models, args.output))
