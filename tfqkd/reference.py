"""Load only unmodified author noise functions; never execute notebook cells."""

import ast
import json
import warnings
import numpy as np

from .config import ROOT


def author_functions():
    notebook = json.loads((ROOT / "sources/data/QKD.ipynb").read_text())
    source = "".join(notebook["cells"][3]["source"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SyntaxWarning)
        parsed = ast.parse(source)
    functions = [node for node in parsed.body if isinstance(node, ast.FunctionDef)]
    module = ast.Module(body=functions, type_ignores=[])
    namespace = {"np": np}
    exec(compile(module, "QKD.ipynb:cell3:unmodified-functions", "exec"), namespace)
    return namespace


if __name__ == "__main__":
    print("Loaded author functions:", ", ".join(k for k, v in author_functions().items() if callable(v)))
