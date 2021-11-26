"""
    Add some information here
"""
from src.func_setup import *
from src.func_util import tnow

# Checking for the dependencies
print(tnow() + " INFO: Checking for the dependencies")

mod_name = ["yaml"]
tool_name = ["bowtie2"]
dep_res = {}
for dep in mod_name:
    dep_res.update({dep: check_dep(dep_name=dep, dep_type="module")})
for dep in tool_name:
    dep_res.update({dep: check_dep(dep_name=dep, dep_type="tool")})

if not all(list(dep_res.values())):
    for dep in list(dep_res.keys()):
        if not dep_res[dep]:
            print(tnow() + " FAIL: " + dep + " is not installed")
    raise ModuleNotFoundError
