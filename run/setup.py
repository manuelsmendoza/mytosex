import os
import sys
import json
from src.func_setup import *
from src.func_util import tnow

# Checking for the dependencies
print(tnow() + " INFO: Checking for the dependencies", file=sys.stdout)

module_name = ["yaml", "Bio"]
tool_name = ["fasterq-dump", "pigz", "bowtie2"]
detection = {}
for dependency in module_name:
    detection.update({dependency: check_dep(dependency, "module")})
for dependency in tool_name:
    detection.update({dependency: check_dep(dependency, "tool")})

if not all(list(detection.values())):
    for dependency in list(detection.keys()):
        if not detection[dependency]:
            print(tnow() + " FAIL: " + dependency + " is not installed", file=sys.stderr)
    raise ModuleNotFoundError("Any module required is no installed")

# Checking the settings
print(tnow() + " INFO: Checking the settings values", file=sys.stdout)
settings = check_settings(os.getenv("MYTOSEX_SETTINGS"))

# Check the output directory
if os.path.exists(settings["output_dir"]) and os.path.isdir(settings["output_dir"]):
    print(tnow() + " WARN: Output directory already exists", file=sys.stdout)
else:
    try:
        os.mkdir(settings["output_dir"])
    except OSError:
        print(tnow() + " FAIL: Creation of the directory %s failed " % settings["output_dir"], file=sys.stderr)
    else:
        print(tnow() + " INFO: The output will be store at %s " % settings["output_dir"], file=sys.stdout)

# Creating subdirectories
subdir_names = ["reads", "reference", "tmp", "data", "figures"]
for subdirectory in [os.path.join(settings["output_dir"], subdir) for subdir in subdir_names]:
    if not os.path.exists(subdirectory):
        try:
            os.mkdir(subdirectory)
        except OSError:
            print(tnow() + " FAIL: Creation of the directory %s failed " % subdirectory, file=sys.stderr)

# Export the values of the configuration
with open(os.path.join(settings["output_dir"], "settings.json"), "w") as settings_out:
    json.dump(settings, settings_out)

# Change the configuration path
os.environ["MYTOSEX_SETTINGS"] = os.path.join(settings["output_dir"], "settings.json")

# Create a milestone
open(os.path.join(settings["output_dir"], ".settings.ok"), "w").close()
print(tnow() + " INFO: Settings values checked correctly", file=sys.stdout)
