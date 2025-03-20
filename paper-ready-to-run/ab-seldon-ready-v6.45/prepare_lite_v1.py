
import os
import subprocess
import time
import re
import glob
import torch
from contextlib import redirect_stdout
from itertools import count
from ImmuneBuilder import ABodyBuilder2
from ImmuneBuilder.refine import refine
predictor = ABodyBuilder2(numbering_scheme="martin")

def gen_script(steps, vers):
    step_list = steps.replace("\n","").replace("\r","").split("|")

    cmd5 = "rm optimize.sh"
    subprocess.run(cmd5, shell=True, capture_output=True)

    opt_script = open("optimize.sh", "a", newline="")
    with redirect_stdout(opt_script):
        for stp in step_list:
            if stp == "h3":
                print("python3 swap_h3_v" + vers + ".py")
            if stp == "rep":
                print("python3 swap_rep_v" + vers + ".py")
            if stp == "bnk":
                print("python3 swap_bnk_v" + vers + ".py")
            if stp == "mut":
                print("python3 mutagenesis_v" + vers + ".py")
            if stp == "fwk":
                print("python3 memory_frame_v" + vers + ".py")

##########
# Reading configuration file
input_name = "cfg:prepare|input_name"
num_struc = "cfg:prepare|num_struc"
code_name = "cfg:prepare|code_name"
vers = "cfg:prepare|seldon_version"
steps = "cfg:steps"
scoring_method = "cfg:scoring_method"
ph = "7"
pymol_command = "pymol"

cfg = open("swap_settings.cfg").readlines()

for setting in cfg:
    if "steps" in setting:
        steps = str(setting.split("=")[1].replace("\n", ""))
    if "minimization_pH" in setting:
        ph = str(setting.split("=")[1].replace("\n", ""))
    if "pymol_command" in setting:
        pymol_command = str(setting.split("=")[1].replace("\n", ""))
    if "scoring_method" in setting:
        scoring_method = str(setting.split("=")[1].replace("\n", ""))
    if "prepare|number_of_structures" in setting:
        num_struc = str(setting.split("=")[1].replace("\n", ""))
    if "prepare|code_name" in setting:
        code_name = str(setting.split("=")[1].replace("\n", ""))
    if "prepare|input_name" in setting:
        input_name = str(setting.split("=")[1].replace("\n", "")).replace(".pdb","").replace(".fasta","").replace(".anarci","")
    if "prepare|seldon_version" in setting:
        vers = str(setting.split("=")[1].replace("\n", ""))


def lite_prep():
    last_log = []
    for log in glob.glob("swap_00" + str(amount) +".log"):
        open_log = open(log).readlines()
        last_log.append(open_log[-1].replace("\n", "").replace("\r", "").replace("FINAL","bs1"))

    save_log = open("swap_00" + str(amount) +".log","w",newline="")
    with redirect_stdout(save_log):
        print(last_log[-1])

    cmd3 = "mkdir prev_step_files/"
    subprocess.run(cmd3, shell=True, capture_output=True)

    for final_file in glob.glob("*FINAL*"):
        cmd2 = "cp " + final_file + " " + final_file.replace("FINAL","bs1")
        subprocess.run(cmd2, shell=True, capture_output=True)
        cmd4 = "mv " + final_file + " prev_step_files/"
        subprocess.run(cmd4, shell=True, capture_output=True)


print("Starting pipeline: Input preparation - CONTINUING FROM PREVIOUS STEP (LITE)")

for amount in range(1,(int(num_struc) + 1)):
    lite_prep()
    gen_script(steps, vers)

print("Input preparation completed - CONTINUING FROM PREVIOUS STEP (LITE)")
