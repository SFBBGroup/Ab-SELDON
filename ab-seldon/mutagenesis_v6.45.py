
import glob
import os
import subprocess
import random
import math
import ast
import time
import re
import torch
from os import path
from contextlib import redirect_stdout
from itertools import count
from ImmuneBuilder import ABodyBuilder2
from ImmuneBuilder.refine import refine
predictor = ABodyBuilder2(numbering_scheme="martin")


class ClearCache:
    def __enter__(self):
        torch.cuda.empty_cache()

    def __exit__(self, exc_type, exc_val, exc_tb):
        torch.cuda.empty_cache()


def model_mutated():
    # Models the modified antibody structure and checks if the CDR conformation was preserved
    sequences = {}

    if chain == "H":
        sequences["H"] = newseq
        sequences["L"] = "".join([str(item) for item in [b for b in other_chain if not b == "-"]])
    if chain == "L":
        sequences["H"] = "".join([str(item) for item in [b for b in other_chain if not b == "-"]])
        sequences["L"] = newseq

    output_file = "MUTATION_" + mutation_id + "_" + st + "_" + file.split("_")[2] + "_" + file.split("_")[-1].replace(".anarci", ".pdb")

    with ClearCache():
        antibody = predictor.predict(sequences)
        antibody.save(output_file, check_for_strained_bonds=False)

    make_fasta = output_file.replace(".pdb", ".fasta")
    if chain == "H":
        with redirect_stdout(open(make_fasta, "w", newline="")):
            print(">" + main_header[0])
            print(newseq)
            print(">" + main_header[-1])
            print("".join([str(item) for item in [b for b in other_chain if not b == "-"]]))
    if chain == "L":
        with redirect_stdout(open(make_fasta, "w", newline="")):
            print(">" + main_header[0])
            print("".join([str(item) for item in [b for b in other_chain if not b == "-"]]))
            print(">" + main_header[-1])
            print(newseq)


def conf_check(latest_ab):
    # Checks if the new CDR has kept the same pre-mutation conformation
    if check_conformation == "yes":
        print("Checking if new CDR conformation is correct...")
        first_res = str(globals()["CDR" + chain + x + "_res"][0])
        last_res = str(globals()["CDR" + chain + x + "_res"][-1])

        last_appr_ab = file.replace(".anarci", ".pdb")

        with redirect_stdout(open("pymol_conf_check.pml", "w", newline="")):
            print("load " + latest_ab)
            print("alter (chain " + chain + "), chain=\"Z\"")
            print("load " + last_appr_ab)
            print("align " + latest_ab.replace(".pdb", "") + ", " + last_appr_ab.replace(".pdb", ""))
            print("select new_cdr, chain Z and resi " + first_res + "-" + last_res + " and name CA")
            print("select old_cdr, chain " + chain + " and resi " + first_res + "-" + last_res + " and name CA")
            print("rms_cur new, old, matchmaker=4")

        run_pml_cmd = str(pymol_command + " -qc pymol_conf_check.pml")
        out = subprocess.run(run_pml_cmd, shell=True, capture_output=True)
        log_pml = out.stdout
        rmsd = float(str(log_pml).split("Executive: RMSD =")[-1].split("(")[0].replace(" ", ""))

        if rmsd <= conf_rmsd_limit:
            print("CDR has correct conformation. RMSD = " + str(rmsd))

        if rmsd > conf_rmsd_limit:
            rejected.append(tested_canon[-1])
            curr_log_rej = "rej_mod_" + curr_struct + ".log"
            updt_log_rej = open(curr_log_rej, "a", newline="")
            with redirect_stdout(updt_log_rej):
                print("\trejmut_" + tested_canon[-1] + "|" + "Wrong CDR conf. RMSD=" + str(rmsd) + "|FAIL")
            del tested_canon[-1]
            print("Wrong CDR conformation. RMSD = " + str(rmsd))

            raise Exception("Wrong CDR conformation")


def ter_fix():
    # Fixing possible incorrect chain terminations introduced by pdb4amber
    for ter in glob.glob("out_pdb4amber-addH-preterfix.pdb"):
        with open(ter, "r") as f_in, open("out_pdb4amber-addH.pdb", "w") as f_out:
            lines = f_in.readlines()
            corrected_lines = []
            for line in lines:
                if not line.startswith("TER"):
                    corrected_lines.append(line)
                if line.startswith("TER"):
                    if ("ATOM" in lines[int(lines.index(line) + 1)] or "HETATM" in lines[int(lines.index(line) + 1)]):
                        chain_id = line[21]
                        next_chain = lines[int(lines.index(line) + 1)][21]
                        if next_chain != chain_id:
                            corrected_lines.append(line)
                        elif next_chain == chain_id:
                            # Skip incorrect "TER" lines in the middle of a chain
                            continue
                    else:
                        corrected_lines.append(line)

            # Print the corrected PDB data
            for line in corrected_lines:
                f_out.write(line)
        f_out.close()


def fix_md_models(chain_ids, latest_ab):
    for rep in glob.glob("out-min.pdb"):
        rep_chain = "CHN_" + rep
        rep_score = "SCO_" + rep

        # Renaming chains
        with open(rep, "r") as f_in, open(rep_chain, "w") as f_out:
            chain_index = 0
            current_chain_id = ""
            prev_line_was_ter = False
            for line in f_in:
                if line.startswith("TER"):
                    prev_line_was_ter = True
                elif prev_line_was_ter:
                    chain_index += 1
                    if chain_index < len(chain_ids):
                        current_chain_id = chain_ids[chain_index]
                        prev_line_was_ter = False
                elif not current_chain_id:
                    current_chain_id = chain_ids[0]
                if line.startswith(("ATOM", "TER", "HETATM")):
                    if "LYN" in line:
                        new_line = line[:21].replace("LYN", "LYS").replace("HETATM", "ATOM  ") + current_chain_id + line[
                                                                                                                    22:]
                        f_out.write(new_line)
                    else:
                        new_line = line[:21] + current_chain_id + line[22:]
                        f_out.write(new_line)
                if line.startswith("END"):
                    f_out.write(line.replace("\n", ""))

            f_out.close()

            # Renumbering the antibody chains
            num_residues = []
            with open(latest_ab, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') and line[21] in ['H', 'L']:
                        residue = line[22:27].strip() + "_" + line[21]
                        if residue not in num_residues:
                            num_residues.append(residue)

            complex_2_residues = []
            with open(rep_chain, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') and line[21] in ['H', 'L']:
                        residue = int(line[22:26].strip())
                        if residue not in complex_2_residues:
                            complex_2_residues.append(residue)

            residue_map = {}
            for i in range(len(num_residues)):
                residue_map[complex_2_residues[i]] = num_residues[i]

            with open(rep_chain, 'r') as f_in, open(rep_score, 'w') as f_out:
                last_residue = -1
                while True:
                    line = f_in.readline()
                    if not line:
                        break
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        if line[21] in ['H', 'L']:
                            residue = int(line[22:26].strip())
                            if residue != last_residue:
                                if residue in residue_map:
                                    new_residue = residue_map[residue].replace("_" + line[21], "")
                                    last_residue = residue
                                else:
                                    new_residue = residue
                                match = re.match(r'^(\d+)(\D*)$', str(new_residue))
                                if match:
                                    digits = match.group(1)
                                    letters = match.group(2)
                                    new_residue = digits.rjust(4) + letters
                            if not any(c.isalpha() for c in new_residue):
                                line = line[:22] + '{:<4}'.format(new_residue) + line[26:]
                            else:
                                line = line[:22] + '{:<4}'.format(new_residue) + line[27:]
                        f_out.write(line)
                    else:
                        if "TER" in line:
                            f_out.write("TER\n")
                        if "END" in line:
                            f_out.write("END")


def new_complex(latest_ab, latest_complex, latest_ab_complex, temporary, renum, scoring_strictness, approval_threshold):
    # Checks the affinity score of the antibody structure against the target protein
    print("Starting evaluation protocol...")
    # Alignment
    with redirect_stdout(open("pymol_align.pml", "w", newline="")):
        print("load " + latest_complex)
        print("select ori_h, chain H")
        print("select ori_l, chain L")
        print("load " + latest_ab)
        print("align " + latest_ab.replace(".pdb", "") + ", " + latest_complex.replace(".pdb", ""))
        print("remove ori_h")
        print("remove ori_l")
        print("save " + latest_ab_complex + ", all")
        print("get_chains")
        print("select new_h, name CA and chain H")
        print("select new_l, name CA and chain L")
        print("select new_complex, name CA")

    run_pml_cmd = str(pymol_command + " -qc pymol_align.pml")
    out = subprocess.run(run_pml_cmd, shell=True, capture_output=True)
    pml_out = out.stdout.decode("utf-8")
    not_chains = ["\'", "[", "]", " ", ",", "\n"]
    chain_ids = [ch for ch in pml_out.split("cmd.get_chains:")[-1] if ch not in not_chains]

    # Refinement with OpenMM
    print("Refining complex side chains...")
    with ClearCache():
        refine(latest_ab_complex, latest_ab_complex)

    try:
        # Removal of insertion codes from pdb complex
        cmd_renum = "pdb_fixinsert " + latest_ab_complex + " > " + temporary + renum + "_" + latest_ab_complex
        subprocess.run(cmd_renum, shell=True, capture_output=True)
        # Protonation adjustment with propka
        cmd_prot = "pdb2pqr30 --with-ph " + ph + " --ff AMBER --ffout AMBER --pdb-output out_pdb2pqr.pdb --titration-state-method propka " + temporary + renum + "_" + latest_ab_complex + " out_pdb2pqr.pqr"
        # pdb4amber
        cmd_pdb4amber = "pdb4amber -i out_pdb2pqr.pdb -o out_pdb4amber-addH-preterfix.pdb"
        # Sampling with Amber MD
        cmd_md = "sh min_imp.sh"

        print("Running minimization...")
        subprocess.run(cmd_prot, shell=True, capture_output=True)
        subprocess.run(cmd_pdb4amber, shell=True, capture_output=True)
        # Fixing incorrect chain terminations
        ter_fix()
        subprocess.run(cmd_md, shell=True, capture_output=True)

        fix_md_models(chain_ids, latest_ab)
        scoring(latest_ab, latest_ab_complex, scoring_strictness, approval_threshold)

    except:
        print("Minimization failure. Check " + tested_canon[-1] + "_error.tar for more information.")
        rejected.append(tested_canon[-1])
        curr_log_rej = "rej_mod_" + curr_struct + ".log"
        updt_log_rej = open(curr_log_rej, "a", newline="")
        with redirect_stdout(updt_log_rej):
            print("\trejmut_" + tested_canon[-1] + "|Minimization error|FAIL")

        for rep in glob.glob("out-min.pdb"):
            rep_chain = "CHN_" + rep
            rep_score = "SCO_" + rep
            cmd_del_int = "rm " + rep_chain + " " + rep + " " + rep_score + " out* *.out *.info *.rst7 *.crd *.leap *.pqr *.parm7 *.mdinfo *.mdout *.nc temp*"
            cmd_bkp_error = "tar -cvf " + tested_canon[
                -1] + "_error.tar " + rep_chain + " " + rep + " " + rep_score + " out* *.out *.info *.rst7 *.crd *.leap *.pqr *.parm7 *.mdinfo *.mdout *.nc temp*"
            subprocess.run(cmd_bkp_error, shell=True, capture_output=True)
            subprocess.run(cmd_del_int, shell=True, capture_output=True)

        del tested_canon[-1]


def scoring(latest_ab, latest_ab_complex, scoring_strictness, approval_threshold):
        # Running CSM-AB or Rosetta (REF15 score) to calculate the interaction score of the new complex
        print("Scoring...")
        new_score = []
        latest_score = 0
        csm_ab_status = False
        starttime = time.time()

        if scoring_method == "csm":
            for rep in glob.glob("out-min.pdb"):
                rep_score = "SCO_" + rep
                while csm_ab_status == False:
                    csm = "curl https://biosig.lab.uq.edu.au/csm_ab/api/prediction_single -X POST -i -F pdb_file=@" + rep_score
                    out_csm = subprocess.run(csm, shell=True, capture_output=True)
                    out_csm_cln = out_csm.stdout.decode("utf-8")

                    if not "Please provide the Ab-Ag complex" in out_csm_cln and not out_csm_cln == "":
                        job_id = str(out_csm_cln.split("{\"job_id\": ")[-1].replace("}", "").replace("\"", "")).replace(
                            "\n",
                            "")
                        print("Starting CSM-AB. Job id is " + job_id)

                        csm_ab_done = False
                        starttime2 = time.time()

                        with redirect_stdout(open("CSM-AB.log", "a", newline="")):
                            print(out_csm_cln)
                            print("job id:" + job_id)

                        while csm_ab_done == False:
                            get_csm = "curl https://biosig.lab.uq.edu.au/csm_ab/api/prediction_single -X GET -F job_id=" + job_id
                            out_get_csm = subprocess.run(get_csm, shell=True, capture_output=True)
                            out_get_csm_cln = out_get_csm.stdout.decode("utf-8")
                            with redirect_stdout(open("CSM-AB.log", "a", newline="")):
                                print(out_get_csm_cln)
                                print("job id:" + job_id)
                            if "\"status\": \"PROCESSING\"" in out_get_csm_cln:
                                print("Running CSM-AB...")
                                time.sleep(10.0 - ((time.time() - starttime2) % 10.0))
                            if not "\"status\": \"PROCESSING\"" in out_get_csm_cln and not out_csm_cln == "":
                                if "\"typeofAb\":" in out_get_csm_cln:
                                    new_score.append(float(out_get_csm_cln.split("\"prediction\": ")[-1].replace("}", "")))
                                    csm_ab_done = True
                                    csm_ab_status = True
                                else:
                                    print("Error while getting results. Trying again in 10s")
                                    time.sleep(10.0 - ((time.time() - starttime2) % 10.0))
                            if out_csm_cln == "":
                                print("Error while getting results. Trying again in 10s")
                                time.sleep(10.0 - ((time.time() - starttime2) % 10.0))

                    else:
                        if out_csm_cln == "":
                            print("Error while sending job. Trying again in 10s")
                            time.sleep(10.0 - ((time.time() - starttime) % 10.0))
                        else:
                            if "Given PDB has no Antibody Sequence" in out_csm_cln:
                                print(
                                    "CSM-AB prediction error. Check the most recent AB/AG complex structure and if it is the file that was sent to CSM-AB. File: " + latest_ab_complex)
                            else:
                                print(
                                    "Unknown CSM-AB error. Skipping current modification. Check CSM-AB.log for more information")
                            with redirect_stdout(open("CSM-AB.log", "a", newline="")):
                                print(out_csm_cln)
                            break

        if scoring_method == "ref15":
            for rep in glob.glob("out-min.pdb"):
                rep_score = "SCO_" + rep
                with redirect_stdout(open("pymol_score.pml", "w", newline="")):
                    print("load " + rep_score)
                    print("select ab, chain H + chain L")
                    print("save antibody_score.pdb, ab")
                    print("remove ab")
                    print("save antigen_score.pdb, all")

                run_pml_cmd = str(pymol_command + " -qc pymol_score.pml")
                subprocess.run(run_pml_cmd, shell=True, capture_output=True)

                scored_complex = pyrosetta.pose_from_pdb(rep_score)
                ab_only = pyrosetta.pose_from_pdb("antibody_score.pdb")
                ag_only = pyrosetta.pose_from_pdb("antigen_score.pdb")

                scorefxn = get_fa_scorefxn()
                complex_en = float(scorefxn(scored_complex))
                ab_en = float(scorefxn(ab_only))
                ag_en = float(scorefxn(ag_only))
                new_score.append(float(complex_en - ab_en - ag_en))
                subprocess.run("rm antibody_score.pdb antigen_score.pdb", shell=True, capture_output=True)

        # Score comparison
        if new_score:
            av_new_score = new_score[-1]
            curr_log = "swap_" + curr_struct + ".log"
            open_curr_log = open(curr_log, "r").readlines()

            for z in open_curr_log:
                if z.startswith(file):
                    latest_score = (z.split("|")[1]).replace("\n", "")

            kT = 0
            if scoring_method == "ref15":
                kT = 1.0
            if scoring_method == "csm":
                kT = 0.593

            approval_condition = False

            if scoring_strictness == "metro":
                approval_prob = math.exp((float(latest_score) - av_new_score) / kT)
                metropolis = random.random()
                if (float(av_new_score) < float(latest_score)) or (
                        metropolis < approval_prob and float(av_new_score) != float(latest_score)):
                    approval_condition = True
                elif (float(av_new_score) > float(latest_score) and metropolis > approval_prob) or float(
                        av_new_score) == float(
                        latest_score):
                    approval_condition = False

            if scoring_strictness == "normal":
                if float(av_new_score) < float(latest_score):
                    approval_condition = True
                elif float(av_new_score) > float(latest_score) or float(av_new_score) == float(latest_score):
                    approval_condition = False

            if scoring_strictness == "strict":
                if float(av_new_score) - float(latest_score) <= approval_threshold:
                    approval_condition = True
                elif float(av_new_score) - float(latest_score) > approval_threshold:
                    approval_condition = False

            if approval_condition:
                file_number = (file.split("_")[4]).replace(".anarci", "")
                new_file_number = '{:06}'.format(int(file_number) + 1)
                file_start = "_".join([str(item) for item in [p for p in file.split("_")[:4]]])
                new_file_name = file_start + "_" + str(new_file_number) + ".anarci"

                updt_log = open(curr_log, "a", newline="")
                with redirect_stdout(updt_log):
                    print("\tMUTATION_" + tested_canon[-1] + "|" + str(av_new_score) + "|PASS")
                    print(new_file_name + "|" + str(av_new_score))

                cmd = str("ANARCI -i " + latest_ab.replace(".pdb", ".fasta") + " -s martin -o " + new_file_name + " --assign_germline -p 8 --use_species human")
                out = subprocess.run(cmd, shell=True, capture_output=True)
                log9 = (out.stderr).decode("utf-8")
                final.append(new_file_name)
                final_sco.append(av_new_score)

                log3 = open("anarci_out.log", "a", newline="")
                with redirect_stdout(log3):
                    print(new_file_name)
                    print(log9)

                print(tested_canon[-1] + "|" + str(av_new_score) + "|Modification approved")
                condition.append("PASS")

                rename_cmd = str("cp " + latest_ab_complex + " complex_" + new_file_name.replace(".anarci", ".pdb"))
                rename_cmd2 = str(
                    "cp " + latest_ab.replace(".pdb", ".fasta") + " " + new_file_name.replace(".anarci", ".fasta"))
                rename_cmd3 = str("cp " + latest_ab + " " + new_file_name.replace(".anarci", ".pdb"))
                subprocess.run(rename_cmd, shell=True, capture_output=True)
                subprocess.run(rename_cmd2, shell=True, capture_output=True)
                subprocess.run(rename_cmd3, shell=True, capture_output=True)

            if not approval_condition:
                print(tested_canon[-1] + "|" + str(av_new_score) + "|Modification rejected")
                updt_log = open(curr_log, "a", newline="")
                with redirect_stdout(updt_log):
                    print("\tMUTATION_" + tested_canon[-1] + "|" + str(av_new_score) + "|FAIL")

            # Removes scoring files from the current cycle
            for rep in glob.glob("out-min.pdb"):
                rep_chain = "CHN_" + rep
                rep_score = "SCO_" + rep
                cmd_del_int = "rm " + rep_chain + " " + rep + " " + rep_score + " out* *.out *.info *.rst7 *.crd *.leap *.pqr *.parm7 *.mdinfo *.mdout *.nc out* temp*"
                subprocess.run(cmd_del_int, shell=True, capture_output=True)


##########
# Reading configuration file to determine how many swap cycles will be attempted, the probabilities of a CDR or chain
# being chosen for modification, whether or not the conformation of each new CDR will be compared to the original
# structure, and the RMSD threshold above which a CDR sequence will be discarded during a conformation check
check_conformation = "cfg:mut|check_conf"
conf_rmsd_limit = "cfg:mut|conf_rmsd_limit"
n_cycles = "cfg:mut|n_cycles"
cdr_prob = "cfg:mut|cdr_prob"
chain_prob = "cfg:mut|chain_prob"
h3_cap = "cfg:mut|h3_cap"
res_prob_mode = "cfg:mut|res_prob_mode"
steps = "cfg:steps"
scoring_strictness = "cfg:scoring_strictness"
approval_threshold = 0
ph = "7"
pymol_command = ""
scoring_method = "cfg:scoring_method"

cfg = open("swap_settings.cfg").readlines()
for setting in cfg:
    if "mut|check_conf" in setting:
        check_conformation = str(setting.split("=")[1].replace("\n", ""))
    if "mut|conf_rmsd_limit" in setting:
        conf_rmsd_limit = float(str(setting.split("=")[1].replace("\n", "")))
    if "mut|n_cycles" in setting:
        n_cycles = int(setting.split("=")[1].replace("\n", ""))
    if "mut|cdr_prob" in setting:
        cdr_prob = str(setting.split("=")[1].replace("\n", ""))
    if "mut|chain_prob" in setting:
        chain_prob = str(setting.split("=")[1].replace("\n", ""))
    if "mut|res_prob_mode" in setting:
        res_prob_mode = str(setting.split("=")[1].replace("\n", ""))
    if "mut|h3_cap" in setting:
        h3_cap = str(setting.split("=")[1].replace("\n", ""))
    if "minimization_pH" in setting:
        ph = str(setting.split("=")[1].replace("\n", ""))
    if "pymol_command" in setting:
        pymol_command = str(setting.split("=")[1].replace("\n", ""))
    if "scoring_method" in setting:
        scoring_method = str(setting.split("=")[1].replace("\n", ""))
    if "steps" in setting:
         steps = str(setting.split("=")[1].replace("\n", "").replace("\r",""))
    if "scoring_strictness" in setting:
        scoring_strictness = str(setting.split("=")[1].replace("\n", "").replace("\r",""))
    if "approval_threshold" in setting:
        if scoring_strictness == "strict":
            approval_threshold = float(setting.split("=")[1].replace("\n", "").replace("\r",""))

print("Starting new step: Mutagenesis")

if scoring_method == "ref15":
    from pyrosetta import *
    init("-mute all")

base_structure = []
# Stores which and how many antibodies are being developed simultaneously (numbered 001, 002...)
structures = []

step_list = steps.split("|")
if "" in step_list:
    step_list.remove("")
bs_index = step_list.index("mut")
base_structure.append(str(bs_index + 1))
base_structure.append(str(bs_index + 2))

for struc in glob.glob("bs" + base_structure[0] + "_*"):
    structure = struc.split("_")[1]
    if not structure in structures:
        structures.append(structure)

for st in structures:
    # Stores which modifications have already been tested
    tested_canon = []
    # Stores the IDs of sequences that were rejected for producing CDRs with a different conformation compared to the original structure
    rejected = []
    # Stores the name of the final structure, after all modifications, and its score
    final = []
    final_sco = []
    # Checks if there have been any errors during the modification process
    error = False
    # Stores the CDR lengths of the first bs4 structure
    H_cdr1_canon = []
    H_cdr2_canon = []
    H_cdr3_canon = []
    L_cdr1_canon = []
    L_cdr2_canon = []
    L_cdr3_canon = []

    print("New structure")

    # Keeps track of the current modification supercycle
    supercycle = 1

    # Main loop; this "while" loop covers the entire modification/scoring cycle.
    # This loop will be repeated for n cycles.
    while len(tested_canon) < n_cycles:
        latest = []
        num_list = []
        run_list = []
        sco_list = []
        condition = []
        # Finds the antibody structure with the best score in the main log file, which will serve as the starting point for the current modification cycle
        for file in glob.glob("swap_" + st + ".log"):
            open_log = open(os.path.join(os.getcwd(), file), "r", newline='')
            pass_count = count(start=2, step=1)
            for line in open_log:
                if line.startswith("bs" + base_structure[0] + "_"):
                    name = (line.split("|")[0]).replace("\n", "")
                    score = (line.split("|")[1]).replace("\n", "")
                    progress = (name.split("_")[-1]).replace(".anarci", "")
                    run = name.split(progress)[0]
                    run_list.append(run)
                    num_list.append(progress)
                    sco_list.append(float(score))
                if line.startswith("\tMUTATION_"):
                    tested = (line.split("|")[0]).replace("\tMUTATION_", "")
                    if "PASS" in line:
                        supercycle = next(pass_count)
                    if not tested in tested_canon:
                        tested_canon.append(tested)

            latest.clear()
            condition.clear()

            newest = run_list[-1] + num_list[-1] + ".anarci"
            best = sco_list[-1]

            latest.append(newest)
            final.append(newest)
            final_sco.append(best)
            run_list.clear()
            num_list.clear()
            open_log.close()
            error = False

        for rej in glob.glob("rej_mod_" + st + ".log"):
            open_log_rej = open(os.path.join(os.getcwd(), rej), "r", newline='')
            for rej_line in open_log_rej:
                if rej_line.startswith("\trejmut_"):
                    rejected_mut = (rej_line.split("|")[0]).replace("\trejmut_", "")
                    if not rejected_mut in rejected:
                        rejected.append(rejected_mut)

        for new in latest:
            for file in glob.glob(new):
                # Opens anarci and log files
                open_num = open(file, "r").readlines()
                curr_struct = file.split("_")[1]
                log = (open("swap_" + curr_struct + ".log", "r")).readlines()
                id = file.split("_")[2]
                # List of chains and CDR numbers.
                Ch_list = ["H", "L", "H3"]
                HL_list = ["1", "2", "3"]
                # Residues that define each CDR, according to the Martin/Enhanced Chothia scheme
                CDRH1_res = [26, 32]
                CDRH2_res = [52, 56]
                CDRH3_res = [95, 102]
                CDRL1_res = [24, 34]
                CDRL2_res = [50, 56]
                CDRL3_res = [89, 97]

                # Stores the residues of each CDR from the newest antibody
                for ch in [cdr for cdr in Ch_list if cdr != "H3"]:
                    for cdr_num in HL_list:
                        curr_cdr_lst = locals()[ch + "_cdr" + cdr_num + "_canon"]
                        curr_cdr_lst.clear()
                    for seq in open_num:
                        if seq.startswith(ch):
                            for num in HL_list:
                                CDR_id = ch + num
                                curr_cdr_res = locals()["CDR" + CDR_id + "_res"]
                                curr_cdr_lst = locals()[ch + "_cdr" + num + "_canon"]
                                num_res = int(seq.split(" ")[1])
                                amino = (seq.split(" ")[-1]).replace("\n", "")
                                full_num_res = (seq.split(" ")[1] + seq.split(" ")[-2]).replace(" ", "")
                                if num_res >= curr_cdr_res[0] and num_res <= curr_cdr_res[-1] and not amino == "-":
                                    curr_cdr_lst.append(full_num_res + "|" + amino)

                # Chain being modified is chosen according to the specified probabilities from the cfg file
                swap_Ch = []
                H3 = False

                if h3_cap == "no":
                    swap_Ch.append(random.choices([chn for chn in Ch_list if chn != "H3"],
                                             weights=(int(chain_prob.split("|")[0]), int(chain_prob.split("|")[1])),
                                             k=1)[-1])
                if h3_cap == "yes":
                    # Approximately half of all mutagenesis cycles are allocated to CDR H3
                    chosen_chain = random.choices(Ch_list,
                                             weights=(int(chain_prob.split("|")[0]), int(chain_prob.split("|")[1]), (int(chain_prob.split("|")[0]) + int(chain_prob.split("|")[1]))),
                                             k=1)

                    if not "H3" in chosen_chain:
                        swap_Ch.append(chosen_chain[-1])
                    if "H3" in chosen_chain:
                        swap_Ch.append("H")
                        H3 = True

                # The pieces that will form the modified antibody sequence are stored here, along with the length of the CDR that is being grafted
                newseq_lst = []
                # Header information for each chain is stored here
                main_header = []
                other_chain = []
                print("Starting new cycle")
                print("Modifying CDR...")

                for chain in swap_Ch:
                    # CDR being modified is chosen according to the specified probabilities in the cfg file
                    swapH = []
                    swapL = random.choices(HL_list, weights=(
                    int(cdr_prob.split("|")[3]), int(cdr_prob.split("|")[4]), int(cdr_prob.split("|")[5])), k=1)

                    if h3_cap == "no":
                        swapH.append(random.choices(HL_list, weights=(
                        int(cdr_prob.split("|")[0]), int(cdr_prob.split("|")[1]), int(cdr_prob.split("|")[2])), k=1)[-1])

                    if h3_cap == "yes":
                        if not H3:
                            swapH.append(random.choices([cdr for cdr in HL_list if cdr != "3"], weights=(
                            int(cdr_prob.split("|")[0]), int(cdr_prob.split("|")[1])), k=1)[-1])
                        if H3:
                            swapH.append("3")

                    for x in locals()["swap" + chain]:
                        # The numbers collected below will determine the probability that each position of the chosen CDR has of being chosen for mutation (diver_list),
                        # and the probability that each of the 20 possible residues has of being introduced in the selected position (freq_list)
                        freq_list = []
                        diver_list = []

                        if res_prob_mode == "probabilistic":
                            if path.exists(os.path.dirname(os.path.abspath(__file__)) + "/mem_mutation_data/" + chain + x + "-" + str(len(locals()[chain + "_cdr" + x + "_canon"])) + "_mem_entropy.csv"):
                                canon_file = os.path.dirname(os.path.abspath(__file__)) + "/mem_mutation_data/" + chain + x + "-" + str(len(locals()[chain + "_cdr" + x + "_canon"])) + "_mem_entropy.csv"
                                open_canon = (open(canon_file, "r")).readlines()

                                for y in open_canon:
                                    if not y.startswith("Index|"):
                                        res_freq = {}
                                        diversity = y.split(";")[1]
                                        frequencies = ast.literal_eval("{" + y.split(";")[2].replace("\n", "") + "}")
                                        diver_list.append(float(diversity))
                                        res_freq.update(frequencies)
                                        freq_list.append(res_freq)

                            else:
                                print("No data for CDR " + chain + x + " of length " + str(
                                    len(locals()[chain + "_cdr" + x + "_canon"])) + ". Selecting random mutation...")
                                # Cysteine was excluded from random mode, but can be added back wihtout issue in this dictionary, if desired
                                res_freq = {"A":1,"D":1,"E":1,"F":1,"G":1,"H":1,"I":1,"K":1,"L":1,"M":1,"N":1,"P":1,"Q":1,"R":1,"S":1,"T":1,"V":1,"Y":1,"W":1}
                                for i in range(1, (int(len(locals()[chain + "_cdr" + x + "_canon"])) + 1)):
                                    diver_list.append(1)
                                    freq_list.append(res_freq)

                        if res_prob_mode == "random":
                            # Cysteine was excluded from random mode, but can be added back wihtout issue in this dictionary, if desired
                            res_freq = {"A": 1, "D": 1, "E": 1, "F": 1, "G": 1, "H": 1, "I": 1, "K": 1, "L": 1, "M":1, "N": 1,
                                        "P": 1, "Q": 1, "R": 1, "S": 1, "T": 1, "V": 1, "Y": 1, "W": 1}
                            for i in range(1, (int(len(locals()[chain + "_cdr" + x + "_canon"])) + 1)):
                                diver_list.append(1)
                                freq_list.append(res_freq)

                        positions = range(1, (int(len(locals()[chain + "_cdr" + x + "_canon"])) + 1))

                        mutation_id = "chosen mutation"
                        orig_res_num = "number"
                        choose_mutation = "other number"
                        saturation_check = count(start=1, step=1)
                        most_mut = (int(len(locals()[chain + "_cdr" + x + "_canon"]))*0.5)*20

                        while mutation_id == "chosen mutation":
                            choose_position = random.choices(positions, weights=(diver_list), k=1)
                            position_index = choose_position[-1] - 1
                            chosen_freq_dict = freq_list[position_index]
                            choose_mutation = random.choices(list(chosen_freq_dict.keys()),
                                                             weights=(list(chosen_freq_dict.values())), k=1)[-1]
                            orig_res_num = locals()[chain + "_cdr" + x + "_canon"][position_index]
                            sel_mut_id = str(supercycle) + "-" + chain + x + "-" + orig_res_num.split("|")[1] + "-" + orig_res_num.split("|")[0] + "-" + choose_mutation[-1]
                            saturation_num = int(next(saturation_check))

                            if not sel_mut_id in tested_canon and not sel_mut_id in rejected and not orig_res_num.split("|")[1] == choose_mutation[-1]:
                                mutation_id = sel_mut_id
                            else:
                                if saturation_num > most_mut:
                                    if saturation_num == most_mut*2:
                                        print("All possible mutations for CDR " + chain + x + " of this length have been tested. Modifying different CDR...")
                                        break

                                    print("Most possible mutations for CDR " + chain + x + " of this length have been tested. Selecting random mutation...")
                                    diver_list.clear()
                                    freq_list.clear()
                                    res_freq = {"A": 1, "D": 1, "E": 1, "F": 1, "G": 1, "H": 1, "I": 1, "K": 1, "L": 1,
                                                "M": 1, "N": 1, "P": 1, "Q": 1, "R": 1, "S": 1, "T": 1, "V": 1, "Y": 1,
                                                "W": 1}

                                    for i in range(1, (int(len(locals()[chain + "_cdr" + x + "_canon"])) + 1)):
                                        diver_list.append(1)
                                        freq_list.append(res_freq)

                        if not mutation_id in tested_canon and not mutation_id in rejected and not mutation_id == "chosen mutation":
                            print("Mutation selected : " + mutation_id)
                            tested_canon.append(mutation_id)

                            for d in open_num:
                                other_chain_id = [c for c in Ch_list if not c == chain and not c == "H3"][-1]
                                # Uses the anarci file to determine where in the original sequence the new CDR sequence should be introduced
                                if d.startswith("# " + id):
                                    ch_header = (d.split("# ")[1]).replace("\n", "")
                                    main_header.append(ch_header)

                                if d.startswith(chain):
                                    num_res = d.split(" ")[1]
                                    amino = (d.split(" ")[-1]).replace("\n", "")
                                    full_num_res = (d.split(" ")[1] + d.split(" ")[-2]).replace(" ", "")

                                    if full_num_res != orig_res_num.split("|")[0]:
                                        newseq_lst.append(amino)
                                    if full_num_res == orig_res_num.split("|")[0]:
                                        newseq_lst.append(choose_mutation)

                                if d.startswith(other_chain_id):
                                    amino = (d.split(" ")[-1]).replace("\n", "")
                                    other_chain.append(amino)

                        if main_header:
                            # If the sequence was sucessfully modified, the length of the new CDR will be checked and the modified sequence will be modeled
                            newseq = "".join(item for item in [w for w in newseq_lst if not w == "-"])

                            latest_ab = "MUTATION_" + mutation_id + "_" + st + "_" + file.split("_")[2] + "_" + file.split("_")[-1].replace(".anarci", ".pdb")
                            latest_complex = "complex_" + file.replace(".anarci", ".pdb")
                            latest_ab_complex = "complex_" + latest_ab
                            temporary = "temp"
                            renum = "_renum"

                            try:
                                print("Modelling new antibody...")
                                model_mutated()
                                conf_check(latest_ab)
                            except:
                                # If an error occurs during modelling, that step will be repeated
                                print("Modelling error. Restarting cycle...")
                                error = True
                                continue

                            if not error:
                                # If the modelling step gets completed sucessfully, the new antibody structure will be scored according to its affinity with the target protein
                                new_complex(latest_ab, latest_complex, latest_ab_complex, temporary, renum, scoring_strictness, approval_threshold)

                if "PASS" in condition:
                    break

    if not error:
        # At the end of all modification cycles, the structure and sequence of the best antibody produced are copied and renamed, in preparation for the next modification step
        print("Finalizing current structure...")
        final_num = ""
        if not step_list[-1] == "mut":
            final_num = str(final[-1]).replace("bs" + base_structure[0], "bs" + base_structure[1])
        if step_list[-1] == "mut":
            final_num = str(final[-1]).replace("bs" + base_structure[0], "FINAL")

        final_structure = final_num.replace(".anarci", ".pdb")
        final_sequence = final_num.replace(".anarci", ".fasta")
        final_complex = "complex_" + final_num.replace(".anarci", ".pdb")
        finalize_num = str("cp " + str(final[-1]) + " " + final_num)
        finalize_str = str("cp " + str((final[-1]).replace(".anarci", ".pdb")) + " " + final_structure)
        finalize_seq = str("cp " + str((final[-1]).replace(".anarci", ".fasta")) + " " + final_sequence)
        finalize_cpx = str("cp complex_" + str((final[-1]).replace(".anarci", ".pdb")) + " " + final_complex)
        subprocess.run(finalize_num, shell=True, capture_output=True)
        subprocess.run(finalize_str, shell=True, capture_output=True)
        subprocess.run(finalize_seq, shell=True, capture_output=True)
        subprocess.run(finalize_cpx, shell=True, capture_output=True)

        curr_log2 = "swap_" + st + ".log"
        updt_log2 = open(curr_log2, "a", newline="")
        with redirect_stdout(updt_log2):
            print(final_num + "|" + str(final_sco[-1]))

        tested_canon.clear()
