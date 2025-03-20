
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


class ClearCache:
    def __enter__(self):
        torch.cuda.empty_cache()

    def __exit__(self, exc_type, exc_val, exc_tb):
        torch.cuda.empty_cache()


def gen_script(steps, vers):
    step_list = steps.replace("\n","").replace("\r","").split("|")
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


def model_ab():
    print("Modelling new antibody...")
    sequences = {}

    header_H = ">" + input_name + "|H|" + str(cdr_lengths[0]) + "|" + str(cdr_lengths[1]) + "|" + str(cdr_lengths[2]) + "|" + germlines[0] + "|"
    header_L = ">" + input_name + "|L|" + str(cdr_lengths[3]) + "|" + str(cdr_lengths[4]) + "|" + str(cdr_lengths[5]) + "|" + \
               germlines[1] + "|"
 
    sequences["H"] = "".join([str(item) for item in [b for b in H_seq if not b == "-"]])
    sequences["L"] = "".join([str(item) for item in [b for b in L_seq if not b == "-"]])
    
    output_file = filename + ".pdb"

    with ClearCache():
        antibody = predictor.predict(sequences)
        antibody.save(output_file, check_for_strained_bonds=False)
    
    make_fasta = output_file.replace(".pdb", ".fasta")
    with redirect_stdout(open(make_fasta, "w", newline="")):
        print(header_H)
        print("".join([str(item) for item in [b for b in H_seq if not b == "-"]]))
        print(header_L)
        print("".join([str(item) for item in [b for b in L_seq if not b == "-"]]))


def ter_fix():
    #Fixing possible incorrect chain terminations introduced by pdb4amber
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


def new_complex(latest_ab, latest_complex, latest_ab_complex, temporary, no_hetatm, renum):
    # Checks the affinity score of the antibody structure against the target protein
    print("Starting evaluation protocol...")
    # Alignment and assembly of new Ab-Ag complex
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

    # Removal of heteroatoms
    pdb_line_list = open(latest_ab_complex).readlines()
    for non_hetatm in pdb_line_list:
        if not non_hetatm.startswith("HETATM"):
            save = open(temporary + no_hetatm + "_" + latest_ab_complex, "a", newline="")
            with redirect_stdout(save):
                print(non_hetatm.replace("\n", "").replace("\r", ""))

    # Refinement with OpenMM
    print("Refining complex side chains...")
    with ClearCache():
        refine(temporary + no_hetatm + "_" + latest_ab_complex, temporary + no_hetatm + "_" + latest_ab_complex)

    try:
        # Removal of insertion codes from pdb complex
        cmd_renum = "pdb_fixinsert " + temporary + no_hetatm + "_" + latest_ab_complex + " > " + temporary + no_hetatm + renum + "_" + latest_ab_complex
        subprocess.run(cmd_renum, shell=True, capture_output=True)
        # Protonation adjustment with propka
        cmd_prot = "pdb2pqr30 --with-ph " + ph + " --ff AMBER --ffout AMBER --pdb-output out_pdb2pqr.pdb --titration-state-method propka " + temporary + no_hetatm + renum + "_" + latest_ab_complex + " out_pdb2pqr.pqr"
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

    except:
        print("Minimization failure. Check amber log files and aligned complex for more information.")
        cmd_bkp_error = "tar -cvf " + latest_ab + " _error.tar out* *.out *.info *.rst7 *.crd *.leap *.pqr *.parm7 *.mdinfo *.mdout *.nc temp*"
        subprocess.run(cmd_bkp_error, shell=True, capture_output=True)

    fix_md_models(chain_ids, latest_ab)


def scoring():
    # Running CSM-AB or Rosetta (REF15 score) to calculate the interaction score of the new complex
    print("Scoring...")
    new_score = []
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
                    job_id = str(out_csm_cln.split("{\"job_id\": ")[-1].replace("}", "").replace("\"", "")).replace("\n",
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
                                "CSM-AB prediction error. Check the most recent AB/AG complex structure and if it is the file that was sent to CSM-AB. File: complex_" + latest_ab)
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
        curr_log = "swap_" + number + ".log"
        updt_log = open(curr_log, "a", newline="")
        with redirect_stdout(updt_log):
            print(filename + ".anarci|" + str(float(av_new_score)))

        # Use minimized file as starting point
        cmd_update = "cp SCO_out-min.pdb " + complex_filename
        subprocess.run(cmd_update, shell=True, capture_output=True)

    # Remove scoring files from the current cycle
    for rep in glob.glob("out-min.pdb"):
        rep_chain = "CHN_" + rep
        rep_score = "SCO_" + rep
        cmd_del_int = "rm " + rep_chain + " " + rep + " " + rep_score + " out* *.out *.info *.rst7 *.crd *.leap *.pqr *.parm7 *.mdinfo *.mdout *.nc out* temp*"
        subprocess.run(cmd_del_int, shell=True, capture_output=True)


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
        input_name = str(setting.split("=")[1].replace("\n", "")).replace(".pdb","").replace(".fasta","")
    if "prepare|seldon_version" in setting:
        vers = str(setting.split("=")[1].replace("\n", ""))

print("Starting pipeline: Input preparation")

if scoring_method == "ref15":
    from pyrosetta import *
    init("-mute all")

CDRH1_res = [26, 32]
CDRH2_res = [52, 56]
CDRH3_res = [95, 102]
CDRL1_res = [24, 34]
CDRL2_res = [50, 56]
CDRL3_res = [89, 97]

counter = ("%03i" % i for i in count(1))

for amount in range(1,(int(num_struc) + 1)):
    fasta = input_name + ".fasta"
    number = next(counter)
    filename = "bs1_" + number + "_" + input_name.replace("_","-") + "_" + code_name + "_000001"
    complex_filename = "complex_" + filename + ".pdb"

    # Run anarci to collect information about the antibody
    cmd = str("ANARCI -i " + fasta + " -s martin -o " + filename + ".anarci -p 8 --assign_germline --use_species human")
    out = subprocess.run(cmd, shell=True, capture_output=True)

    log_anar = (out.stderr).decode("utf-8")
    log_an_prt = open("anarci_out.log", "a", newline="")
    with redirect_stdout(log_an_prt):
        print(log_anar)

    open_num = open(os.path.join(os.getcwd(), filename + ".anarci"), "r", newline='')
    germlines = []
    cdr_seqs = {"H1":[],"H2":[],"H3":[],"L1":[],"L2":[],"L3":[]}
    cdr_lengths = []
    H_seq = []
    L_seq = []

    for num in open_num:
        if num.startswith("#|species|v_gene|v_identity|j_gene|j_identity|"):
            germlines.append(next(open_num).split("|")[2])

        if num.startswith("H") or num.startswith("L"):
            ch = num.split(" ")[0]
            num_res = num.split(" ")[1]
            amino = (num.split(" ")[-1]).replace("\n", "")
            
            locals()[ch + "_seq"].append(amino)

            for i in range(1,4):
                if int(num_res) >= locals()["CDR" + ch + str(i) + "_res"][0] and int(num_res) <= locals()["CDR" + ch + str(i) + "_res"][-1]:
                    cdr_seqs.get(ch + str(i)).append(amino)

    for key, value in cdr_seqs.items():
        cdr_lengths.append(len(value))

    model_ab()

    latest_ab = filename + ".pdb"
    latest_complex = input_name + ".pdb"
    latest_ab_complex = "complex_" + latest_ab
    temporary = "temp"
    no_hetatm = "nohetatm"
    renum = "_renum"

    new_complex(latest_ab, latest_complex, latest_ab_complex, temporary, no_hetatm, renum)
    scoring()
    gen_script(steps, vers)

    # Updating the .anarci file
    cmd = str("ANARCI -i " + filename + ".fasta -s martin -o " + filename + ".anarci -p 8 --assign_germline --use_species human")
    out = subprocess.run(cmd, shell=True, capture_output=True)

print("Input preparation completed")
