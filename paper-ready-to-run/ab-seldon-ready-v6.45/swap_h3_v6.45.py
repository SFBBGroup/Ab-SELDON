
import glob
import os
import subprocess
import random
import math
import linecache
import time
import re
import torch
from contextlib import redirect_stdout
from ImmuneBuilder import ABodyBuilder2
from ImmuneBuilder.refine import refine
predictor = ABodyBuilder2(numbering_scheme="martin")


class ClearCache:
    def __enter__(self):
        torch.cuda.empty_cache()

    def __exit__(self, exc_type, exc_val, exc_tb):
        torch.cuda.empty_cache()


def mod_header():
    # checks the length of the CDRs and stores them on a modified header
    lengths = {"H1_len": (main_header[0]).split("|")[2], "H2_len": (main_header[0]).split("|")[3],
               "H3_len": (main_header[0]).split("|")[4], "L1_len": (main_header[1]).split("|")[2],
               "L2_len": (main_header[1]).split("|")[3], "L3_len": (main_header[1]).split("|")[4]}

    headers = {"start_H": str("|".join([str(item) for item in [m for m in (main_header[0]).split("|")[:2]]])), "start_L": str("|".join([str(item) for item in [m for m in (main_header[1]).split("|")[:2]]])), "end_H": str("|".join([str(item) for item in [g for g in (main_header[0]).split("|")[5:]]])), "end_L": str("|".join([str(item) for item in [g for g in (main_header[1]).split("|")[5:]]]))}

    lengths["H3_len"] = str(new_CDR_len[-1])

    new_header_h = headers.get("start_H") + "|" + lengths.get("H1_len") + "|" + lengths.get("H2_len") + "|" + lengths.get(
        "H3_len") + "|" + headers.get("end_H")
    new_header_l = headers.get("start_L") + "|" + lengths.get("L1_len") + "|" + lengths.get("L2_len") + "|" + lengths.get(
        "L3_len") + "|" + headers.get("end_L")

    new_header_light.append(new_header_l)
    new_header_heavy.append(new_header_h)


def model_ab():
    # models the modified antibody structure
    newseq_raw = start + new_CDR + end
    newseq = "".join([str(item) for item in [b for b in newseq_raw if not b == "-"]])

    sequences = {}

    sequences["H"] = newseq
    sequences["L"] = "".join([str(item) for item in [b for b in other_chain if not b == "-"]])

    output_file = file.replace(".anarci", "") + "_" + tested_canon[-1] + ".pdb"

    with ClearCache():
        antibody = predictor.predict(sequences)
        antibody.save(output_file, check_for_strained_bonds=False)

    make_fasta = output_file.replace(".pdb", ".fasta")
    with redirect_stdout(open(make_fasta, "w", newline="")):
        print(">" + new_header_heavy[0])
        print(newseq)
        print(">" + new_header_light[0])
        print("".join([str(item) for item in [b for b in other_chain if not b == "-"]]))


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


def qual_ctrl(latest_ab, latest_complex, latest_ab_complex, temporary, renum, scoring_strictness, approval_threshold):
    # Checks the quality of the produced structure, looking at the CDR H3 RMSD (deviation estimated by AbodyBuilder2);
    # The structure gets discarded if the RMSD is above the threshold specified on the cfg file
    # Otherwise, the new complex is assembled and scored
    read_structure = open(file.replace(".anarci", "") + "_" + tested_canon[-1] + ".pdb", "r").readlines()

    res_list = []
    rmsd_list = []
    rmsd_level = []

    for atom in read_structure:
        if atom.startswith("ATOM"):
            res_num = atom[22:27].replace(" ", "")
            res_num_raw = atom[22:26].replace(" ", "")
            chain = atom[21]
            res_type = atom[17:21].replace(" ", "")
            rmsd = atom[61:66].replace(" ", "")
            res_id = chain + "_" + res_type + "_" + res_num

            if chain == "H" and int(res_num_raw) >= CDRH3_res[0] and int(res_num_raw) <= CDRH3_res[-1]:
                if not res_id in res_list:
                    res_list.append(res_id)
                    rmsd_list.append(rmsd)

    for deviation in rmsd_list:
        if float(deviation) > h3_rmsd_limit:
            rmsd_level.append("high")
        else:
            rmsd_level.append("low")

    if rmsd_level and not "high" in rmsd_level:
        new_complex(latest_ab, latest_complex, latest_ab_complex, temporary, renum, scoring_strictness, approval_threshold)

    if "high" in rmsd_level:
        rejected.append(tested_canon[-1])
        curr_log_rej = "rej_mod_" + curr_struct + ".log"
        updt_log_rej = open(curr_log_rej, "a", newline="")
        with redirect_stdout(updt_log_rej):
            print("\trejh3_" + tested_canon[-1] + "|LQ model|FAIL")
        del tested_canon[-1]
        print("Low H3 quality. Model discarded")


def make_new_seq():
    #chooses a random CDR H3 sequence among those with the selected length in the database
    if number_tested < canon_numbers:
        go = False
        while not go:
            line_range = int(len(open_canon))
            file_line = random.randrange(1, line_range, 2)
            find_header = (linecache.getline(canon_file, file_line)).replace("\n", "")
            find_seq = (linecache.getline(canon_file, file_line + 1)).replace("\n", "")

            line_id = find_header.split("|")[2]

            canon = str(choose_length + "_" + line_id.replace("\n", ""))
            mem = find_seq.replace("\n", "")

            if canon not in tested_canon and canon not in rejected:
                cdr_header.append(canon)
                new_CDR.append(mem)
                new_CDR_len.append(len(mem))
                tested_canon.append(canon)
                go = True
            else:
                new_CDR.clear()
                pass

        if new_CDR:
            for d in open_num:
                # uses the anarci file to determine where in the original sequence the new CDR sequence should be introduced
                if d.startswith("# " + id):
                    ch_header = (d.split("# ")[1]).replace("\n", "")
                    main_header.append(ch_header)

                if d.startswith("H"):
                    num_res = d.split(" ")[1]
                    amino = (d.split(" ")[-1]).replace("\n", "")

                    if int(num_res) < CDRH3_res[0]:
                        start.append(amino)

                    if int(num_res) > CDRH3_res[-1]:
                        end.append(amino)

                if d.startswith("L"):
                    amino = (d.split(" ")[-1]).replace("\n", "")
                    other_chain.append(amino)

    else:
        print("All database sequences for CDR " + choose_length + " have been tested. Trying another length...")


def read_log_init():
    # finds the antibody structure with the best score in the main log file, which will serve as the starting point for the current modification cycle
    for file in glob.glob("swap_" + st + ".log"):
        open_log = open(os.path.join(os.getcwd(), file), "r", newline='')
        for line in open_log:
            if line.startswith("bs" + base_structure[0] + "_"):
                name = (line.split("|")[0]).replace("\n", "")
                score = (line.split("|")[1]).replace("\n", "")
                progress = (name.split("_")[-1]).replace(".anarci", "")
                run = name.split(progress)[0]
                run_list.append(run)
                num_list.append(progress)
                sco_list.append(float(score))
            if line.startswith("\t") and "H3-" in line:
                tested = (line.split("|")[0]).replace("\t", "")
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
        error.clear()

    for rej in glob.glob("rej_mod_" + st + ".log"):
        open_log_rej = open(os.path.join(os.getcwd(), rej), "r", newline='')
        for rej_line in open_log_rej:
            if rej_line.startswith("\trejh3_"):
                rejected_seq = (rej_line.split("|")[0]).replace("\trejh3_", "")
                if not rejected_seq in rejected:
                    rejected.append(rejected_seq)


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


def contact_count(chain_ids, latest_complex, latest_ab_complex):
    # Checks is the number of atomic contacts on the modified chain has grown too much. This would suggest that the model of the new complex is innacurate, so the modification must be discarded.
    chain = "H"
    initial_count = ""
    new_count = ""
    ag_chain_num = 1

    with redirect_stdout(open("pymol_count.pml", "w", newline="")):
        print("load " + latest_complex)
        print("select ori_" + chain + ", chain " + chain + " and name CA")
        for ag_chn in chain_ids:
            if not ag_chn == "L" and not ag_chn == "H":
                print("select ag_" + str(ag_chain_num) + ", chain " + ag_chn + " and name CA")
                print("select ag_" + str(ag_chain_num) + " near_to 4.5 of ori_" + chain)
                ag_chain_num = ag_chain_num + 1

        ag_chain_num = 1
        print("delete " + latest_complex.replace(".pdb", ""))
        print("load " + latest_ab_complex)
        print("select new_" + chain + ", chain " + chain + " and name CA")
        for ag_chn in chain_ids:
            if not ag_chn == "L" and not ag_chn == "H":
                print("select ag_" + str(ag_chain_num) + ", chain " + ag_chn + " and name CA")
                print("select ag_" + str(ag_chain_num) + " near_to 4.5 of new_" + chain)
                ag_chain_num = ag_chain_num + 1

    run_pml_cmd = str(pymol_command + " -qc pymol_count.pml")
    out = subprocess.run(run_pml_cmd, shell=True, capture_output=True)
    pml_cnt_out = out.stdout.decode("utf-8").split("\n")

    partial_old = []
    partial_new = []

    ag_chns = sum("Selector: selection \"sele\" defined with" in c for c in pml_cnt_out) / 2
    oldvnew = 1

    for pml_line in pml_cnt_out:
        if "Selector: selection \"sele\" defined with" in pml_line:
            if oldvnew <= ag_chns:
                partial_old.append(int(pml_line.split("with ")[-1].split(" atoms")[0]))
                initial_count = sum(partial_old)
            if oldvnew > ag_chns:
                partial_new.append(int(pml_line.split("with ")[-1].split(" atoms")[0]))
                new_count = sum(partial_new)
            oldvnew = oldvnew + 1

    if initial_count == 0:
        initial_count = 1

    if new_count > 2 * initial_count:
        rejected.append(tested_canon[-1])
        curr_log_rej = "rej_mod_" + curr_struct + ".log"
        updt_log_rej = open(curr_log_rej, "a", newline="")
        with redirect_stdout(updt_log_rej):
            print("\trejh3_" + tested_canon[-1] + "|Complex assembly error|FAIL")
        del tested_canon[-1]

        print("Too many new contacts; probable modelling error. Discarding modification.")
        raise Exception()


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
    chain_ids = [ch for ch in pml_out.split("cmd.get_chains:")[-1].split("\r")[0] if ch not in not_chains]

    # Refinement with OpenMM
    print("Refining complex side chains...")
    with ClearCache():
        refine(latest_ab_complex,latest_ab_complex)

    contact_ok = True

    try:
        contact_count(chain_ids, latest_complex, latest_ab_complex)
    except:
        contact_ok = False

    if contact_ok:
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
                print("\trejh3_" + tested_canon[-1] + "|Minimization error|FAIL")

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

        kT = 0
        if scoring_method == "ref15":
            kT = 1.0
        if scoring_method == "csm":
            kT = 0.593

        for z in open_curr_log:
            if z.startswith(file):
                latest_score = (z.split("|")[1]).replace("\n", "")

        approval_condition = False

        if scoring_strictness == "metro":
            approval_prob = math.exp((float(latest_score) - av_new_score) / kT)
            metropolis = random.random()
            if (float(av_new_score) < float(latest_score)) or (
                        metropolis < approval_prob and float(av_new_score) != float(latest_score)):
                approval_condition = True
            elif (float(av_new_score) > float(latest_score) and metropolis > approval_prob) or float(av_new_score) == float(
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
                print("\t" + tested_canon[-1] + "|" + str(av_new_score) + "|PASS")
                print(new_file_name + "|" + str(av_new_score))

            cmd = str("ANARCI -i " + file.replace(".anarci", "") + "_" + tested_canon[-1] + ".fasta -s martin -o " + new_file_name + " --assign_germline --use_species human -p 8")
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

            if not length_prob_mode == "restricted":
                if int(choose_length.replace("H3-", "")) >= 21:
                    long_h3.append("TEST_LONG")

            rename_cmd = str("cp " + latest_ab_complex + " complex_" + new_file_name.replace(".anarci", ".pdb"))
            rename_cmd2 = str("cp " + latest_ab.replace(".pdb", ".fasta") + " " + new_file_name.replace(".anarci", ".fasta"))
            rename_cmd3 = str("cp " + latest_ab + " " + new_file_name.replace(".anarci", ".pdb"))
            subprocess.run(rename_cmd, shell=True, capture_output=True)
            subprocess.run(rename_cmd2, shell=True, capture_output=True)
            subprocess.run(rename_cmd3, shell=True, capture_output=True)

        if not approval_condition:
            print(tested_canon[-1] + "|" + str(av_new_score) + "|Modification rejected")
            updt_log = open(curr_log, "a", newline="")
            with redirect_stdout(updt_log):
                print("\t" + tested_canon[-1] + "|" + str(av_new_score) + "|FAIL")

        # Remove scoring files from the current cycle
        for rep in glob.glob("out-min.pdb"):
            rep_chain = "CHN_" + rep
            rep_score = "SCO_" + rep
            cmd_del_int = "rm " + rep_chain + " " + rep + " " + rep_score + " out* *.out *.info *.rst7 *.crd *.leap *.pqr *.parm7 *.mdinfo *.mdout *.nc out* temp*"
            subprocess.run(cmd_del_int, shell=True, capture_output=True)


##########
# Reading configuration file to determine how many swap cycles will be attempted, the mode of choice of CDR H3 for
# grafting (random or probabilistic), and RMSD threshold above which a modelled H3 is deemed "low quality" and discarded
n_cycles = "cfg:swap_h3|n_cycles"
length_prob_mode = "cfg:swap_h3|length_prob_mode"
h3_rmsd_limit = "cfg:swap_h3|h3_rmsd_limit"
steps = "cfg:steps"
scoring_strictness = "cfg:scoring_strictness"
approval_threshold = 0
scoring_method = "cfg:scoring_method"
ph = "7"
pymol_command = "pymol"

cfg = open("swap_settings.cfg").readlines()
for setting in cfg:
    if "swap_h3|n_cycles" in setting:
        n_cycles = int(setting.split("=")[1].replace("\n", ""))
    if "swap_h3|length_prob_mode" in setting:
        length_prob_mode = str((setting.split("=")[1].replace("\n", "")))
    if "swap_h3|h3_rmsd_limit" in setting:
        h3_rmsd_limit = float(setting.split("=")[1].replace("\n", ""))
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

print("Starting new step: CDR H3 grafting")

if scoring_method == "ref15":
    from pyrosetta import *
    init("-mute all")

base_structure = []
# Stores which and how many antibodies are being developed simultaneously (numbered 001, 002...)
structures = []
distribution = {}
seqs_by_size = {}
long_distribution = {}
read_lengths = open("ref_seq/h3_length_distribution.csv", "r").readlines()
read_long_lengths = open("ref_seq/long_h3_length_distribution.csv", "r").readlines()

step_list = steps.split("|")
if "" in step_list:
    step_list.remove("")
bs_index = step_list.index("h3")
base_structure.append(str(bs_index + 1))
base_structure.append(str(bs_index + 2))

for struc in glob.glob("bs" + base_structure[0] + "_*"):
    structure = struc.split("_")[1]
    if not structure in structures:
        structures.append(structure)

# Stores the probabilities of each H3 length being chosen for grafting; probabilities depend on the chosen modification mode (random or probabilistic)
for length in read_lengths:
    size = (length.split(";")[0].replace("ï»¿", "")).replace("\ufeff", "")
    frequency = int(length.split(";")[1].replace("\n", ""))
    seqs_by_size[size] = frequency
    if length_prob_mode == "probabilistic":
        distribution[size] = frequency
    if length_prob_mode == "random":
        distribution[size] = 1

for long_length in read_long_lengths:
    size = (long_length.split(";")[0].replace("ï»¿", "")).replace("\ufeff", "")
    frequency = int(long_length.split(";")[1].replace("\n", ""))
    if length_prob_mode == "probabilistic":
        long_distribution[size] = frequency
    if length_prob_mode == "random":
        long_distribution[size] = 1

for st in structures:
    # Stores which modifications have already been tested
    tested_canon = []
    # Stores the IDs of sequences that were rejected for producing models with low H3 confidence
    rejected = []
    # Stores the name of the final structure, after all modifications, and its score
    final = []
    final_sco = []
    # Checks if there have been any errors during the modification process
    error = []
    # If a long H3 CDR gets grafted and approved, the 10 following cycles will be used to test other long sequences
    long_h3 = []
    initial_h3_len = "size"


    if length_prob_mode == "restricted":
        for first_h3_log in glob.glob("swap_" + st + ".log"):
            open_first_h3_log = open(os.path.join(os.getcwd(), first_h3_log), "r", newline='')
            all_h3 = []
            for line in open_first_h3_log:
                if line.startswith("bs" + base_structure[0]):
                    name = line.split("|")[0]
                    all_h3.append(name)
            open_first_h3 = open(all_h3[0]).readlines()
            initial_h3_len = open_first_h3[0].split("|")[4]
        lowest_size = int(initial_h3_len) - 1
        largest_size = int(initial_h3_len) + 1
        for i in range(lowest_size, (largest_size + 1)):
            distribution["H3-" + str(i)] = 1

    print("New structure")

    while len(tested_canon) < n_cycles:
        latest = []
        num_list = []
        run_list = []
        sco_list = []
        condition = []
        CDRH3_res = [95, 102]

        if not long_h3:
            # Finds the antibody structure with the best score in the main log file, which will serve as the starting point for the current modification cycle
            read_log_init()

            for new in latest:
                for file in glob.glob(new):
                    # Opens anarci and log files
                    open_num = open(file, "r").readlines()
                    curr_struct = file.split("_")[1]
                    log = (open("swap_" + curr_struct + ".log", "r")).readlines()
                    id = file.split("_")[2]
                    # Chooses a CDR H3 length for the new grafted sequence, based on the specified swap mode
                    choose_length = random.choices(list(distribution.keys()), weights=list(distribution.values()), k=1)[0]
                    print(choose_length)

                    canon_numbers = int(seqs_by_size[choose_length])/2
                    # The pieces that will form the modified antibody sequence are stored here, along with the length of the CDR that is being grafted
                    start = []
                    new_CDR = []
                    new_CDR_len = []
                    end = []
                    # Header information for each chain is stored here
                    cdr_header = []
                    main_header = []
                    other_chain = []
                    print("Starting new cycle")
                    print("Modifying CDR H3...")

                    number_tested = sum(s.count(choose_length) for s in tested_canon) + sum(r.count(choose_length) for r in rejected)
                    canon_file = os.path.dirname(
                        os.path.abspath(__file__)) + "/h3_bank/" + choose_length + "_cdr.fasta"
                    open_canon = (open(canon_file, "r")).readlines()
                    # Chooses a random CDR H3 sequence among those with the selected length in the database
                    make_new_seq()

                    if main_header:
                        # If the sequence was sucessfully modified, the length of the new CDR will be checked and the modified sequence will be modeled
                        new_header_heavy = []
                        new_header_light = []
                        mod_header()
                        print("Modelling new antibody...")

                        try:
                            model_ab()
                        except:
                            print("Modelling error. Restarting cycle...")
                            error.append("ERROR")
                            # If an error occurs during modelling, that step will be repeated
                            rejected.append(tested_canon[-1])
                            curr_log_rej = "rej_mod_" + curr_struct + ".log"
                            updt_log_rej = open(curr_log_rej, "a", newline="")
                            with redirect_stdout(updt_log_rej):
                                print("\trejh3_" + tested_canon[-1] + "|Modelling error|FAIL")
                            del tested_canon[-1]

                            continue

                        if not "ERROR" in error:
                            latest_ab = file.replace(".anarci", "") + "_" + tested_canon[-1] + ".pdb"
                            latest_complex = "complex_" + file.replace(".anarci", ".pdb")
                            latest_ab_complex = "complex_" + latest_ab
                            temporary = "temp"
                            renum = "_renum"

                            qual_ctrl(latest_ab, latest_complex, latest_ab_complex, temporary, renum, scoring_strictness, approval_threshold)

                    if "PASS" in condition:
                        break
        # If a long H3 CDR gets grafted and approved, the 10 following cycles will be used to test other long sequences
        if "TEST_LONG" in long_h3:
            for i in range(1,11):
                print("Long H3 approved. Testing other long H3 sequences. Cycle number = " + str(i))
                # Finds the antibody structure with the best score in the main log file, which will serve as the starting point for the current modification cycle
                read_log_init()

                for new in latest:
                    for file in glob.glob(new):
                        open_num = open(file, "r").readlines()
                        curr_struct = file.split("_")[1]
                        log = (open("swap_" + curr_struct + ".log", "r")).readlines()
                        id = file.split("_")[2]

                        choose_length = random.choices(list(long_distribution.keys()), weights=list(long_distribution.values()), k=1)[
                            0]
                        print(choose_length)
                        canon_numbers = int(long_distribution[choose_length]) / 2

                        start = []
                        new_CDR = []
                        new_CDR_len = []
                        end = []

                        cdr_header = []
                        main_header = []
                        other_chain = []
                        print("Starting new cycle")
                        print("Modifying CDR H3...")

                        number_tested = sum((choose_length) in s for s in tested_canon)
                        canon_file = os.path.dirname(
                            os.path.abspath(__file__)) + "/h3_bank/" + choose_length + "_cdr.fasta"
                        open_canon = (open(canon_file, "r")).readlines()

                        make_new_seq()

                        if main_header:
                            new_header_heavy = []
                            new_header_light = []
                            mod_header()
                            print("Modelling new antibody...")

                            try:
                                model_ab()
                            except:
                                print("Modelling error. Restarting cycle...")
                                error.append("ERROR")
                                # If an error occurs during modelling, that step will be repeated
                                rejected.append(tested_canon[-1])
                                curr_log_rej = "rej_mod_" + curr_struct + ".log"
                                updt_log_rej = open(curr_log_rej, "a", newline="")
                                with redirect_stdout(updt_log_rej):
                                    print("\trejh3_" + tested_canon[-1] + "|Modelling error|FAIL")
                                del tested_canon[-1]

                                continue

                            if not "ERROR" in error:
                                latest_ab = file.replace(".anarci", "") + "_" + tested_canon[-1] + ".pdb"
                                latest_complex = "complex_" + file.replace(".anarci", ".pdb")
                                latest_ab_complex = "complex_" + latest_ab
                                temporary = "temp"
                                renum = "_renum"

                                qual_ctrl(latest_ab, latest_complex, latest_ab_complex, temporary, renum, scoring_strictness, approval_threshold)

                        if "PASS" in condition:
                            break
            print("End of long H3 testing cycles")
            long_h3.clear()

    if not "ERROR" in error:
        # At the end of all modification cycles, the structure and sequence of the best antibody produced are copied and renamed, in preparation for the next modification step
        print("Finalizing current structure...")
        final_num = ""
        if not step_list[-1] == "h3":
            final_num = str(final[-1]).replace("bs" + base_structure[0], "bs" + base_structure[1])
        if step_list[-1] == "h3":
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
