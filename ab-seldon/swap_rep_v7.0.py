import glob
import os
import subprocess
import random
import math
import time
import re
import torch
import ast
import numpy as np
import csv
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


def calculate_DSSP(pdb_file):
    """Run DSSP analysis and return list of dictionaries with residue secondary structures"""
    filename = os.path.splitext(os.path.basename(pdb_file))[0]
    dssp_file = f'{filename}.dssp'

    try:
        subprocess.run(
            f"mkdssp {pdb_file} --output-format dssp > {dssp_file}",
            shell=True,
            check=True,
            stderr=subprocess.DEVNULL
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error running DSSP on {pdb_file}") from e

    dssp_data = []
    with open(dssp_file, 'r') as f:
        lines = f.readlines()
        start_idx = None
        for i, line in enumerate(lines):
            if 'RESIDUE AA' in line:
                start_idx = i + 1
                break

        if start_idx is None:
            raise RuntimeError("Could not find RESIDUE AA header in DSSP file")

        for line in lines[start_idx:]:
            if len(line) < 20:
                continue
            dssp_data.append({
                'residue_num': line[5:10].strip(),
                'chain': line[11].strip(),
                'secondary_structure': line[16].strip() if line[16].strip() else ' '
            })

    if os.path.exists(dssp_file):
        os.remove(dssp_file)

    return dssp_data


def get_ss_description(code):
    """Convert DSSP single-letter code to human-readable description"""
    ss_map = {
        'H': 'Alpha helix', 'B': 'Beta bridge', 'E': 'Strand',
        'G': 'Helix_3', 'I': 'Helix_5', 'P': 'Helix_PPII',
        'T': 'Turn', 'S': 'Bend', ' ': 'Loop'
    }
    return ss_map.get(code, 'Unknown')


def get_structure_category(ss_code):
    """Classify DSSP code into low_order, beta, or helix category"""
    structure_types = {
        'low_order': ['T', 'S', ' '],
        'beta': ['B', 'E'],
        'helix': ['H', 'G', 'I', 'P']
    }
    for category, codes in structure_types.items():
        if ss_code in codes:
            return category
    return 'unknown'


def is_change_acceptable(ref_ss, comp_ss):
    """Determine if secondary structure change meets acceptability criteria"""
    ref_cat = get_structure_category(ref_ss)
    comp_cat = get_structure_category(comp_ss)

    if ref_cat == comp_cat:
        return True
    if ref_cat == 'low_order' and comp_cat in ['beta', 'helix']:
        return True

    # Note: B -> low_order is technically a loss of structure, so it returns False here.
    # The tolerance for the "first" B->low change is handled in compare_epitope_DSSP.
    return False


def compare_epitope_DSSP(reference_csv, comparison_pdb, output_csv=None, tolerance=3):
    """Compare epitope secondary structures between reference CSV and comparison PDB"""
    if not os.path.exists(reference_csv):
        raise FileNotFoundError(f"Reference CSV not found: {reference_csv}")

    # Read reference CSV
    ref_data = []
    with open(reference_csv, 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        columns = reader.fieldnames

        if 'chain/resnumber' not in columns or 'DSSP_code' not in columns:
            raise ValueError("CSV must contain 'chain/resnumber' and 'DSSP_code' columns")

        for row in reader:
            ref_data.append(row)

    print("Checking epitope secondary structure with DSSP...")
    dssp_compare = calculate_DSSP(comparison_pdb)

    if len(dssp_compare) == 0:
        raise RuntimeError("DSSP analysis produced no data for comparison structure")

    changes = []
    b_to_low_changes = []

    for row in ref_data:
        residue = row['chain/resnumber']
        ref_ss = row['DSSP_code']
        chain, res_num = residue.split('/')

        # Find matching entry in dssp_compare
        comp_ss = None
        for entry in dssp_compare:
            if entry['chain'] == chain and entry['residue_num'] == res_num:
                comp_ss = entry['secondary_structure']
                break

        if comp_ss is None:
            continue

        if ref_ss != comp_ss:
            changes.append({
                'residue': residue,
                'reference_DSSP': ref_ss,
                'comparison_DSSP': comp_ss,
                'reference_description': get_ss_description(ref_ss),
                'comparison_description': get_ss_description(comp_ss),
                'reference_category': get_structure_category(ref_ss),
                'comparison_category': get_structure_category(comp_ss)
            })

            if ref_ss == 'B' and get_structure_category(comp_ss) == 'low_order':
                b_to_low_changes.append(residue)

    b_to_low_count = len(b_to_low_changes)

    for change in changes:
        change['acceptable'] = is_change_acceptable(
            change['reference_DSSP'],
            change['comparison_DSSP']
        )

    if len(changes) > 0:
        # Calculate raw penalty count (all unacceptable changes)
        raw_unacceptable_count = sum(1 for change in changes if not change['acceptable'])

        # Calculate adjustment: The first B->Low change is "free" (subtract 1 from penalty if it exists)
        adjustment = 1 if b_to_low_count > 0 else 0
        final_penalty_count = raw_unacceptable_count - adjustment

        if final_penalty_count < tolerance:
            print(f"DSSP Passed (Penalty Score: {final_penalty_count})")
        else:
            # Export CSV only when there are unacceptable changes exceeding tolerance
            if output_csv is not None:
                with open(output_csv, 'w', newline='') as csvfile:
                    fieldnames = ['residue', 'reference_DSSP', 'comparison_DSSP',
                                  'reference_description', 'comparison_description',
                                  'reference_category', 'comparison_category', 'acceptable']
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    writer.writerows(changes)

            print(
                f"Loss of epitope secondary structure detected. Penalty Score: {final_penalty_count} (Threshold: {tolerance})")
            raise RuntimeError(f"Loss of epitope secondary structure detected.")
    else:
        print("DSSP Passed")

    return changes


def mod_header():
    # Checks the length of the CDRs and stores them on a modified header
    lengths = {"H1_len": (main_header[0]).split("|")[2], "H2_len": (main_header[0]).split("|")[3],
               "H3_len": (main_header[0]).split("|")[4], "L1_len": (main_header[1]).split("|")[2],
               "L2_len": (main_header[1]).split("|")[3], "L3_len": (main_header[1]).split("|")[4]}

    headers = {"start_H": str("|".join([str(item) for item in [m for m in (main_header[0]).split("|")[:2]]])),
               "start_L": str("|".join([str(item) for item in [m for m in (main_header[1]).split("|")[:2]]])),
               "end_H": str("|".join([str(item) for item in [g for g in (main_header[0]).split("|")[5:]]])),
               "end_L": str("|".join([str(item) for item in [g for g in (main_header[1]).split("|")[5:]]]))}

    lengths[str(chain + x) + "_len"] = str(new_CDR_len)

    new_header_h = headers.get("start_H") + "|" + lengths.get("H1_len") + "|" + lengths.get(
        "H2_len") + "|" + lengths.get(
        "H3_len") + "|" + headers.get("end_H")
    new_header_l = headers.get("start_L") + "|" + lengths.get("L1_len") + "|" + lengths.get(
        "L2_len") + "|" + lengths.get(
        "L3_len") + "|" + headers.get("end_L")

    new_header_light.append(new_header_l)
    new_header_heavy.append(new_header_h)


def model_ab():
    # Models the modified antibody structure
    sequences = {}

    if chain == "H":
        sequences["H"] = newseq
        sequences["L"] = "".join([str(item) for item in [b for b in other_chain if not b == "-"]])
    if chain == "L":
        sequences["H"] = "".join([str(item) for item in [b for b in other_chain if not b == "-"]])
        sequences["L"] = newseq

    output_file = file.replace(".anarci", "") + "_" + tested_canon[-1] + ".pdb"

    with ClearCache():
        antibody = predictor.predict(sequences)
        antibody.save(output_file, check_for_strained_bonds=False)

    make_fasta = output_file.replace(".pdb", ".fasta")
    if chain == "H":
        with redirect_stdout(open(make_fasta, "w", newline="")):
            print(">" + new_header_heavy[0])
            print(newseq)
            print(">" + new_header_light[0])
            print("".join([str(item) for item in [b for b in other_chain if not b == "-"]]))
    if chain == "L":
        with redirect_stdout(open(make_fasta, "w", newline="")):
            print(">" + new_header_heavy[0])
            print("".join([str(item) for item in [b for b in other_chain if not b == "-"]]))
            print(">" + new_header_light[0])
            print(newseq)


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
                        new_line = line[:21].replace("LYN", "LYS").replace("HETATM",
                                                                           "ATOM  ") + current_chain_id + line[
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


def calculate_translation_vector(paratope_com, epitope_com, distance=0.5):
    # Calculates the translation vector to move the antibody away from the antigen.
    ab_com = np.array(paratope_com)
    ag_com = np.array(epitope_com)
    direction_vector = ab_com - ag_com

    magnitude = np.linalg.norm(direction_vector)
    unit_vector = direction_vector / magnitude
    translation_vector = unit_vector * distance

    return translation_vector


def contact_count(chain_ids, chain, latest_ab_complex):
    # Checks if the structural alignment has produced a complex with unrealistically entangled chains. If so, the modification is discarded.
    ab_chain = chain
    no_ext = latest_ab_complex.replace(".pdb", "")
    with redirect_stdout(open("pymol_com.pml", "w", newline="")):
        print(f"load {latest_ab_complex}")
        for ag_chn in chain_ids:
            if not ag_chn == "L" and not ag_chn == "H":
                print(f"select epitope, chain {ag_chn} near_to 10.0 of chain {ab_chain}")
                print(f"select paratope, chain {ab_chain} near_to 10.0 of chain {ag_chn}")
        print("centerofmass epitope")
        print("centerofmass paratope")

    run_com_cmd = str(pymol_command + " -qc pymol_com.pml")
    out = subprocess.run(run_com_cmd, shell=True, capture_output=True)
    pml_com_out = out.stdout.decode("utf-8").split("\n")

    epitope_com = tuple(ast.literal_eval(pml_com_out[-4].replace("Center of Mass: ", "")))
    paratope_com = tuple(ast.literal_eval(pml_com_out[-2].replace("Center of Mass: ", "")))

    translation = calculate_translation_vector(paratope_com, epitope_com, 0.5)
    distance_around = 2.7

    with redirect_stdout(open("pymol_count.pml", "w", newline="")):
        print("load " + latest_ab_complex)
        for i in range(15):
            for ag_chn in chain_ids:
                if not ag_chn == "L" and not ag_chn == "H":
                    print(
                        f"select chain {ag_chn} and name CA near_to {distance_around} of chain {ab_chain} and name CA")
            print(
                f"translate [{translation[0]:.6f},{translation[1]:.6f},{translation[2]:.6f}], {no_ext} and chain {ab_chain}")

    run_pml_cmd = str(pymol_command + " -qc pymol_count.pml")
    out = subprocess.run(run_pml_cmd, shell=True, capture_output=True)
    pml_cnt_out = out.stdout.decode("utf-8").split("\n")

    counts = []

    for pml_line in pml_cnt_out:
        if "Selector: selection \"sele\" defined with" in pml_line:
            contacts = pml_line.split("Selector")[1:]
            for cn in contacts:
                counts.append(int(cn.split("defined with ")[-1].split(" atoms")[0]))

    if any(c > 0 for c in counts):
        rejected.append(tested_canon[-1])
        curr_log_rej = "rej_mod_" + curr_struct + ".log"
        updt_log_rej = open(curr_log_rej, "a", newline="")
        with redirect_stdout(updt_log_rej):
            print("\trejrep_" + tested_canon[-1] + "|Chain entanglement|FAIL")
        del tested_canon[-1]

        raise Exception("Inter-chain entanglement detected. Discarding modification.")


def new_complex(latest_ab, latest_complex, latest_ab_complex, temporary, renum, scoring_strictness, approval_threshold,
                chain, dssp_tolerance):
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
    chain_ids = [ch for ch in pml_out.split("cmd.get_chains:")[-1].split("\r")[0].split("\n")[0] if ch not in not_chains]

    # Refinement with OpenMM
    print("Refining complex side chains...")
    with ClearCache():
        refine(latest_ab_complex, latest_ab_complex)

    contact_ok = True
    DSSP_ref_ok = True

    try:
        contact_count(chain_ids, chain, latest_ab_complex)
    except:
        print("Inter-chain entanglement detected. Discarding modification.")
        contact_ok = False

    # DSSP check after refinement
    reference_csv = latest_complex.split("_")[3] + "_" + latest_complex.split("_")[4] + "_epitope-DSSP.csv"
    try:
        compare_epitope_DSSP(reference_csv, latest_ab_complex,
                             "DSSP_rep_refinement_" + latest_ab_complex.replace(".pdb", ".csv"), dssp_tolerance)
    except:
        print("Discarding modification. Check DSSP_rep_refinement_" + latest_ab_complex.replace(".pdb",
                                                                                                ".csv") + " for more information.")
        DSSP_ref_ok = False

    if contact_ok and DSSP_ref_ok:
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

        except:
            print("Minimization failure. Check " + tested_canon[-1] + "_error.tar for more information.")
            rejected.append(tested_canon[-1])
            curr_log_rej = "rej_mod_" + curr_struct + ".log"
            updt_log_rej = open(curr_log_rej, "a", newline="")
            with redirect_stdout(updt_log_rej):
                print("\trejrep_" + tested_canon[-1] + "|Minimization error|FAIL")

            for rep in glob.glob("out-min.pdb"):
                rep_chain = "CHN_" + rep
                rep_score = "SCO_" + rep
                cmd_del_int = "rm " + rep_chain + " " + rep + " " + rep_score + " out* *.out *.info *.rst7 *.crd *.leap *.pqr *.parm7 *.mdinfo *.mdout *.nc temp*"
                cmd_bkp_error = "tar -cvf " + tested_canon[
                    -1] + "_error.tar " + rep_chain + " " + rep + " " + rep_score + " out* *.out *.info *.rst7 *.crd *.leap *.pqr *.parm7 *.mdinfo *.mdout *.nc temp*"
                subprocess.run(cmd_bkp_error, shell=True, capture_output=True)
                subprocess.run(cmd_del_int, shell=True, capture_output=True)

            del tested_canon[-1]

        # DSSP check after minimization
        DSSP_min_ok = True
        try:
            compare_epitope_DSSP(reference_csv, "SCO_out-min.pdb",
                                 "DSSP_rep_minimization_" + latest_ab_complex.replace(".pdb", ".csv"), dssp_tolerance)
        except:
            print("Discarding modification. Check DSSP_rep_minimization_" + latest_ab_complex.replace(".pdb",
                                                                                                      ".csv") + " for more information.")
            DSSP_min_ok = False

        if DSSP_min_ok:
            scoring(latest_ab, latest_ab_complex, scoring_strictness, approval_threshold)
        else:
            rejected.append(tested_canon[-1])
            curr_log_rej = "rej_mod_" + curr_struct + ".log"
            updt_log_rej = open(curr_log_rej, "a", newline="")
            with redirect_stdout(updt_log_rej):
                print("\trejrep_" + tested_canon[-1] + "|Post-minimization error|FAIL")

            for rep in glob.glob("out-min.pdb"):
                rep_chain = "CHN_" + rep
                rep_score = "SCO_" + rep
                cmd_del_int = "rm " + rep_chain + " " + rep + " " + rep_score + " out* *.out *.info *.rst7 *.crd *.leap *.pqr *.parm7 *.mdinfo *.mdout *.nc temp*"
                subprocess.run(cmd_del_int, shell=True, capture_output=True)

            del tested_canon[-1]

    else:
        print("Modification discarded during structure refinement.")
        rejected.append(tested_canon[-1])
        curr_log_rej = "rej_mod_" + curr_struct + ".log"
        updt_log_rej = open(curr_log_rej, "a", newline="")
        with redirect_stdout(updt_log_rej):
            print("\trejrep_" + tested_canon[-1] + "|Complex assembly error|FAIL")

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
                print("\t" + tested_canon[-1] + "|" + str(av_new_score) + "|PASS")
                print(new_file_name + "|" + str(av_new_score))

            cmd = str("ANARCI -i " + file.replace(".anarci", "") + "_" + tested_canon[
                -1] + ".fasta -s martin -o " + new_file_name + " --assign_germline -p 8 --use_species human")
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
                print("\t" + tested_canon[-1] + "|" + str(float(av_new_score)) + "|FAIL")

        # Removes scoring files from the current cycle
        for rep in glob.glob("out-min.pdb"):
            rep_chain = "CHN_" + rep
            rep_score = "SCO_" + rep
            cmd_del_int = "rm " + rep_chain + " " + rep + " " + rep_score + " out* *.out *.info *.rst7 *.crd *.leap *.pqr *.parm7 *.mdinfo *.mdout *.nc out* temp*"
            subprocess.run(cmd_del_int, shell=True, capture_output=True)


##########
# Reading configuration file to determine the maximum number of swap cycles allowed; and the
# probabilities of a CDR or chain being chosen for modification
max_cycles = "cfg:swap_rep|max_cycles"
cdr_prob = "cfg:swap_rep|cdr_prob"
chain_prob = "cfg:swap_rep|chain_prob"
steps = "cfg:steps"
scoring_strictness = "cfg:scoring_strictness"
approval_threshold = 0
scoring_method = "cfg:scoring_method"
ph = "7"
pymol_command = "pymol"
dssp_tolerance = 3

cfg = open("swap_settings.cfg").readlines()
for setting in cfg:
    if "swap_rep|max_cycles" in setting:
        max_cycles = int(str(setting.split("=")[1].replace("\n", "")))
    if "swap_rep|cdr_prob" in setting:
        cdr_prob = str(setting.split("=")[1].replace("\n", ""))
    if "swap_rep|chain_prob" in setting:
        chain_prob = str(setting.split("=")[1].replace("\n", ""))
    if "minimization_pH" in setting:
        ph = str(setting.split("=")[1].replace("\n", ""))
    if "pymol_command" in setting:
        pymol_command = str(setting.split("=")[1].replace("\n", ""))
    if "scoring_method" in setting:
        scoring_method = str(setting.split("=")[1].replace("\n", ""))
    if "steps" in setting:
        steps = str(setting.split("=")[1].replace("\n", "").replace("\r", ""))
    if "scoring_strictness" in setting:
        scoring_strictness = str(setting.split("=")[1].replace("\n", "").replace("\r", ""))
    if "approval_threshold" in setting:
        if scoring_strictness == "strict":
            approval_threshold = float(setting.split("=")[1].replace("\n", "").replace("\r", ""))
    if "DSSP_changes_tolerance" in setting:
        dssp_tolerance = int(setting.split("=")[1].replace("\n", "").replace("\r", ""))

print("Starting new step: Representative non-H3 CDR grafting")

if scoring_method == "ref15":
    from pyrosetta import *

    init("-mute all")

base_structure = []
# Stores which and how many antibodies are being developed simultaneously (numbered 001, 002...)
structures = []

step_list = steps.split("|")
if "" in step_list:
    step_list.remove("")
bs_index = step_list.index("rep")
base_structure.append(str(bs_index + 1))
base_structure.append(str(bs_index + 2))

for struc in glob.glob("bs" + base_structure[0] + "_*"):
    structure = struc.split("_")[1]
    if not structure in structures:
        structures.append(structure)

for st in structures:
    # Stores which modifications have already been tested
    tested_canon = []
    # Stores the IDs of sequences that were rejected for producing erroneous models
    rejected = []
    # Stores the name of the final structure, after all modifications, and its score
    final = []
    final_sco = []
    # Checks if there have been any errors during the modification process
    error = False
    print("New structure")

    # Keeps track of the current modification supercycle and of how many consensus sequences have been tested during it
    supercycle = 1
    curr_sc_tested = 0

    # Main loop; this "while" loop covers the entire modification/scoring cycle.
    # This loop will be repeated 53 times, to test all possible canonical structures for each CDR (except H3)
    while curr_sc_tested < 53 and len(tested_canon) < max_cycles:
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
                if line.startswith("bs" + base_structure[0]):
                    name = (line.split("|")[0]).replace("\n", "")
                    score = (line.split("|")[1]).replace("\n", "")
                    progress = (name.split("_")[-1]).replace(".anarci", "")
                    run = name.split(progress)[0]
                    run_list.append(run)
                    num_list.append(progress)
                    sco_list.append(float(score))
                if line.startswith("\tsc"):
                    tested = (line.split("|")[0]).replace("\t", "")
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
                if rej_line.startswith("\trejrep_"):
                    rejected_seq = (rej_line.split("|")[0]).replace("\trejrep_", "")
                    if not rejected_seq in rejected:
                        rejected.append(rejected_seq)

        curr_sc_tested = sum(s.count('sc' + str(supercycle)) for s in tested_canon) + sum(
            s.count('sc' + str(supercycle)) for s in rejected)

        for new in latest:
            for file in glob.glob(new):
                # Opens anarci and log files
                open_num = open(file, "r").readlines()
                curr_struct = file.split("_")[1]
                log = (open("swap_" + curr_struct + ".log", "r")).readlines()
                id = file.split("_")[2]
                # List of chains and CDR numbers. H_list does not include CDR H3.
                Ch_list = ["H", "L"]
                H_list = ["1", "2"]
                L_list = ["1", "2", "3"]
                # Number of canonical structures for each CDR (Kelow,2022)
                canon_numbers = {"H1": 8, "H2": 8, "L1": 17, "L2": 3, "L3": 17}
                # Residues that define each CDR and DE loop, according to the Martin/Enhanced Chothia scheme. H1 and H2 altered to fit the consensus sequences.
                CDRH1_res = [23, 35]
                CDRH2_res = [50, 58]
                CDRL1_res = [24, 34]
                CDRL2_res = [50, 56]
                CDRL3_res = [89, 97]
                H4_res = [71, 75]
                L4_res = [66, 71]
                # Chain being modified is chosen according to the specified probabilities from the cfg file
                swap_Ch = []
                if curr_sc_tested < 40:
                    swap_Ch = random.choices(Ch_list,
                                             weights=(int(chain_prob.split("|")[0]), int(chain_prob.split("|")[1])),
                                             k=1)
                if curr_sc_tested >= 40:
                    print("Most canonical structures have been tested. Selecting a random CDR.")
                    swap_Ch = random.choices(Ch_list, weights=(1, 1), k=1)
                # The pieces that will form the modified antibody sequence are stored here, along with the length of the CDR that is being grafted
                start = []
                new_CDR = "Sequence"
                new_CDR_len = 0
                mid = []
                new_DE = "Loop DE"
                end = []
                # Header information for each chain is stored here
                cdr_header = []
                main_header = []
                other_chain = []
                print("Starting new cycle")
                print("Modifying CDR...")

                for chain in swap_Ch:
                    # CDR being modified is chosen according to the specified probabilities in the cfg file
                    swapH = []
                    swapL = []

                    if curr_sc_tested < 40:
                        swapH = random.choices(H_list,
                                               weights=(int(cdr_prob.split("|")[0]), int(cdr_prob.split("|")[1])), k=1)
                        swapL = random.choices(L_list, weights=(
                            int(cdr_prob.split("|")[2]), int(cdr_prob.split("|")[3]), int(cdr_prob.split("|")[4])), k=1)
                    if curr_sc_tested >= 40:
                        swapH = random.choices(H_list, weights=(1, 1), k=1)
                        swapL = random.choices(L_list, weights=(1, 1, 1), k=1)

                    DE_co_graft = 0
                    if chain == "H":
                        DE_co_graft = 2
                    if chain == "L":
                        DE_co_graft = 1

                    for x in locals()["swap" + chain]:
                        # Checks if all canonical structures of the chosen CDR have already been tested or not
                        number_tested = sum(
                            s.count('sc' + str(supercycle) + "_" + chain + x) for s in tested_canon) + sum(
                            r.count('sc' + str(supercycle) + "_" + chain + x) for r in rejected)

                        if number_tested < canon_numbers[str(chain + x)]:
                            while not cdr_header:
                                random_cdr = random.choice(list(open("ref_seq/pyigclassify2_consensus_v4.txt")))
                                if random_cdr.startswith(str(chain + x)):
                                    # Randomly chooses one of the consensus sequences for the chosen CDR to graft into the antibody
                                    canon = random_cdr.split("|")[0]
                                    rep = (random_cdr.split("|")[1])
                                    loop_DE = ((random_cdr.split("|")[2]).replace("\n", "")).replace("\t", "")
                                    sc_canon = "sc" + str(supercycle) + "_" + canon
                                    if sc_canon not in tested_canon and sc_canon not in rejected:
                                        cdr_header.append(canon)
                                        new_CDR = rep
                                        new_CDR_len = int((random_cdr.split("|")[4]))
                                        new_DE = loop_DE
                                        tested_canon.append(sc_canon)
                                    else:
                                        pass

                            for d in open_num:
                                other_chain_id = [c for c in Ch_list if not c == chain][-1]
                                # Uses the anarci file to determine where in the original sequence the new CDR sequence should be introduced
                                if d.startswith("# " + id):
                                    ch_header = (d.split("# ")[1]).replace("\n", "")
                                    main_header.append(ch_header)

                                if d.startswith(chain):
                                    num_res = d.split(" ")[1]
                                    amino = (d.split(" ")[-1]).replace("\n", "")

                                    # If the CDR being subjected to the grafting process is H2, the DE loop will be co-grafted.
                                    # The same applies to CDR L1, if that CDR is chosen for modification.
                                    # Based on conclusions from Kelow,2020
                                    if not int(x) == DE_co_graft:
                                        if int(num_res) < (locals()["CDR" + chain + x + "_res"])[0]:
                                            start.append(amino)

                                        if int(num_res) > (locals()["CDR" + chain + x + "_res"])[-1]:
                                            end.append(amino)

                                    if int(x) == DE_co_graft:
                                        if int(num_res) < (locals()["CDR" + chain + x + "_res"])[0]:
                                            start.append(amino)

                                        if int(num_res) > locals()["CDR" + chain + x + "_res"][-1] and int(
                                                num_res) < locals()[chain + "4_res"][0]:
                                            mid.append(amino)

                                        if int(num_res) > locals()[chain + "4_res"][-1]:
                                            end.append(amino)

                                if d.startswith(other_chain_id):
                                    amino = (d.split(" ")[-1]).replace("\n", "")
                                    other_chain.append(amino)

                        else:
                            print(
                                "All consensus sequences for CDR " + chain + x + " have been tested. Trying another CDR...")

                        if main_header:
                            # If the sequence was sucessfully modified, the length of the new CDR will be checked and the modified sequence will be modeled
                            new_header_heavy = []
                            new_header_light = []

                            newseq = "sequence"

                            if int(x) == DE_co_graft:
                                newseq = "".join(
                                    item for item in [w for w in start if not w == "-"]) + new_CDR + "".join(
                                    item for item in [w for w in mid if not w == "-"]) + new_DE + "".join(
                                    item for item in [w for w in end if not w == "-"])

                            if not int(x) == DE_co_graft:
                                newseq = "".join(
                                    item for item in [w for w in start if not w == "-"]) + new_CDR + "".join(
                                    item for item in [w for w in end if not w == "-"])

                            mod_header()
                            print("Modelling new antibody...")

                            try:
                                model_ab()
                            except:
                                # If an error occurs during modelling, that step will be repeated
                                rejected.append(tested_canon[-1])
                                curr_log_rej = "rej_mod_" + curr_struct + ".log"
                                updt_log_rej = open(curr_log_rej, "a", newline="")
                                with redirect_stdout(updt_log_rej):
                                    print("\trejrep_" + tested_canon[-1] + "|Modelling error|FAIL")
                                del tested_canon[-1]

                                print("Modelling error. Restarting cycle...")
                                error = True
                                continue

                            if not error:
                                # If the modelling step gets completed sucessfully, the new antibody structure will be scored according to its affinity with the target protein
                                latest_ab = file.replace(".anarci", "") + "_" + tested_canon[-1] + ".pdb"
                                latest_complex = "complex_" + file.replace(".anarci", ".pdb")
                                latest_ab_complex = "complex_" + latest_ab
                                temporary = "temp"
                                renum = "_renum"

                                new_complex(latest_ab, latest_complex, latest_ab_complex, temporary, renum,
                                            scoring_strictness, approval_threshold, chain, dssp_tolerance)

                if "PASS" in condition:
                    break

    if not error:
        # At the end of all modification cycles, the structure and sequence of the best antibody produced are copied and renamed, in preparation for the next modification step
        print("Finalizing current structure...")
        final_num = ""
        if not step_list[-1] == "rep":
            final_num = str(final[-1]).replace("bs" + base_structure[0], "bs" + base_structure[1])
        if step_list[-1] == "rep":
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