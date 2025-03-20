
import glob
import os
import subprocess
import time
import re
import random
import math
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


def model_hyb_h():
    # Model the memory/naive hybrid structures
    print("Modelling hybrid (memory H chain) antibodies...")
    counter = ("%04i" % i for i in count(1))
    for key, value in hybrids_h.items():
        number = next(counter)
        sequences = {"H": str(key.split("-|")[1]), "L": str(value.split("-|")[1])}
        output_file = (fasta.replace(fasta.split("_")[-1], "") + "hyb-H_" + number + ".pdb").replace("bs" + base_structure[0] + "_", "")

        try:
            with ClearCache():
                antibody = predictor.predict(sequences)
                antibody.save(output_file, check_for_strained_bonds=False)
        except:
            print("Modelling failure in hybrid stage (H). Trying again.")
            with ClearCache():
                antibody = predictor.predict(sequences)
                antibody.save(output_file, check_for_strained_bonds=False)

        make_fasta = output_file.replace(".pdb", ".fasta")
        with redirect_stdout(open(make_fasta, "w", newline="")):
            print(">" + str(key.split("-|")[0]))
            print(str(key.split("-|")[1]))
            print(">" + str(value.split("-|")[0]))
            print(str(value.split("-|")[1]))


def model_hyb_l():
    # Model the memory/naive hybrid structures
    print("Modelling hybrid (memory L chain) antibodies...")
    counter = ("%04i" % i for i in count(1))
    for key, value in hybrids_l.items():
        number = next(counter)
        sequences = {"H": str(value.split("-|")[1]), "L": str(key.split("-|")[1])}
        output_file = (fasta.replace(fasta.split("_")[-1], "") + "hyb-L_" + number + ".pdb").replace("bs" + base_structure[0] + "_", "")

        try:
            with ClearCache():
                antibody = predictor.predict(sequences)
                antibody.save(output_file, check_for_strained_bonds=False)
        except:
            print("Modelling failure in hybrid stage (L). Trying again.")
            with ClearCache():
                antibody = predictor.predict(sequences)
                antibody.save(output_file, check_for_strained_bonds=False)

        make_fasta = output_file.replace(".pdb", ".fasta")
        with redirect_stdout(open(make_fasta, "w", newline="")):
            print(">" + str(value.split("-|")[0]))
            print(str(value.split("-|")[1]))
            print(">" + str(key.split("-|")[0]))
            print(str(key.split("-|")[1]))


def fullmem_rec():
    # Produces the selected amount of "fullmemory" antibody structures, using the memory chains of each of the hybrid
    # antibodies with the lowest RMSD in relation to the original structure
    print("Modelling full memory antibodies...")
    for low in lowest_h:
        fasta_name = low.replace(".pdb", ".fasta")
        heavy_seq_mem = "heavyseq"
        heavy_head_mem = "heavyhead"
        for heavy_fasta in glob.glob(fasta_name):
            open_heavy_fasta = open(os.path.join(os.getcwd(), heavy_fasta), "r", newline='')
            for seq in open_heavy_fasta:
                if seq.startswith(">") and "|H|" in seq:
                    heavy_seq_mem = next(open_heavy_fasta).replace("\n", "")
                    heavy_head_mem = seq.replace("\n", "").replace(">", "")

            for low_l in lowest_l:
                light_fasta_name = low_l.replace(".pdb", ".fasta")
                open_light_fasta = open(os.path.join(os.getcwd(), light_fasta_name), "r", newline='')
                for seq_li in open_light_fasta:
                    if seq_li.startswith(">") and "L" in seq_li:
                        light_seq_mem = next(open_light_fasta).replace("\n", "")
                        light_head_mem = seq_li.replace("\n", "").replace(">", "")

                        h_id = fasta_name.split("_")[3:5]
                        l_id = light_fasta_name.split("_")[3:5]

                        full_mem = "fullmem_" + h_id[0] + h_id[1].replace(".fasta", "") + "_" + l_id[0] + l_id[
                            1].replace(
                            ".fasta", "")
                        full_mem_fasta = (fasta.replace(fasta.split("_")[-1], full_mem)).replace("bs" + base_structure[0] + "_", "") + ".pdb"

                        sequences = {"H": heavy_seq_mem, "L": light_seq_mem}

                        output_file = full_mem_fasta

                        try:
                            with ClearCache():
                                antibody = predictor.predict(sequences)
                                antibody.save(output_file, check_for_strained_bonds=False)
                        except:
                            print("Modelling failure in full memory stage. Trying again.")
                            with ClearCache():
                                antibody = predictor.predict(sequences)
                                antibody.save(output_file, check_for_strained_bonds=False)

                        make_fasta = output_file.replace(".pdb", ".fasta")
                        with redirect_stdout(open(make_fasta, "w", newline="")):
                            print(">" + heavy_head_mem)
                            print(heavy_seq_mem)
                            print(">" + light_head_mem)
                            print(light_seq_mem)


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

def mem_nai_swap():
    # Take the memory sequence selected and swap its CDRs with those of the original antibody
    sequence = next(open_db).replace("\n", "").replace("\r", "")
    CDRH1_res = [26, 32]
    CDRH2_res = [52, 56]
    CDRH3_res = [95, 102]
    H4_res = [71, 75]
    CDRL1_res = [24, 34]
    CDRL2_res = [50, 56]
    CDRL3_res = [89, 97]
    L4_res = [66, 71]

    start = []
    new_CDR1 = []
    fw2 = []
    new_CDR2 = []
    fw3_1 = []
    new_DE = []
    fw3_2 = []
    new_CDR3 = []
    end = []

    if "FULL_H" in header:
        run_anarci_mem = "ANARCI -i " + sequence + " -s martin -o heavy_mem.anarci -p 8"
        subprocess.run(run_anarci_mem, shell=True, capture_output=True)

        open_num_mem = open("heavy_mem.anarci", "r").readlines()
        open_num_nai = open(fasta.replace(".fasta", ".anarci"), "r").readlines()

        for e in open_num_nai:
            if e.startswith("H"):
                num_res = int(e.split(" ")[1])
                amino = (e.split(" ")[-1]).replace("\n", "")

                if num_res >= CDRH1_res[0] and num_res <= CDRH1_res[-1]:
                    new_CDR1.append(amino)
                if num_res >= CDRH2_res[0] and num_res <= CDRH2_res[-1]:
                    new_CDR2.append(amino)
                if num_res >= H4_res[0] and num_res <= H4_res[-1]:
                    new_DE.append(amino)
                if num_res >= CDRH3_res[0] and num_res <= CDRH3_res[-1]:
                    new_CDR3.append(amino)

        for d in open_num_mem:
            if d.startswith("H"):
                num_res = int(d.split(" ")[1])
                amino = (d.split(" ")[-1]).replace("\n", "")
                if num_res < CDRH1_res[0]:
                    start.append(amino)
                if num_res > CDRH1_res[1] and num_res < CDRH2_res[0]:
                    fw2.append(amino)
                if num_res > CDRH2_res[1] and num_res < H4_res[0]:
                    fw3_1.append(amino)
                if num_res > H4_res[1] and num_res < CDRH3_res[0]:
                    fw3_2.append(amino)
                if num_res > CDRH3_res[1]:
                    end.append(amino)

    if "FULL_L" in header:
        run_anarci_mem = "ANARCI -i " + sequence + " -s martin -o light_mem.anarci -p 8"
        subprocess.run(run_anarci_mem, shell=True, capture_output=True)

        open_num_mem = open("light_mem.anarci", "r").readlines()
        open_num_nai = open(fasta.replace(".fasta", ".anarci"), "r").readlines()

        for e in open_num_nai:
            if e.startswith("L"):
                num_res = int(e.split(" ")[1])
                amino = (e.split(" ")[-1]).replace("\n", "")

                if num_res >= CDRL1_res[0] and num_res <= CDRL1_res[-1]:
                    new_CDR1.append(amino)
                if num_res >= CDRL2_res[0] and num_res <= CDRL2_res[-1]:
                    new_CDR2.append(amino)
                if num_res >= L4_res[0] and num_res <= L4_res[-1]:
                    new_DE.append(amino)
                if num_res >= CDRL3_res[0] and num_res <= CDRL3_res[-1]:
                    new_CDR3.append(amino)

        for d in open_num_mem:
            if d.startswith("L"):
                num_res = int(d.split(" ")[1])
                amino = (d.split(" ")[-1]).replace("\n", "")
                if num_res < CDRL1_res[0]:
                    start.append(amino)
                if num_res > CDRL1_res[1] and num_res < CDRL2_res[0]:
                    fw2.append(amino)
                if num_res > CDRL2_res[1] and num_res < L4_res[0]:
                    fw3_1.append(amino)
                if num_res > L4_res[1] and num_res < CDRL3_res[0]:
                    fw3_2.append(amino)
                if num_res > CDRL3_res[1]:
                    end.append(amino)

    newseq = start + new_CDR1 + fw2 + new_CDR2 + fw3_1 + new_DE + fw3_2 + new_CDR3 + end
    new_sequence.append("".join([str(item) for item in [w for w in newseq if not w == "-"]]))


def pymol_aln():
    # Each of the hybrid and fullmemory structures are aligned with the original; the RMSD is stored in a log file
    CDRH1_res = [26, 32]
    CDRH2_res = [52, 56]
    CDRH3_res = [95, 102]
    H4_res = [71, 75]
    CDRL1_res = [24, 34]
    CDRL2_res = [50, 56]
    CDRL3_res = [89, 97]
    L4_res = [66, 71]

    regions = ["CDRH1_res", "CDRH2_res", "CDRH3_res", "CDRL1_res", "CDRL2_res", "CDRL3_res"]

    with redirect_stdout(open("pymol_align_mem.pml", "w", newline="")):
        print("load " + fasta.replace(".fasta", ".pdb"))
        print("alter (chain H), chain=\"Z\"")
        print("alter (chain L), chain=\"Y\"")
        print("load " + pdb)
        print("align " + fasta.replace(".fasta", "") + ", " + pdb.replace(".pdb", ""))
        for region in regions:
            chn = "chn"
            alt_chn = "altchn"

            if "H" in region:
                chn = "H"
                alt_chn = "Z"
            if "L" in region:
                chn = "L"
                alt_chn = "Y"
            print("select ab1_" + region + ", chain " + alt_chn + " and resi " + str(locals()[region][0]) + "-" + str(
                locals()[region][-1]) + " and name CA")
            print(
                "select ab2_" + region + ", chain " + chn + " and resi " + str(locals()[region][0]) + "-" + str(locals()[region][
                                                                                                          -1]) + " and name CA")
        for region in regions:
            print("rms_cur ab1_" + region + ", ab2_" + region + ", matchmaker=4")

    run_pml_cmd = str(pymol_command + " -qc pymol_align_mem.pml")
    out = subprocess.run(run_pml_cmd, shell=True, capture_output=True)
    log_pml = out.stdout
    total_rmsd = float(str(log_pml).split("Executive: RMSD =")[1].split("(")[0].replace(" ", ""))
    h1_rmsd = float(str(log_pml).split("Executive: RMSD =")[2].split("(")[0].replace(" ", ""))
    h2_rmsd = float(str(log_pml).split("Executive: RMSD =")[3].split("(")[0].replace(" ", ""))
    h3_rmsd = float(str(log_pml).split("Executive: RMSD =")[4].split("(")[0].replace(" ", ""))
    l1_rmsd = float(str(log_pml).split("Executive: RMSD =")[5].split("(")[0].replace(" ", ""))
    l2_rmsd = float(str(log_pml).split("Executive: RMSD =")[6].split("(")[0].replace(" ", ""))
    l3_rmsd = float(str(log_pml).split("Executive: RMSD =")[7].split("(")[0].replace(" ", ""))
    rmsd_all = str(total_rmsd) + "|" + str(h1_rmsd) + "|" + str(h2_rmsd) + "|" + str(h3_rmsd) + "|" + str(
        l1_rmsd) + "|" + str(l2_rmsd) + "|" + str(l3_rmsd)
    rmsd_list = [h1_rmsd, h2_rmsd, h3_rmsd, l1_rmsd, l2_rmsd, l3_rmsd]

    log_pml_prt = open("pymol_out_" + fasta.replace(".fasta", ".log").replace("bs" + base_structure[0] + "_", ""), "a", newline="")
    with redirect_stdout(log_pml_prt):
        if not any(rm > cdr_rmsd_limit for rm in rmsd_list) and not total_rmsd > conf_rmsd_limit:
            print(pdb + ":" + rmsd_all)
        else:
            print(pdb + ":" + rmsd_all + "|HIGH_RMSD")


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
    print("Starting evaluation protocol for structure: " + low)
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
        print("Minimization failure. Check " + low + "_error.tar for more information.")
        for rep in glob.glob("out-min.pdb"):
            rep_chain = "CHN_" + rep
            rep_score = "SCO_" + rep
            cmd_del_int = "rm " + rep_chain + " " + rep + " " + rep_score + " out* *.out *.info *.rst7 *.crd *.leap *.pqr *.parm7 *.mdinfo *.mdout *.nc temp*"
            cmd_bkp_error = "tar -cvf " + low + "_error.tar " + rep_chain + " " + rep + " " + rep_score + " out* *.out *.info *.rst7 *.crd *.leap *.pqr *.parm7 *.mdinfo *.mdout *.nc temp*"
            subprocess.run(cmd_bkp_error, shell=True, capture_output=True)
            subprocess.run(cmd_del_int, shell=True, capture_output=True)


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
        curr_log = "swap_" + st + ".log"
        open_curr_log = open(curr_log, "r").readlines()
        new_file_name = latest_ab.replace("bs" + base_structure[0] + "_", "")

        for z in open_curr_log:
            if z.startswith(fasta.replace(".fasta", ".anarci")):
                latest_score = (z.split("|")[1]).replace("\n", "")

        kT = 0
        if scoring_method == "ref15":
            kT = 1.0
        if scoring_method == "csm":
            kT = 0.593

        updt_log = open(curr_log, "a", newline="")

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
            print(new_file_name + "|" + str(float(av_new_score)) + "|Modification approved")
            with redirect_stdout(updt_log):
                print("\t" + new_file_name + "|" + str(float(av_new_score)) + "|PASS")
        if not approval_condition:
            print(new_file_name + "|" + str(float(av_new_score)) + "|Modification rejected")
            with redirect_stdout(updt_log):
                print("\t" + new_file_name + "|" + str(float(av_new_score)) + "|FAIL")

        # Remove scoring files from the current cycle
        for rep in glob.glob("out-min.pdb"):
            rep_chain = "CHN_" + rep
            rep_score = "SCO_" + rep
            cmd_del_int = "rm " + rep_chain + " " + rep + " " + rep_score + " out* *.out *.info *.rst7 *.crd *.leap *.pqr *.parm7 *.mdinfo *.mdout *.nc temp*"
            subprocess.run(cmd_del_int, shell=True, capture_output=True)


##########
# Reading configuration file to determine number of hybrid produced for each chain (total number of hybrid structures is always 2x the specified number);
# and how many of the lowest RMSD memory chains will be combined to produce full memory antibodies. The total number of full memory structures will always be
# the square of the specified number (fullmem = 5 means 5*5=25 fullmem structures will be produced). Note: fullmem must always be equal to or smaller than hyb.
hyb = "cfg:memory_frame|num_hyb_struc"
fullmem = "cfg:memory_frame|num_fullmemory"
fwrk_mode = "cfg:memory_frame|fwrk_mode"
conf_rmsd_limit = "cfg:memory_frame|conf_rmsd_limit"
cdr_rmsd_limit = "cfg:memory_frame|cdr_rmsd_limit"
scoring_method = "cfg:scoring_method"
steps = "cfg:steps"
scoring_strictness = "cfg:scoring_strictness"
approval_threshold = 0
ph = "7"
pymol_command = ""

# Making sure that fullmem is never bigger than hyb
if fullmem > hyb:
    fullmem = hyb

cfg = open("swap_settings.cfg").readlines()

for setting in cfg:
    if "memory_frame|num_hyb_struc" in setting:
        hyb = int(setting.split("=")[1].replace("\n", ""))
    if "memory_frame|num_fullmemory" in setting:
        fullmem = int(setting.split("=")[1].replace("\n", ""))
    if "memory_frame|fwrk_mode" in setting:
        fwrk_mode = str(setting.split("=")[1].replace("\n", ""))
    if "memory_frame|conf_rmsd_limit" in setting:
        conf_rmsd_limit = float(str(setting.split("=")[1].replace("\n", "")))
    if "memory_frame|cdr_rmsd_limit" in setting:
        cdr_rmsd_limit = float(str(setting.split("=")[1].replace("\n", "")))
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

print("Starting new step: Framework maturation")

if scoring_method == "ref15":
    from pyrosetta import *
    init("-mute all")

base_structure = []
# Stores which and how many antibodies are being developed simultaneously (numbered 001, 002...)
structures = []

step_list = steps.split("|")
if "" in step_list:
    step_list.remove("")
bs_index = step_list.index("fwk")
base_structure.append(str(bs_index + 1))
base_structure.append(str(bs_index + 2))

for struc in glob.glob("bs" + base_structure[0] + "_*.anarci"):
    structure = struc.split("_")[1]
    if not structure in structures:
        structures.append(structure)

for st in structures:
    orig = "original structure"
    # Stores the sequence and header of each chain of each hybrid antibody produced
    hybrids_h = {}
    hybrids_l = {}
    # Opens log file with scoring results of each antibody being developed
    current_log = "swap_" + st + ".log"
    open_log = open(current_log, "r").readlines()

    # Stores the name of the starting antibody sequence
    for log in open_log:
        if log.startswith("bs" + base_structure[0]):
            orig = log.split("|")[0].replace(".anarci", ".fasta")

    for fasta in glob.glob(orig):
        open_fasta = open(os.path.join(os.getcwd(), fasta), "r", newline='')

        H_len = []
        L_len = []
        best_h = []
        best_l = []

        nai_H = []
        nai_L = []
        # Stores the lengths of each of the original antibody chain's CDRs, their sequences, and headers
        for x in open_fasta:
            if x.startswith(">"):
                chain = x.split("|")[1]
                lengths = "|".join([str(item) for item in [p for p in x.split("|")[2:5]]])
                chain_lengths = chain + "|" + lengths + "|"
                if chain == "H":
                    nai_line_id = x.replace(">", "").replace("\n", "").replace("\r", "")
                    H_len.append(chain_lengths)
                    nai_H.append(nai_line_id + "-|" + next(open_fasta).replace("\n", ""))
                if chain == "L":
                    nai_line_id = x.replace(">", "").replace("\n", "").replace("\r", "")
                    L_len.append(chain_lengths)
                    nai_L.append(nai_line_id + "-|" + next(open_fasta).replace("\n", ""))

        # Runs blast search to find the most similar (to the original antibody) memory cell antibodies in the chosen database. Results are divided by chain.
        if fwrk_mode == "free":
            blast = str("blastp -query " + fasta + " -db ref_seq/mem_all_mut_100 -evalue 1e-6 -num_threads 6 -out " + fasta.replace(".fasta", ".blast"))
            subprocess.run(blast, shell=True, capture_output=True)
        if fwrk_mode == "directed":
            blast = str("blastp -query " + fasta + " -db ref_seq/mem_all_conf_mut_100 -evalue 1e-6 -num_threads 6 -out " + fasta.replace(".fasta", ".blast"))
            subprocess.run(blast, shell=True, capture_output=True)

        open_blast = open(fasta.replace(".fasta", ".blast"), "r").readlines()
        # Finds and stores the identifiers of the database sequences with highest identity
        for line in open_blast:
            if line.startswith("FULL_H"):
                line_id = line.split("|")[4]
                best_h.append(line_id)

        for line in open_blast:
            if line.startswith("FULL_L"):
                line_id = line.split("|")[4]
                best_l.append(line_id)

        for mem_data in glob.glob("ref_seq/mem_full_all_mut.fasta"):
            open_db = open(os.path.join(os.getcwd(), mem_data), "r", newline='')
            for header in open_db:
                if header.startswith(">"):
                    for id in best_h[:hyb]:
                        if id in header:
                            new_sequence = []
                            # Takes the memory sequence selected and swap its CDRs with those of the original antibody
                            mem_nai_swap()
                            # Stores sequence and header information for both the memory and the naive (original) chains
                            hybrids_h[nai_H[-1].split("-|")[0] + id + "-|" + new_sequence[-1]] = nai_L[-1]
                            new_sequence.clear()
                    for id in best_l[:hyb]:
                        if id in header:
                            new_sequence = []
                            mem_nai_swap()
                            hybrids_l[nai_L[-1].split("-|")[0] + id + "-|" + new_sequence[-1]] = nai_H[-1]
                            new_sequence.clear()
        # Model the memory/naive hybrid structures
        # One chain will have a memory antibody framework + the CDRs of the original antibody; the other will be the same as the original
        model_hyb_h()
        model_hyb_l()
        # Each of the hybrid structures are aligned with the original; the RMSD is stored in a log file
        pdb_chains = [st + "*_hyb-L_*.pdb", st + "*_hyb-H_*.pdb"]
        for chains in pdb_chains:
            for pdb in glob.glob(chains):
                pymol_aln()

        pml_rmsd_h = {}
        pml_rmsd_l = {}
        open_pyml = open("pymol_out_" + fasta.replace(".fasta", ".log").replace("bs" + base_structure[0] + "_", ""), "r").readlines()
        # The name of each hybrid structure is stored, along with their RMSD value
        for lines in open_pyml:
            name = lines.split(":")[0]
            rmsd_value = lines.split(":")[1].replace("\n", "").split("|")[0]
            conf_check = lines.split(":")[1].replace("\n", "").split("|")[-1]
            if "_hyb-H_" in name and not conf_check == "HIGH_RMSD":
                pml_rmsd_h[name] = rmsd_value
            if "_hyb-L_" in name and not conf_check == "HIGH_RMSD":
                pml_rmsd_l[name] = rmsd_value
        # The hybrid structures with the lowest RMSD are selected and their memory chains are combined to form "fullmemory" antibodies
        lowest_h = [y[0] for y in sorted(pml_rmsd_h.items(), key=lambda x: x[1])[:fullmem]]
        lowest_l = [y[0] for y in sorted(pml_rmsd_l.items(), key=lambda x: x[1])[:fullmem]]

        # Model fullmemory antibodies
        fullmem_rec()

        # The fullmemory structures are aligned to the original
        for pdb in glob.glob(st + "*fullmem*.pdb"):
            pymol_aln()

        pml_rmsd_fmem = {}
        open_pyml2 = open("pymol_out_" + fasta.replace(".fasta", ".log").replace("bs" + base_structure[0] + "_", ""), "r").readlines()
        # The selected amount of fullmemory structures with the lowest RMSD are scored, to check their affinity with the target protein
        for fmem in open_pyml2:
            name = fmem.split(":")[0]
            rmsd_value = fmem.split(":")[1].replace("\n", "").split("|")[0]
            conf_check = fmem.split(":")[1].replace("\n", "").split("|")[-1]
            if "fullmem" in fmem and not conf_check == "HIGH_RMSD":
                pml_rmsd_fmem[name] = rmsd_value

        lowest_fmem = [y[0] for y in sorted(pml_rmsd_fmem.items(), key=lambda x: x[1])]

        for low in lowest_fmem:
            latest_ab = low
            latest_ab_complex = "complex_" + latest_ab
            latest_complex = "complex_" + fasta.replace(".fasta", ".pdb")
            temporary = "temp"
            renum = "_renum"

            new_complex(latest_ab, latest_complex, latest_ab_complex, temporary, renum, scoring_strictness, approval_threshold)

        # The fullmemory structures with the best (lowest) scores are stored
        prev_name = "bs" + base_structure[0]
        prev_score = 999999
        any_approval = "apr"
        best_mem = {}
        num_list = []
        run_list = []
        sco_list = []
        open_final_log = open(current_log, "r").readlines()
        for sc in open_final_log:
            if sc.startswith("bs" + base_structure[0] + "_"):
                name = (sc.split("|")[0]).replace("\n", "")
                progress = (name.split("_")[-1]).replace(".anarci", "")
                new_progress = '{:06}'.format(int(progress) + 1)
                run = ""
                if step_list[-1] == "fwk":
                    run = (name.split(progress)[0]).replace("bs" + base_structure[0], "FINAL")
                if not step_list[-1] == "fwk":
                    run = (name.split(progress)[0]).replace("bs" + base_structure[0], "bs" + base_structure[-1])
                run_list.append(run)
                num_list.append(progress)
                num_list.append(new_progress)
                prev_name = name
                prev_score = float((sc.split("|")[-1]).replace("\n", ""))
            if sc.startswith("\t" + st + "_"):
                name = (sc.split("|")[0]).replace("\n", "").replace("\t", "")
                score = float((sc.split("|")[1]).replace("\n", ""))
                approval = sc.split("|")[-1].replace("\n", "")
                if approval == "PASS":
                    any_approval = "PASS"
                    best_mem[name] = score

        if not any(cond.split("|")[-1].replace("\n", "") == "PASS" in cond for cond in [p for p in open_final_log if p.startswith("\t" + st + "_")]):
            best_mem[prev_name] = prev_score
            any_approval = "FAIL"
            num_list.remove(num_list[-1])

        best = min(best_mem, key=best_mem.get)
        best_score = best_mem[best]

        print("Finalizing current structure...")
        # The fullmemory structure with the best score registered in the main log file
        write_final_log = open(current_log, "a", newline="")
        with redirect_stdout(write_final_log):
            print(run_list[-1] + str(num_list[-1]) + ".anarci|" + str(best_score))

        # ANARCI is used to renumber the best fullmemory structure, using the martin scheme
        cmd = str("ANARCI -i " + best.replace(".pdb", ".fasta") + " -s martin -o " + run_list[-1] + str(num_list[-1]) + ".anarci" + " --assign_germline -p 8 --use_species human")
        out = subprocess.run(cmd, shell=True, capture_output=True)
        log9 = (out.stderr).decode("utf-8")

        log3 = open("anarci_out.log", "a", newline="")
        with redirect_stdout(log3):
            print(best)
            print(log9)

        # Renamed copies of the final structure, sequence and numbering file are created, for use in the next steps of the modification process
        rename_cmd = str("cp complex_" + best.replace(".anarci", ".pdb") + " complex_" + run_list[-1] + str(num_list[-1]) + ".pdb")
        rename_cmd2 = str("cp " + best.replace(".anarci", ".fasta").replace(".pdb", ".fasta") + " " + run_list[-1] + str(num_list[-1]) + ".fasta")
        rename_cmd3 = str("cp " + best.replace(".anarci", ".pdb") + " " + run_list[-1] + str(num_list[-1]) + ".pdb")
        subprocess.run(rename_cmd, shell=True, capture_output=True)
        subprocess.run(rename_cmd2, shell=True, capture_output=True)
        subprocess.run(rename_cmd3, shell=True, capture_output=True)

    hybrids_h.clear()
    hybrids_l.clear()
