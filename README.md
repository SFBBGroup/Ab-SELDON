# Ab-SELDON
Antibody Structural Enhancement Leveraging Diversity for Optimization of iNteractions

Paper: Ab-SELDON: Leveraging Diversity Data for an Efficient Automated Computational Pipeline for Antibody Design

## Abstract
The utilization of predictive tools has become increasingly prevalent in the development of biopharmaceuticals, reducing the time and cost of research. However, most methods for computational antibody design are hampered by their reliance on scarcely available antibody structures, potential for immunogenic modifications, and a restricted exploration of the paratope's potential chemical and conformational space. We propose Ab-SELDON, a modular and easily customizable antibody design pipeline capable of iteratively optimizing an antibody-antigen (Ab-Ag) interaction in five different modification steps, including CDR and framework grafting, and mutagenesis. The optimization process is guided by diversity data collected from millions of publicly available human antibody sequences. This approach enhanced the exploration of the chemical and conformational space of the paratope during computational tests involving the optimization of an anti-HER2 antibody. Optimization of another antibody against Gal-3BP stabilized the Ab-Ag interaction in molecular dynamics simulations. Tests with SKEMPI’s Ab-Ag mutations also demonstrated the pipeline’s ability to correctly identify the effect of most mutations.

## Installation

### Requirements

This pipeline requires the following software:
-	[ImmuneBuilder](https://github.com/oxpig/ImmuneBuilder)
-	[AmberMD and AmberTools](https://ambermd.org/GetAmber.php)
-	[PyRosetta](https://www.pyrosetta.org/downloads#h.6vttn15ac69d)
-	[PDB2PQR](https://pdb2pqr.readthedocs.io/en/latest/getting.html#python-package-installer-pip)
-	[ANARCI](https://github.com/oxpig/ANARCI)
-	[BLAST](https://www.ncbi.nlm.nih.gov/books/NBK52640/)
-	[pdb-tools](https://github.com/haddocking/pdb-tools)
-	[PyMOL](https://github.com/schrodinger/pymol-open-source)

### Setup

Once the required softwares are installed, enter the ab-seldon-v6.4/ folder and extract the pipeline's databases with:

` $ unrar x ab-seldon-databases.rar `

## Usage
The pipeline takes as input:
1) A fasta file [NAME].fasta with ONLY the sequence of the initial antibody that will be optimized. The heavy chain MUST come before the light chain;
2) A PDB file [NAME].pdb (same name as the fasta) with the initial antibody-antigen complex whose interaction will be optimized. It must contain only one antibody molecule and its antigen. The antibody chains must be named H and L. It must not contain heteroatoms, only proteins.

These files must be put into the pipeline's main folder, ab-seldon-v6.4/

To configure your optimization run, you must edit the configuration file (swap_settings.cfg). If you wish to run the pipeline with its default settings, simply edit the first parameter (prepare|input_name=) to replace the [NAME] with the name of your fasta/pdb file (eg. prepare|input_name=[NAME] becomes prepare|input_name=6PHB-clean if your inputs are 6PHB-clean.fasta and 6PHB-clean.pdb)

By default, it is assumed that PyMOL can be executed on the terminal with the `pymol` command. If not, change the command in pymol_command= (eg. to `pymol.exe`)

After these steps, run the pipeline by simply executing the main script:

` $ sh seldon.sh `

After the optimization process is concluded, the output files with the optimized antibody and antibody-antigen complex will be identified wth a `FINAL` prefix:
- complex_FINAL(...).pdb
- FINAL(...).fasta

If you wish to customize other aspects of your optimization run (such as the optimization steps you wish to run, the number of cycles executed in each step, weightings of the probabilistic modification site selection, etc.), please read the documentation of the pipeline [TBA].
