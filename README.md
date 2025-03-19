# Ab-SELDON
Antibody Structural Enhancement Leveraging Diversity for Optimization of iNteractions

Paper: Ab-SELDON: Leveraging Diversity Data for an Efficient Automated Computational Pipeline for Antibody Design

## Abstract
The utilization of predictive tools has become increasingly prevalent in the development of biopharmaceuticals, reducing the time and cost of research. However, most methods for computational antibody design are hampered by their reliance on scarcely available antibody structures, potential for immunogenic modifications, and a restricted exploration of the paratope's potential chemical and conformational space. We propose Ab-SELDON, a modular and easily customizable antibody design pipeline capable of iteratively optimizing an antibody-antigen (Ab-Ag) interaction in five different modification steps, including CDR and framework grafting, and mutagenesis. The optimization process is guided by diversity data collected from millions of publicly available human antibody sequences. This approach enhanced the exploration of the chemical and conformational space of the paratope during computational tests involving the optimization of an anti-HER2 antibody. Optimization of another antibody against Gal-3BP stabilized the Ab-Ag interaction in molecular dynamics simulations. Tests with SKEMPI’s Ab-Ag mutations also demonstrated the pipeline’s ability to correctly identify the effect of most mutations.

## Installation

### Requirements

This pipeline requires the following software:
- Python (tested versions: 3.11.6; 3.12.3)
-	[ImmuneBuilder](https://github.com/oxpig/ImmuneBuilder) (tested versions: 1.0.1; 1.1.1)
-	[AmberMD and AmberTools](https://ambermd.org/GetAmber.php) (tested versions: Amber 22; Amber 24)
-	[PyRosetta](https://www.pyrosetta.org/downloads#h.6vttn15ac69d) (tested versions: pyrosetta-2023.36; pyrosetta-2024.19)
-	[PDB2PQR](https://pdb2pqr.readthedocs.io/en/latest/getting.html#python-package-installer-pip) (tested version: 3.6.1)
-	[ANARCI](https://github.com/oxpig/ANARCI) 
-	[BLAST](https://www.ncbi.nlm.nih.gov/books/NBK52640/) (tested versions: 2.12.0; 2.15.0)
-	[pdb-tools](https://github.com/haddocking/pdb-tools) (tested version: 2.5.0)
-	[PyMOL](https://github.com/schrodinger/pymol-open-source) (tested versions: 2.6.0a0; 3.0)

### Setup

Ensure all requirements are installed correctly. In particular, use `python -m openmm.testInstallation` to see if OpenMM (on which ImmuneBuilder depends) is successfully computing forces with CUDA.

Once the required softwares are installed, enter the ab-seldon/ folder and extract the pipeline's databases with:

` $ unrar x ab-seldon-databases.rar `

## Usage
The pipeline takes as input:
1) A fasta file [NAME].fasta with ONLY the sequence of the initial antibody that will be optimized. The heavy chain MUST come before the light chain;
2) A PDB file [NAME].pdb (same name as the fasta) with the initial antibody-antigen complex whose interaction will be optimized. It must contain only one antibody molecule and its antigen. The antibody chains must be named H and L. It must not contain heteroatoms, only proteins.

These files must be put into the pipeline's main folder, ab-seldon/

To configure your optimization run, you must edit the configuration file (swap_settings.cfg). If you wish to run the pipeline with its default settings, simply edit the first parameter (prepare|input_name=) to replace the [NAME] with the name of your fasta/pdb file (eg. prepare|input_name=[NAME] becomes prepare|input_name=6phb if your inputs are named 6phb.fasta and 6phb.pdb)

By default, it is assumed that PyMOL can be executed on the terminal with the `pymol` command. If not, change the command in pymol_command= (eg. to `pymol.exe`)

For more information about how to customize other aspects of your optimization run (such as the optimization steps you wish to run, the number of cycles executed in each step, weightings of the probabilistic modification site selection, etc.), please read the deatils of each of the settings parameters [here](LINK).

After these steps, run the pipeline by simply executing the main script:

` $ sh seldon.sh `

After the optimization process is concluded, the output files with the optimized antibody and antibody-antigen complex will be identified wth a `FINAL` prefix:
- complex_FINAL(...).pdb
- FINAL(...).fasta

Additionally, a file called `swap_001.log` documents which modifications were tested, approved or rejected, along with their associated Ab-Ag predicted interaction energy. 


