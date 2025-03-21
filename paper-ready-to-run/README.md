# Ready-to-run examples of the pipeline

## Ab-SELDON versions
To quickly test the pipeline, two versions are available in separate folders. 
To run an optimization example, simply enter either of these folders, extract the `ab-seldon-databases.rar` file and execute the main script with the command `sh seldon.sh`

### ab-seldon-paper-v6.2
Version used in the paper's optimization experiments. 
- While the individual modules can be run independently and in any order, this requires manual setup. 
- The inputs must be named `complex_[NAME].pdb` and `[NAME].fasta`.
- Uses the same settings used in the optimization tests of the paper. Folders with settings files for diversity-guided (`diversity/`) or random (`random/`) optimizations are available.

#
### ab-seldon-ready-v6.45
Most current version of the pipeline as of march 2025. This new version has a few differences compared to the paper version.
- The use of the optimization modules in non-default order and number has become much simpler with the introduction of the `steps` parameter. For more information, visit the [instructions for the configuration file](https://github.com/SFBBGroup/Ab-SELDON/blob/main/configuration_file_instructions.md).
- Customizable modification approval thresholds have been added to the pipeline, through the `scoring_strictness` and `approval_threshold` parameters in the configuration file.
- The input files must have the same name ([NAME].pdb and [NAME].fasta)
- A bug has been corrected in the OAS CDR grafting module (`swap_bnk_v6.45.py`), where CDRs whose germline was represented in the database by too few sequences could cause an infinite loop.
- Nevertheless, these differences are not relevant to the tests made with the previous version.

Alternative .pdb complexes and associated .fasta files are also available in the `example-complexes/` subfolders. To optimize any of these complexes, simply copy the fasta and pdb files to the main pipeline folder and edit the `prepare|input_name` parameter accordingly, before running the main script.
To test the basic functionality of the pipeline's modules, you can set the `n_cycles`/`max_cycles` and `num_fullmemory` parameters to reduce the number of optimization cycles. Suggested values are:
- `n_cycles`/`max_cycles` = 2
- `num_fullmemory` = 1

For more information about these parameters and how to edit the configuration file, visit the [instructions for the configuration file].(https://github.com/SFBBGroup/Ab-SELDON/blob/main/configuration_file_instructions.md).
