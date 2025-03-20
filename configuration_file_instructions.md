## Settings file parameters
The `swap_settings.cfg` file contains the parameters of the optimization run, allowing the user to tailor the pipeline to their use case. Here, these parameters are set to their default values.

### prepare|input_name=[NAME]
- Specifies the input files. Replace the `[NAME]` with the name of your fasta/pdb file (eg. `prepare|input_name=[NAME]` becomes `prepare|input_name=6phb` if your inputs are named `6phb.fasta` and `6phb.pdb`.
### prepare|code_name(optional)=ab
- An additional identifier for the optimization run. Useful to differentiate replicates. If not needed, you can leave it on default (`ab`), but not empty. 
### prepare|number_of_structures=1
- How many optimization runs will be done sequentially for the same input files. If you wish to run more than one optimization in parallel instead of sequentially, leave this parameter on default (`1`), copy the main folder and execute the other optimization from the new location.
### prepare|seldon_version=6.45
- Version of the pipeline being executed.
### steps=h3|rep|bnk|mut|fwk|
- Specifies the optimization modules that will be used on the run and their order. By default (`h3|rep|bnk|mut|fwk|`), it executes CDR H3 grafting (`h3`), representative CDR grafting (`rep`), OAS CDR grafting (`bnk`), mutagenesis (`mut`), and framework grafting (`fwk`), in that order. 
- Changing this parameter to `mut|` will only run the mutagenesis step, for example, while `fwk|h3|mut|` will run the framework grafting, CDR H3 grafting and mutagenesis steps, in that order. Any combination, number and order of optimization modules is allowed by the pipeline.
### scoring_strictness(metro/normal/strict)=metro
- Specifies the strictness of the scoring procedure. 
- `metro` will approve any modifications that reduce the predicted interaction energy, while using the metropolis criterion to decide whether to approve modifications that increase the predicted interaction energy.
- `normal` will approve any modifications that reduce the predicted interaction energy and reject any that don't.
- `strict` will only approve modifications that reduce the predicted interaction energy by a specified amount (the `approval_threshold`) and reject any that don't.
### approval_threshold=-2
- The minimum predicted interaction energy reduction required to approve a modification. This parameter is only used if `scoring_strictness` is set to `strict`.
### minimization_pH=7
- The pH at which the Amber minimization will be done during each cycle.
### pymol_command=pymol
- The command used to run PyMOL on shell. If typing `pymol` on your terminal doesn't activate PyMOL, you must change this parameter to the correct command, or the modification cycle will fail right after the modelling procedure, during Ab-Ag complex assembly.
### scoring_method(ref15/csm)=ref15
- The scoring method used to predict the Ab-Ag interaction energy. Please note that REF15 (`ref15`) runs locally as part of pyrosetta, while CSM-AB (`csm`) uses an [external server](https://biosig.lab.uq.edu.au/csm_ab/api) to calculate the score, taking a longer time to get the result. 
### germ_mode(free/restricted)=restricted
- Defines whether grafted sequences can come from any antibody (free) or only from antibodies with the same germline as the input (restricted).

### swap_h3|n_cycles=150
- Number of optimization cycles dedicated to CDR H3 grafting, if this module is executed.
### swap_h3|length_prob_mode(probabilistic/random/restricted)=restricted
- Determines the possible lengths of H3 CDRs that will be grafted for testing in this step. 
- `restricted` allows only H3 sequences with the same length as the input's CDR H3, or those with one more or one less residue than the input's CDR H3.
- `probabilistic` allows CDR H3 sequences of any length available in the dataset. It probabilistically selects a length based on the distribution of CDR H3 lengths in the dataset (i.e. very short and very long H3 sequences are less likely to be selected, while H3 sequences with 12-15 residues are the most likely to be selected. The weights can be found at ab-seldon/h3_bank/h3_length_distribution.csv.
- `random` allows CDR H3 sequences of any length available in the dataset. The length is chosen at random.
### swap_h3|h3_rmsd_limit=3.0
- Discards any modifications whose resulting antibody model has any residue with an RMSD above the specified value, as measured by ImmuneBuilder's own built-in error prediction. The lower the value, the better the minimum quality of the models will be, but this will also increase the number of modifications rejected without evaluation by the scoring procedure. 

### swap_rep|max_cycles=150
- Number of optimization cycles dedicated to representative CDR grafting, if this module is executed.
### swap_rep|cdr_prob(H1|H2|L1|L2|L3)=393|954|8957|359|241825
- The weights for probabilities of any of the non-H3 CDRs being selected for modification in each cycle. The default values are based on CDR diversity data calculated from human OAS sequences of naïve antibodies of healthy donors. Note that these probabilities are only used after the chain is selected, so H1 will only compete for selection with H2, while L1, L2 and L3 will only compete between themselves.
- To make the CDR selection random, simply change this parameter to `1|1|1|1|1`
- Any probabilities other than 0 are allowed. If, for example, modifications on CDR L1 must be prioritized, this parameter could be set to `1|1|999|1|1`. 
### swap_rep|chain_prob(H|L)=22|16
- The weights for probabilities of either antibody chain being selected for modification.


### swap_bnk|check_conf(yes/no)=yes
- Whether or not the conformation of the new grafted CDR should be checked, to allow only CDRs with the same conformation that existed prior to the beginning of this step (by default, this is the conformation resulting from the representative CDR swapping step). 
### swap_bnk|conf_rmsd_limit=1.5
- Maximum RMSD deviation allowed for the modified CDR, if the `conf_check` parameter is set to `yes`. 
### swap_bnk|n_cycles=150
- Number of optimization cycles dedicated to OAS CDR grafting, if this module is executed.
### swap_bnk|cdr_prob(H1|H2|L1|L2|L3)=393|954|8957|359|241825
- The weights for probabilities of any of the non-H3 CDRs being selected for modification in each cycle. Works the same way as the `swap_rep|cdr_prob` parameter.
### swap_bnk|chain_prob(H|L)=22|16
- The weights for probabilities of either antibody chain being selected for modification. Works the same way as the `swap_rep|chain_prob` parameter.
### swap_bnk|indel_prob(no_indel|del|ins)=955|26|19
- The weights for probabilities of the new CDR sequence having the same length as before (`no_indel`), having one less residue than before (`del`), or having one more residue than before (`ins`). This parameters simulates the occurrence of indels in CDR sequences during the antibody maturation process.

### mut|check_conf(yes/no)=yes
- Whether or not the conformation of the new grafted CDR should be checked. Works the same way as the `swap_bnk|check_conf` parameter.
### mut|conf_rmsd_limit=1.5
- Maximum RMSD deviation allowed for the modified CDR, if the `conf_check` parameter is set to `yes`.
### mut|n_cycles=150
- Number of optimization cycles dedicated to mutagenesis, if this module is executed.
### mut|cdr_prob(H1|H2|H3|L1|L2|L3)=377|2062|141994|5305|183|26120
- The weights for probabilities of any of the CDRs being selected for modification in each cycle. Note that this includes all six CDRs, instead of excluding CDR H3. The default weights reflect the diversities of these CDRs in memory antibodies, instead of naive antibodies like in the other modules.
### mut|chain_prob(H|L)=22|16
- The weights for probabilities of either antibody chain being selected for modification. Works the same way as the `swap_rep|chain_prob` parameter.
### mut|h3_cap(yes/no)=yes
- Whether or not the number of mutations directed to CDR H3 is set to 50% of all tested mutations. This option is available because of the extremely high diversity of CDR H3 compared to H1 and H2, which can lead to situations where these other two CDRs are never selected for mutagenesis. If `no` is selected, the CDR selection probabilities will simply follow the weights.
### mut|res_prob_mode(probabilistic/random)=probabilistic
- Whether or not the mutations on the CDRs will be made at random (random position and random new residue) or will be weighted be the diversities of the different positions and the residues found in them in memory antibodies.


### memory_frame|num_hyb_struc=50
- Number of hybrid structures (half memory, half naïve) that will be generated for each chain. That is, the number of memory heavy and light chain frameworks that will be grafted and modelled. For example, setting this parameter to 50 will result in the modelling of 50 hybrid antibodies with a memory framework in the heavy chain, and 50 antibodies with a memory framework in the light chain. Note that these hybrid antibodies are only modelled and have their conformations evaluated, no scoring procedure is performed. Therefore, this parameter can be set to higher numbers at little computational cost, if desired.
### memory_frame|num_fullmemory=5
- Number of memory frameworks of each chain that will be combined to form variable domains where both framework regions come from memory antibodies. For example, setting this parameter to 5 will result in the modelling and scoring of (5 heavy * 5 light) = 25 antibodies whose framework comes from memory sequences. 
### memory_frame|fwrk_mode(free/directed)=directed
- Whether the dataset of all memory antibody sequences (`free`), or the dataset of antibodies with a minimum number of mutations and mutations in important positions (`directed`) will be used as a source of framework regions for the grafting process.
### memory_frame|conf_rmsd_limit=1.5
- Maximum RMSD deviation allowed for the modelled antibodies (whole variable region). Above this threshold, the modification is discarded for excessively changing the conformation of the antibody.
### memory_frame|cdr_rmsd_limit=1.5
- Maximum RMSD deviation allowed for the individual CDRs. If any of the six CDRs reaches above this threshold, the modification is discarded for excessively changing the conformation of the paratope.
