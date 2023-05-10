Use generate_SNP_pileups.py to generate training features for SNP calling, and then run model_run.py to training a model.
Use generate_indel_pileups.py or generate_indel_pileups_hifi.py (for PacBio HiFi)  to generate training features for indel calling, and then run model_run_indels.py to training a model.

You can print help for each script using `--help`, e.g. `python generate_SNP_pileups.py --help`
This code uses tensorflow 1.13 for training a model.
