# ARG_Benchmarking

This repository contains code associated with the paper "Evaluating ARG-based methods for Estimating Population-mean Polygenic Score Histories".

To simulate haplotypes and ARGs according to the haplotypes. Please run `generate_true_samples.R` under `small_example/example`:
```
Rscript generate_true_samples.R --n_chromss 2000 --n_loci 100 --out mssel_out --temp mssel_temp --iter_start 1 --iter_end 100
```
The details of the arguments are provided here:
| flag | details |
|-------| -------|
| --n_chromss | number of chromosomes |
| --n_loci | number of loci |
| --out | directory to store output files|
| --temp | directory to store temporary files |
| --iter_start | start number of trait |
|--iter_end | end number of trait|

To sample and input simulated haplotypes to ARG-estimation methods, run `loop_pheno_sim_reps.R` under `small_example/example`:
```
Rscript loop_pheno_sims_reps.R --n_chromss 200 --n_loci 100 --out singer_output --temp temp_out --sample_seq_num 200 --singer yes --iter_start 1 --iter_end 100
```
The above command samples 200 out of 2000 simulated haplotypes and runs SINGER with these haplotypes for 100 iterations. It saves SINGER estimated trees in a list.

To apply Edge & Coop (2019) estimators on the estimated marginal tree, run `run_analysis2.R` under `small_example/example`:
```
Rscript run_analysis2.R --n_chromss 200 --n_loci 100 --out singer_output --iter_start 1 --iter_end 100 --singer yes --save singer_analyzed
```

