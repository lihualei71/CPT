# Paper Repository for Cyclic Permutation Test (CPT)

This repository contains the code to implement all examples in our paper: [An Assumption-Free Exact Test For Fixed-Design Linear Models With Exchangeable Errors](https://arxiv.org/abs/). 
The folder `code/` contains all R files for implementation and bash files for submitting the jobs to a cluster with [Slurm system](https://slurm.schedmd.com/overview.html). The total computation load is $\sim$20000 CPU hours to reproduce all experimental results in this paper. I used 256 cores for $\sim$7 days. 

- `CPT.R` contains the implementation of cyclic permutation test. It depends on the "gaoptim" package for the genetic algorithm which has been removed from CRAN on 2018-06-17. But it has better performance than other CRAN packages for our purpose so we keep it. The following code can be used to download the package.
```
devtools::install_url("https://cran.r-project.org/src/contrib/Archive/gaoptim/gaoptim_1.1.tar.gz")
```
- `CPT_expr_prepareX.R` contains the template to search for good ordering for CPT using Genetic Algorithm. It takes four external inputs: `Xdist` for the type of design matrices, `ratio` for the ratio between n and p, `ninds` for the number of indices to be tested, `seed` for random seed. The output includes the desing matrix (`X`), a weak ordering by running the genetic algorithm with 1000 samples (`ordering1`) and a strong ordering by running the genetic algorithm with 10000 samples (`ordering`). To reproduce the matrices used in the paper, run the following code on Shell.
```p
cd CPT
mkdir results
cd code
./gen_jobs_prepareX.sh params_prepareX.txt
```
Each row in the file `params_prepareX.txt` contains the four external inputs, namely `Xdist`, `ratio`, `ninds` and `seed`, corresponding to an experiment. The bash file `gen_jobs_prepareX.sh` will automatically generate bash files to call `CPT_expr_prepareX.R` in `results/` folder for each row of `params_prepareX.txt`, submit each to a core, and save the output in the `data/` folder. If the user wants to reduce the number of experiments or add other experiments, the simplest way is to create another .txt file mimicing `params_prepareX.txt` and including the desired experimental setups and run `./gen_jobs_prepareX.sh NEW_TXT_FILE.txt`.
- `CPT_expr.R` contains the code template for all numerical experiments in Section 3 and Appendix A. It takes five external inputs: `Xdist` for the type of design matrices, `epsdist` for the error distribution, `ratio` for the ratio between n and p, `ninds` for the number of indices to be tested, `seed` for random seed. For given inputs, it required `CPT_expr_prepareX.R` to be ran at least once with `Xdist`, `ratio`, `ninds` and `seed` so that its output is stored in `data/` folder; otherwise the job will be halted. The output of `CPT_expr.R` is a data frame containing the inputs as well as the relative signal strength and the power for each method. To reproduce the results in Section 3 and Appendix A, first run the following code on Shell.
```
./gen_jobs_CPT.sh params_CPT_expr.txt
```
Each row in the file `params_CPT_expr.txt` contains the five external inputs, namely `Xdist`, `epsdist`, `ratio`, `ninds` and `seed`, corresponding to an experiment. The bash file `gen_jobs_CPT.sh` will automatically generate bash files to call `CPT_expr.R` in `results/` folder for each row of `params_CPT_expr.txt`, submit each to a core, and save the output in the `data/` folder. If the user wants to reduce the number of experiments or add other experiments, the simplest way is to create another .txt file mimicing `params_CPT_expr.txt` and including the desired experimental setups and run `./gen_jobs_CPT.sh NEW_TXT_FILE.txt`. After the jobs are done, run the following code to aggregate the results from each core.
```
R CMD BATCH agglomerate.R
```
This will generate a file `power_CPT_expr.RData` in `data/` folder. It is currently included for users who do not want to reproduce the experiments.
- `CPT_GA_SS.R` contains the code template to compare genetic algorithm with stochastic search as in Figure 2. It takes four external inputs: `Xdist` for the type of design matrices, `algo` for the algorithm (GA or SS), `popSize` for the population size of GA (not matter for SS), and `seed` for random seed. The output of `CPT_GA_SS.R` is a vector of O* values. To reproduce the results in Figure 2, first run the following code on Shell.
```
./gen_jobs_GA_SS.sh params_GA_SS.txt
```
Each row in the file `params_GA_SS.txt` contains the four external inputs, namely `Xdist`, `algo`, `popSize` and `seed`, corresponding to an experiment. The bash file `gen_jobs_GA_SS.sh` will automatically generate bash files to call `CPT_GA_SS.R`  in `results/` folder for each row of `params_GA_SS.txt`, submit each to a core, and save the output in the `data/` folder. If the user wants to reduce the number of experiments or add other experiments, the simplest way is to create another .txt file mimicing `params_GA_SS.txt` and including the desired experimental setups and run `./gen_jobs_GA_SS.sh NEW_TXT_FILE.txt`. After the jobs are done, run the following code to aggregate the results from each core.
```
R CMD BATCH agglomerate.R
```
This will generate a file `CPT_GA_SS.RData` in `data/` folder. It is currently included for users who do not want to reproduce the experiments.
- `illustrate_perm.R` produces Figure 1.
- `CPT_GA_plot.R` produces Figure 2.
- `CPT_plot.R` produces all Figure 3-16.
- `expr_functions.R` contains all helper functions.

