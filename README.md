# Introduction
This is an anonymized and reduced repository for the paper 
Quantifying Instance Hardness of Protein Folding within the HP model. It offers
the bare minimum code required for reproducing the required results and 
generating the figures used in the paper.

The code uses a protein folding framework called prospr. We wrote this framework
and made it open source available on GitHub. Due to the required anonymity, we
decided not to link it in the paper. This is an explicit note to NOT look up the 
framework for maintaining the anonymity.


# Directory structure.
All experiments were executed on the Lisa computing cluster of SURF, which is
also linked in the paper. The algorithms can be found in the `/code/experiments`
folder. Note that there is a boolean argument for each algorithm to specify 
whether the code is executed on the cluster or locally. There are two versions
of the depth-first branch-and-bounds algorithm, one serial and one parallel. The
parallel version is used to fold 16 proteins at the same time.

The experiments are executed through jobs, which are located in the `/jobs`
folder. Each job has a unique name, with in the folder a shell script with the
same name. There is a job for each length of protein, and one folder with all
results combined (`/job/s8_1000_total`). Note that the shell script for this
folder is empty and just functions as a placeholder.

A job is executed via the `run_program.py` program, which also allows for local
execution. This script is located in the main directory. Next to it, you will
find the `run_statistics.py` file, which can be used to compute statistics from
the collected data. There is also `run_visualizer.py`, which can be used to
create all figures from the paper. Both `run_statistics.py` and 
`run_visualizer.py` use helper code, which is located in the `/code/helpers` 
directory.


# Collecting data
A job can be locally executed through:
```shell
python run_program.py dfs_bnb_parallel -j s8_1000_<length> -le <length>
```
The code will overwrite the data currently in the jobs folder. We combined all
results in the `/jobs/s8_1000_total` folder for computing the statistics and
generating the figures.

TODO: SLURM!
Running on the Lisa cluster can be done by simply executing a job's shell script.


# Reproducing statistics
The `run_statistics.py` script compute the paper's statistics values for the: 
 - goodness-of-fit tests for distributions
 - goodness-of-fit tests for the parameters
 - structural statistics (P-ratio, bonds 1st half, and bonds with last quarter)
 - the statistics on recursions

The goodness-of-fit tests for the distributions can be executed through:
```shell
python run_statistics.py comp_dist_gof -j s8_1000_total -le 10 15 20 25
```

The goodness-of-fit tests for the parameters can be executed through:
```shell
python run_statistics.py comp_params_gof -j s8_1000_total -le 10 15 20 25
```

The structural statistics and number of recursions tests can be executed through:
```shell
python run_statistics.py comp_stats_extrema -j s8_1000_total -le 10 15 20 25 
```

# Reproducing figures
The `run_visualizer.py` script generates the papers used in the paper. Figures
1a and 1b are created by hand. The rest of figures are generated with the 
script, which include: 
 - the most right front figure (Fig. 1c)
 - the instance hardness distributions (Fig. 2a)
 - the fitted parameters (Fig 2b)
 - the conformation plot (Fig. 3)

The right front figure can be generated through:
```shell
python run_visualizer.py front_figure
```

The instance hardness distributions can be generated through:
```shell
python run_visualizer.py fit_group -j s8_1000_total -s --no_show
```

The fitted parameters can be generated through:
```shell
python run_visualizer.py params -j s8_1000_total -s --no_show
```

The fitted parameters can be generated through:
```shell
python run_visualizer.py plot_hardest -j s8_1000_total -s --no_show
```
