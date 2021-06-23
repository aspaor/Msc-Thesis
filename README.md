# Msc-Thesis
Here you can find all the scripts were used to my thesis and reproduce similar studies.

# Pipeline 1
For pipeline 1 the command is 
``` python3 wrapper_try.py 'param1' 'param2' 'num_of_pops' 'merge_ranges' 'merge_changes' 'migration_changes' 'bottleneck_changes' ```

Wrapper_try.py script contains CLASS.py and score_plotting.py. CLASS.py is the main code of pipeline 1 with the production and the simulations of the data. The score_plotting.py produces all the diagrams based on the output score files. 
In the command above:
- **param1** and **param2**: are the parameters that are going to change during simulations and they can be: length, samples, recomb, mut, N,R, recent_merge, oldest_merge, bottleneck and migration. 
- **num_of_pops**: is the number of populations that are going to be simulated.
- **merge_ranges**: are the time values when populations merged backwards in time.
- **merge_changes**: are the altered values time that populations merged.
- **migration_changes**: are the new migration ratio between pops.
- **bottleneck changes**: are the changes in size of a population during a bottleneck.

Migration_runs.sh, bottleneck_runs.sh, recombination_runs.sh and recent_merge_bottleneck_runs.sh are examples of commands to simulate scenarios with 3 pops. Recent_merge is the time when p2~p3 coalesced.


# Pipeline 2
For pipeline 2 abc.R, plot.R and run_all.R are needed. For reproducing our study ms.sims.stat and migmat_prior are prerequisite. Data were generated via CoMuStat. The abc.R is the main script for the analysis of pipeline 2, plot.R is used for the illustration of the results and run_all.R is the input of the changing parameters.
