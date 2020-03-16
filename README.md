# Trip-Centrality
This repository contains the code used for the data analysis in the paper 
Zaoli, S., Mazzarisi, P.; Lillo, F. Trip Centrality: walking on a temporal multiplex with non-instantaneous link travel time. Sci Rep 9, 10570 (2019). https://doi.org/10.1038/s41598-019-47115-6

The code is written in Matlab. 
The following files are found in the repository:
day.mat: Sample of the data used in the article
adjacency.m : script to compute the adjacency matrices A_sched and A_real
compute_TC.m : script to compute Trip centrality (scheduled and realised) from the adjacency matrices
airports.mat : list of airports
airlines.mat : list of airlines

The input file 'day.mat' contains a table of flights (1 day). For each flight, the table contains scheduled and realised departure and arrival times, departure and arrival airport and airline. Airports and airlines are indicated by numbers which corresponds to their position in the lists 'airports.mat' and 'airlines.mat'. 

If you use this code, please cite the airticle as: 
Zaoli, S., Mazzarisi, P.; Lillo, F. Trip Centrality: walking on a temporal multiplex with non-instantaneous link travel time. Sci Rep 9, 10570 (2019). https://doi.org/10.1038/s41598-019-47115-6
