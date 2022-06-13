function [C, tmp, popSizeData] = loadData(cont_matrix_fname, pop_dist1_fname, pop_dist2_fname)
% Function to load data required in the SimLeaky simulations from files.
% INPUTS:
% - cont_matrix_fname: contact matrix
% - pop_dist1_fname: benchmark population distribution for contact matrix
% - pop_dist2_fname: population distribution
% OUTPUTS:


C = readmatrix(cont_matrix_fname); % Get Prem et al contact matrix from data folder
tmp = readmatrix(pop_dist1_fname); % Load NZ population structure from data folder
popSizeData = readmatrix(pop_dist2_fname); % Load NZ population structure from data folder

end