function transition_paths_driver()
%file to load in user-defined parameters for calculating transition paths in
%'transition_paths.m'

%regions: array where each column corresponds to a data point(D-dim
%vector) in the ambient space R^D. Each data point should define an attractor of
%the dynamical system of interest.

%dist: distance which defines the region of state space assigned to the columns
%of 'regions' array, i.e the open ball with radius dist around a given
%attractor should define it's basin 

%Xzero: desired intial starting point in ambient space for transition paths
%simulation, column vector in R^D

%Example: SDE dXt = -grad U dt + dWt, U has minima at
% p1 = [-1 2 1]^T, and  p2 = [5 9 2]^T, with each having a basin of attraction
% greater than radius 100 around it

%p1 = [-1;2; 1];
%p2 = [5;9;2];
%regions = [p1 p2];
%dist = 100;
%Xzero = p1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Re-write the values of the below lines to input parameters
regions = [0 1];
dist = 0.25;
Xzero = 0.2;
T=10000;
%%save parameters into .mat file for use with 'transition_paths.m'
save('transition_params.mat');


end

