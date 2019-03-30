# ATLAS readme:
## Work in Progress!

This repository is an implementation of the ATLAS algorithm from Crosskey and
Maggioni's paper [ATLAS: A Geometric Approach to Learning High-Dimensional
Stochastic Systems Near Manifolds](https://epubs.siam.org/doi/abs/10.1137/140970951). For a more in depth description of the code, please visit the 

## Prerequisites

The codes are all written in MATLAB, and running the code just requires downloading
the src directory.

## Overview

The main files to run the code are:

atlas_driver.m - choose from parameters, set up examples from Crosskey,
Maggioni paper

construction.m - construct the ATLAS for the given example,

learned_simulator.m - Euler-Maruyama simulator for ATLAS learned simulator

The current files for testing output are:

run_atlas_tpaths2d.m - (For 2D potential example below) Computes average
transition times between the potential wells of the original simulator and the
ATLAS simulator

makestats_tpaths2d.m - (For 2D potential example below) Plots the bar graph of
transition averages from run_atlas_2paths2d.m

## Running ATLAS

### 2D, 3-Well Potential example

The line
```
atlas_driver(3);
```
by default:

-sets ATLAS parameters to the Crosskey, Maggioni parameters for ex. 5.3.1 in
the paper 

-constructs the delta-net

-runs construction.m, which constructs the ATLAS and local SDE simulators

To see transition path comparisons for the ATLAS simulator against the
original:
```
run_atlas_tpaths2D.m
makestats_tpaths2D.m
```

## Simulators

In this code, simulators of SDEs generally have the form

```
function [path_end,X] = simulator(Xzero,m,T,f,dt)
%Simulator for linear SDE, using Euler-Maruyama method
%SDE is of form dx = f(x) dt + dW, X(0) = Xzero

%Inputs:
	%Xzero: initial condition for SDE
	%m: number of simulations to run
	%T: Time limit for simulation
	%f: deterministic portion of SDE above
	%dt: timestep 
%Outputs:	
	%path_end: collection of endpoints of all m simulated paths
	%X: full path of the last simulated
```
The simulators should have the capability of:
	
	1)Running and storing one long trajectory
	2)Running a collection of trajectories and storing each trajectory endpoint

<!--  ## Acknowledgments-->

<!--  * Hat tip to anyone whose code was used
  * Inspiration
  * etc
 --> 

