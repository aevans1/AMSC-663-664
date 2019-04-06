# ATLAS readme:
## Work in Progress!

This repository is an implementation of the ATLAS algorithm from Crosskey and
Maggioni's paper [ATLAS: A Geometric Approach to Learning High-Dimensional
Stochastic Systems Near Manifolds](https://epubs.siam.org/doi/abs/10.1137/140970951). For a more in depth description of the code, see the [User Guide](doc/ATLAS_UserGuide.pdf).

## Prerequisites

The codes are all written in MATLAB, and running the code just requires downloading
the src directory and adding the src directory to your matlab path.

## Overview

The main files to run the code are:

**atlas_driver.m** - choose from parameters, set up examples from Crosskey,
Maggioni paper

**construction.m** - construct the ATLAS for the given example

**learned_simulator.m** - Euler-Maruyama simulator for ATLAS learned simulator

The current files for testing output are:

**run_atlas_tpaths2d.m** - (For 2D potential example below) Computes average
transition times between the potential wells of the original simulator and the
ATLAS simulator

**makestats_tpaths2d.m** - (For 2D potential example below) Plots the bar graph of
transition averages from run_atlas_2paths2d.m

## Running ATLAS

### Parameters

#### Problem-Specific inputs

The user should have an input SDE of form dx = f(x) dt + dW in R^D with a known timestep and metric 'rho' for the ambient space. In particular, the user must input:

* **f** - function with inputs as column vectors in R^D, outputs as column vectors in
R^D

* **rho** - metric in the space of th input simulator. Example: euclidean distance, rho = @(x,y) norm(x -
y)

* **d** - intrinsic dimension of the manifold. All ATLAS outputs will be in R^d.

* **dt_sim** - original time-step of the input SDE


#### ATLAS inputs

* **delta** - spatial homogenization scale. ATLAS is more accurate and more
computationally intensive for smaller delta values.

* **t_0** - simulation time for the short paths run by ATLAS. Should be chosen such
trajectories running the original simulator for time t_0 are delta away from
the starting location on average.

* **m** - the number of landmarks for each net point. m should be >= d, accuracy and
computational complexity increase as m increases.

* **p** - the number of sample paths for each net point, should be O(delta^-4).

* **dt** - time step for ATLAS simulator, should be O(delta/ln(1/delta)).

### 2D, 3-Well Potential example

The line
```
atlas_driver(3);
```
by default:

* sets ATLAS parameters to the Crosskey, Maggioni parameters for ex. 5.3.1 in
the paper (see page 31 of Crosskey and Maggioni's [ATLAS paper](https://epubs.siam.org/doi/abs/10.1137/140970951).

* constructs the delta-net

* runs construction.m, which constructs the ATLAS and local SDE simulators

To see transition path comparisons for the ATLAS simulator against the
original, run:

```
run_atlas_tpaths2D.m
makestats_tpaths2D.m
```

The file run_atlas_tpaths2d simulates one long trajectory for the original
simulator, tracks transition paths between the potential wells, and then does
the same for the learned atlas simluator. Makestats generates the bargraphs
with the comparison of the average transition times.


The file atlas_driver.m contains two other examples from the [ATLAS paper](https://epubs.siam.org/doi/abs/10.1137/140970951), which correspond to examples 5.2.1,5.2.2.

<!--  ## Acknowledgments-->

<!--  * Hat tip to anyone whose code was used
  * Inspiration
  * etc
 --> 

