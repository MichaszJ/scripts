# Scripts

This repo is a collection of various scripts, mainly related to my schoolwork in aerospace engineering. Examples on how to use the provided scripts can be found in the `Usage` Jupyter Notebooks in each folder.

# Showcase
## Orbital Dynamics - N-Body Problem
The various orbital dynamics scripts that solve the two-body, circular restricted three-body, and general three-body all use the `rk_45` ODE solver (see [solve_ode_ivp.py](https://github.com/MichaszJ/scripts/blob/main/Numerical-Methods/solve_ode_ivp.py)). 

Two-body problem solution using `two_body_propagator` function:
![](images/animation1.gif)

Three-body problem solution using `three_body_propagator` function:

![](images/animation2.gif)

Lunar transfer and circularization simulation, made using a modified version of the `three_body_cr_propagator`, which included thrust terms to simulate the spacecraft firing its engines:

![](images/animation3.gif)