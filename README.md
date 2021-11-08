# Turing-Model-of-pPf
A Turing model of the transition from laminar to turbulence in wall-bounded shear flows. It solves a system of 2 PDE with the inclution of stastic multiplicative noise. 
The scripts written allows you to :
1) Draw phase portraits of the model
2) Explore the kinetics of the model equations in the time domain alone i.e temporal dynamcis without space 
   - Noise-free
   - Stochastic system with noise
3) Simulate the full PDE as a a nonlinear inital value problem with options : 
   - System with OR without nonlinear advection
   - Noise-free or stohchastic
4) A bunch of initial conditions can be constructed inluding pulses, array of pulses, pattern or a chimera of both pulses and patterns with the same domain
5) The model can be tuned by changing the parameters suitably
6) A weakly nonlinear analysis is included to ascertaint he nature of the Turing instability : Supercritical / subcritical.

The script is written in Python 3 with PyQt5 rendering the necessary GUI. 
The python modules required are available in the requirements.txt file which can be used to setup a python environment. 
The script utilizes the module Dedalus to solve the PDE. Click [here](https://dedalus-project.readthedocs.io/en/latest/) to read their documentation for installing Dedalus.
To run the script, copy all the files into a single folder, activate the python environment and then type :
'''
python3 Model_Turbulent_transition.py
'''
