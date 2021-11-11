# Turing-Model-of-pPf
A Turing model of the transition from laminar to turbulence in wall-bounded shear flows is encoded in the files. It solves a system of 2 PDE with the inclusion of stochastic multiplicative noise. 
The scripts written allows you to :
1) Draw phase portraits of the model
2) Perfom linear instability analysis
3) Explore the kinetics of the model equations in the time domain alone i.e temporal dynamcis without space for the
   - Noise-free system
   - Stochastic system with noise
4) Simulate the full PDE as a nonlinear inital value problem for : 
   - System with OR without nonlinear advection
   - Noise-free or stohchastic
5) A bunch of initial conditions can be constructed inluding pulses, array of pulses, pattern or a chimera of both pulses and patterns within the same domain
6) The model can be tuned by changing the parameters suitably
7) A weakly nonlinear analysis is included to ascertain the nature of the Turing instability : Supercritical / subcritical.

The script is written in Python 3 with PyQt5 rendering the necessary GUI. The following python modules are necessary:
- Dedalus : Click [here](https://dedalus-project.readthedocs.io/en/latest/) to read their documentation for installing Dedalus.
- numpy
- sympy
- scipy
- PyQt5
- matplotlib
- glob
- os

To run the script, copy all the files into a single folder, activate the python environment and then type :
```
python3 Model_Turbulent_transition.py
```
Click [here](https://github.com/PavanVKashyap/Turing-Model-of-pPf/blob/main/user_manual.pdf) for the user manual.
