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

The script is written in Python 3 with PyQt5 rendering the necessary GUI.

<strong>Note: This can execute only on a Linux or Mac operating system</strong>.

The simulation needs python3 with the following modules:
- Anaconda python distribution
   - Please follow the installation instruction [here](https://www.anaconda.com/download) for the specific OS
- Dedalus v2: Follow the steps below for installing Dedalus v2 after Anaconda is installed
   - Activate the base anaconda environment
   ```
   conda activate
   ```
   - Create a python environment for Dedalus v2
   ```
   conda create -n dedalus2
   conda activate dedalus2
   ```
   - If you have MacOS run the following command
   ```
   conda config --env --set subdir osx-64
   ```
   - Disable multi-threading
   ```
   conda env config vars set OMP_NUM_THREADS=1
   conda env config vars set NUMEXPR_MAX_THREADS=1
   ```
   - Install Dedalus v2 from conda-forge
   ```
   conda install -c conda-forge dedalus=2.2207.3
   ```
   - Activate the environment
   ```
   conda activate dedalus2
   ```
- Sympy module
   - Execute the following commands in the dedalus2 environment created
   ```
   conda install sympy
   ```
- PyQt5 module
   - Execute the following command for installing PyQt5
   ```
   conda install pyqt
   ```
Once these modules are insatlled. You are set to go !!

To run the script, copy all the files into a single folder, activate the python environment and then type :
```
python3 Model_Turbulent_transition.py
```
Click [here](https://github.com/PavanVKashyap/Turing-Model-of-pPf/blob/main/user_manual.pdf) for the user manual.
