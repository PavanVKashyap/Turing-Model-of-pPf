import numpy as np
import dedalus.public as de
import h5py as hp
import Configure_Model_functions as cnfg

def solver_initialize(self):
    
    #Parameters for simulation
    Lx=self.Lx_sim_E.value();
    Nx=self.Nx_sim_E.value();
    dt=self.dt_sim_E.value();

    #Control Parameter
    R=self.R_sim_E.value()
    
    #Waleffe model parameters
    wlm = 2.47
    wlu = 5.2
    wlv = 7.67
    wlw = 7.13
    
    wsm = 0.31
    wsu = 1.29
    wsv = 0.22
    wsw = 0.68
    
    M0 = 1
    
    # Waleffe Parameter relations
    lm = wlm
    bm = (wsu**2 * wsv**2) / (wlu * wlv**2)
    sm = wsm
    
    lw = wlw
    bw = (wsu * wsv * wsw) / (wlu * wlv) 
    gw = wsv**2 / wlv
    
    M1 = ( wlv * wsw) / (wsu * wsv)
    
    # Scaling factors
    f1=self.f1.value()
    f2=self.f2.value()
    f3=self.f3.value()
    f4=self.f4.value()
    f5=self.f5.value()
    f6=self.f6.value()
    
    
    # Model parameters
    alpha_m = f1 * lm / R
    beta_m = f2 * bm * R**3
    sigma_m = f3 * sm
    
    alpha_w = f4 * lw / R
    beta_w = f5 * bw * R**2
    gamma_w = f6 * gw * R
    
    Mbar = M1 / R
    
    
    #Define the basis
    xbasis=de.Fourier('x',Nx,interval=(0,Lx),dealias=3/2)
    
    #Create the domain
    domain=de.Domain([xbasis],grid_dtype=np.float64)
    
    #define the problem
    problem=de.IVP(domain, variables=['m','mx','w','wx'])
    problem.parameters['m0']=M0;
    problem.parameters['mb']=Mbar;
    
    problem.parameters['am']=alpha_m
    problem.parameters['bm']=beta_m
    problem.parameters['sm']=sigma_m
    problem.parameters['Dm']=self.dm.value()
    
    problem.parameters['aw']=alpha_w
    problem.parameters['bw']=beta_w
    problem.parameters['gw']=gamma_w
    problem.parameters['Dw']=self.dw.value()
        
    if self.NS_sim_check.isChecked() and not self.adv_sim_check.isChecked() :
            
        #Define the noise term
        noise=domain.new_field()
        noise.set_scales(1)
        noise['g']=np.random.normal(loc=0.0,scale=1.0,size=Nx)
        problem.parameters['N']=noise
        problem.parameters['S']=self.NS_sim_E.value()
        
        #Define the equations
        problem.add_equation("dt(m) - Dm*dx(mx) = -m*mx + am*(m0-m) - bm*(m-mb)*w**4 + sm*w**2")
        problem.add_equation("mx - dx(m) = 0")
        problem.add_equation(" dt(w) - Dw*dx(wx) = -m*wx - aw*w + bw*(m-mb)*w**3 - gw*w**3 - sm*m*w + S*N*w")
        problem.add_equation("wx-dx(w) = 0")
    
    if self.adv_sim_check.isChecked() and not self.NS_sim_check.isChecked() :
            
        
        problem.add_equation("dt(m) - Dm*dx(mx) = am*(m0-m) - bm*(m-mb)*w**4 + sm*w**2")
        problem.add_equation("mx - dx(m) = 0")
        problem.add_equation(" dt(w) - Dw*dx(wx) = - aw*w + bw*(m-mb)*w**3 - gw*w**3 - sm*m*w")
        problem.add_equation("wx-dx(w) = 0")
        
    if self.NS_sim_check.isChecked() and self.adv_sim_check.isChecked():

        #Define the noise term
        noise=domain.new_field()
        noise.set_scales(1)
        noise['g']=np.random.normal(loc=0.0,scale=1.0,size=Nx)
        problem.parameters['N']=noise
        problem.parameters['S']=self.NS_sim_E.value()
        
        #Define the equations
        problem.add_equation("dt(m) - Dm*dx(mx) = am*(m0-m) - bm*(m-mb)*w**4 + sm*w**2")
        problem.add_equation("mx - dx(m) = 0")
        problem.add_equation(" dt(w) - Dw*dx(wx) = - aw*w + bw*(m-mb)*w**3 - gw*w**3 - sm*m*w + S*N*w")
        problem.add_equation("wx-dx(w) = 0")
        
    if (not self.NS_sim_check.isChecked()) and (not self.adv_sim_check.isChecked()):
        
        problem.add_equation("dt(m) - Dm*dx(mx) = -m*mx + am*(m0-m) - bm*(m-mb)*w**4 + sm*w**2")
        problem.add_equation("mx - dx(m) = 0")
        problem.add_equation(" dt(w) - Dw*dx(wx) = -m*wx - aw*w + bw*(m-mb)*w**3 - gw*w**3 - sm*m*w")
        problem.add_equation("wx-dx(w) = 0")
        

    #Build the solver
    solver=problem.build_solver(de.timesteppers.SBDF1)
    
    #Initial condition
    txt=self.IC_select.currentText()
    
    if txt=='Noise':
        
        self.Msim,self.Wsim=cnfg.noise_plot(self)
        
    if txt=='Pulse':
        
        self.Msim,self.Wsim=cnfg.pulse_plot(self)
        
    if txt=='Pattern':
        
        self.Msim,self.Wsim=cnfg.pattern_plot(self)
        
    if txt=='Chimera':
        
        self.Msim,self.Wsim=cnfg.chimera_plot(self)
        
    ## decalre the variable
    self.m = solver.state['m']
    mx = solver.state['mx']
    self.w = solver.state['w']
    wx = solver.state['wx']
    
    #Set the initial condition
    self.m['g'] = self.Msim
    self.m.differentiate(0, out=mx)
    self.w['g'] = self.Wsim
    self.w.differentiate(0, out=wx)
    
    ## Save the data
    self.svst=solver.evaluator.add_file_handler("case", sim_dt=self.spt_dt.value())
    self.svst.add_system(solver.state)

    return solver

def simulation_time_step(self):
    
    ## Simulate for 10 iterations befre live update
    for k in range(0,10):
            
        if self.NS_sim_check.isChecked():
                
            #Update the noise term
            ns=np.random.normal(loc=0.0,scale=np.sqrt(self.dt_sim_E.value()),size=self.Nx_sim_E.value()-1)
            ns=np.block([ns,ns[0]])
            self.solver.problem.namespace['N'].set_scales(1)
            self.solver.problem.namespace['N'].value=ns
            
        #Time stepper
        self.solver.step(self.dt_sim_E.value())
        
    
    ## Update the values
    self.m.set_scales(1)
    self.w.set_scales(1)
    
    self.Msim=self.m['g']
    self.Wsim=self.w['g']
    
def spt_plot(self):

    file=hp.File("case/case_s1/case_s1_p0.h5",'r')
    
    x=file['scales']['x']['1.0'][:]
    t=file['scales']['sim_time'][:]
    m=file['tasks']['m'][:]
    w=file['tasks']['w'][:]
    
    self.spt_window.sptaxs[0].cla()
    self.spt_window.sptaxs[1].cla()
    
    self.spt_window.sptaxs[0].pcolormesh(x,t,m,shading='gouraud',vmin=0,vmax=1,cmap='jet')
    self.spt_window.sptaxs[1].pcolormesh(x,t,w,shading='gouraud',vmin=0,vmax=0.1,cmap='jet')
    
    self.spt_window.sptaxs[0].set_title("M")
    self.spt_window.sptaxs[0].set_xlabel("x")
    self.spt_window.sptaxs[0].set_ylabel("Time")
    
    self.spt_window.sptaxs[1].set_title("W")
    self.spt_window.sptaxs[1].set_xlabel("x")
    self.spt_window.sptaxs[1].set_ylabel("Time")
    
    self.spt_window.sptfig.tight_layout()
    self.spt_window.spt_canvas.draw()
    
    
        
    
        