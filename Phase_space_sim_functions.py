# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 00:44:48 2021

@author: pavan
"""
import numpy as np
import sympy as smp
import Phase_plane_analysis_functions as ppf
import Configure_Model_functions as cnfg


### Define the right hand sides
def RHS(self,r,m,w):
    
    alpha_m=cnfg.prm(self,r)[0]
    beta_m=cnfg.prm(self,r)[1]
    sigma_m=cnfg.prm(self,r)[2]
    alpha_w=cnfg.prm(self,r)[3]
    beta_w=cnfg.prm(self,r)[4]
    gamma_w=cnfg.prm(self,r)[5]
    Mbar=cnfg.prm(self,r)[6]
    M0=cnfg.prm(self,r)[7]
    
    dydt=np.zeros(2)    
    dydt[0]=alpha_m*(M0 - m) - beta_m*(m-Mbar)*w**4 + sigma_m*w**2
    dydt[1]=-alpha_w*w + beta_w*(m-Mbar)*w**3 - gamma_w*w**3 - sigma_m*m*w

    return dydt  

### Compute the noise free fixed point
def compute_fixed_points(self,r):
    
   
    ## Get teh parameters
    alpha_m=cnfg.prm(self,r)[0]
    beta_m=cnfg.prm(self,r)[1]
    sigma_m=cnfg.prm(self,r)[2]
    alpha_w=cnfg.prm(self,r)[3]
    beta_w=cnfg.prm(self,r)[4]
    gamma_w=cnfg.prm(self,r)[5]
    Mbar=cnfg.prm(self,r)[6]
    M0=cnfg.prm(self,r)[7]
    
    #declare the symbols
    m,w=smp.symbols('m w', negative=False)
    
    ## Write the equations
    f = alpha_m*(M0 - m) - beta_m*(m-Mbar)*w**4 + sigma_m*w**2
    g = -alpha_w*w + beta_w*(m-Mbar)*w**3 - gamma_w*w**3 - sigma_m*m*w
    
    # Compute the solutions
    sol=smp.nonlinsolve([f,g], [m,w])
    
    ## Sort through the solutions to get the turbulent fixed point
    for k in range(0,len(sol)):
        
        if complex(sol.args[k][0]).imag==0 and complex(sol.args[k][0]).real >= 0 and complex(sol.args[k][1]).imag==0 and complex(sol.args[k][1]).real >= 0:
            
            ms=complex(sol.args[k][0]).real
            ws=complex(sol.args[k][1]).real
            
            ## Construct the Jacobian J
            fm = -alpha_m - beta_m*ws**4
            fw = -4*beta_m*(ms-Mbar)*ws**3 + 2*sigma_m*ws
            gm = beta_w*ws**3 - sigma_m*ws
            gw = -alpha_w + 3*beta_w*(ms-Mbar)*ws**2 - 3*gamma_w*ws**2 - sigma_m*ms
            
            trJ=fm+gw
            detJ=(fm*gw)-(fw*gm)
            
            rlam=np.sqrt((ms-1)**2+(ws-0)**2)
                
            if trJ<0 and detJ>0:

                if rlam>0:
                    Mt=ms
                    Wt=ws
                    
                    
    return Mt,Wt
    

### Runge kutta integration for noise free simulation
def Runge_Kutta(self,i):
    
   
    m0=self.Mph[i]
    w0=self.Wph[i]
    
    ## Time loop
    k1=self.dt_ph_E.value()*RHS(self,self.R_ph_E.value(),m0,w0)
    k2=self.dt_ph_E.value()*RHS(self,self.R_ph_E.value(),m0+k1[0]/2,w0+k1[1]/2)
    k3=self.dt_ph_E.value()*RHS(self,self.R_ph_E.value(),m0+k2[0]/2,w0+k2[1]/2)
    k4=self.dt_ph_E.value()*RHS(self,self.R_ph_E.value(),m0+k3[0],w0+k3[1])
    
    self.Mph.append( m0 + ( (self.dt_ph_E.value()/6) * (k1[0] + 2*k2[0] +2*k3[0] + k4[0]) ) )
    self.Wph.append( w0 + ( (self.dt_ph_E.value()/6) * (k1[1] + 2*k2[1] +2*k3[1] + k4[1]) ) )
    
    
### Euler Maruyama time step for simulation with Noise
def Euler_Maruyama(self,i):
        
    m0=self.Mph[i]
    w0=self.Wph[i]
    
    nst=w0*self.NS_ph_E.value()*np.random.normal(loc=0,scale=np.sqrt(self.dt_ph_E.value()))
    
    self.Mph.append( m0 + RHS(self,self.R_ph_E.value(),m0,w0)[0]*self.dt_ph_E.value() )
    self.Wph.append( w0 + RHS(self,self.R_ph_E.value(),m0,w0)[1]*self.dt_ph_E.value() + nst )
    



### Plot the results of the full trajectory simulations
def plot_results(self,y):
    
    ## Control parameter
    r=float(self.R_ph_E.value())
    npts=int(self.Np_pp_E.value())
    
        
    #Plot the null clines
    m_wnull, wsam_w, m_mnull, wsam_m = ppf.plot_parameters(self,r,npts)[5:]
    self.Ph_phaxs.plot(m_wnull,wsam_w,color='r',label=r'$W_n$')
    self.Ph_phaxs.plot(m_mnull,wsam_m,color='b',label=r'$M_n$')
    self.Ph_phaxs.set_xlabel("M")
    self.Ph_phaxs.set_ylabel("W")
    self.Ph_phaxs.set_xlim([0,self.PP_Mmax.value()])
    self.Ph_phaxs.set_ylim([0,self.PP_Wmax.value()])
    
    ## Plot the results
    ## Plot M(t)
    self.Ph_axs[0,0].plot(np.arange(0,int(self.iter_ph_E.value())+1),y[0,:],color='blue')
    self.Ph_axs[0,0].set_xlabel("Iterations")
    self.Ph_axs[0,0].set_ylabel("M")
    self.Ph_axs[0,0].set_ylim([0,self.PP_Mmax.value()])
    
    ## Plot W(t)
    self.Ph_axs[1,0].plot(np.arange(0,int(self.iter_ph_E.value())+1),y[1,:],color='orange')
    self.Ph_axs[1,0].set_xlabel("Iterations")
    self.Ph_axs[1,0].set_ylabel("W")
    self.Ph_axs[1,0].set_ylim([0,self.PP_Wmax.value()])
    
    ## Plot the fixed point
    if self.R_ph_E.value()>self.Rsn:
        Ms,Ws=compute_fixed_points(self,r)
    else:
        Ms=1.0
        Ws=0.0
    self.Ph_axs[0,0].plot(np.arange(0,int(self.iter_ph_E.value())+1),np.ones(int(self.iter_ph_E.value())+1)*Ms,color='k',lw=0.5)
    self.Ph_axs[1,0].plot(np.arange(0,int(self.iter_ph_E.value())+1),np.ones(int(self.iter_ph_E.value())+1)*Ws,color='k',lw=0.5)
        
    ## Plot trajectory in phase space
    self.Ph_phaxs.scatter(y[0,:],y[1,:],c=np.arange(int(self.iter_ph_E.value())+1),cmap='jet',vmin=0,vmax=int(self.iter_ph_E.value())+1,s=10,marker='X',zorder=0)
    self.Ph_phaxs.set_xlabel("M")
    self.Ph_phaxs.set_ylabel("W")

    self.Ph_fig.tight_layout()
    self.Ph_canvas.draw()
    
    
### Runge kutta integration for noise free simulation
def Runge_Kutta_full_traj(self):
    
    ## Control parameter
    r=float(self.R_ph_E.value())
    
    ## Define the functions
    dt=self.dt_ph_E.value()
    y=np.zeros([2,int(self.iter_ph_E.value())+1])
    
    ## Initial conditions
    y[0,0]=self.M0_ph_E.value()
    y[1,0]=self.W0_ph_E.value()
     
    ## Time loop
    for i in range(0,int(self.iter_ph_E.value())):
        
        k1=dt*RHS(self,r,y[0,i],y[1,i])
        k2=dt*RHS(self,r,y[0,i]+k1[0]/2,y[1,i]+k1[1]/2)
        k3=dt*RHS(self,r,y[0,i]+k2[0]/2,y[1,i]+k2[1]/2)
        k4=dt*RHS(self,r,y[0,i]+k3[0],y[1,i]+k3[1])
        
        y[:,i+1] = y[:,i] + ( (dt/6) * (k1 + 2*k2 +2*k3 + k4) )
        
        
    ## Plot the data if evolution is to be seen
    plot_results(self,y)
    
    
### Euler Maruyama time step for simulation with Noise
def Euler_Maruyama_full_traj(self):
    
    ## Control parameter
    r=float(self.R_ph_E.value())
    
    ## Define the functions
    dt=self.dt_ph_E.value()
    y=np.zeros([2,int(self.iter_ph_E.value())+1])
    
    ## Initial conditions
    y[0,0]=self.M0_ph_E.value()
    y[1,0]=self.W0_ph_E.value()
    
    ## Noise term
    noise_amp=self.NS_ph_E.value()
    nst=np.zeros(2)
    
    ## Actual time loop
    for i in range(0,int(self.iter_ph_E.value())):
        
        nst[1]=y[1,i]*noise_amp*np.random.normal(loc=0,scale=np.sqrt(dt))
        
        y[:,i+1]=y[:,i]+RHS(self,r,y[0,i],y[1,i])*dt + nst
    
    ## Plot the results
    plot_results(self,y)

    