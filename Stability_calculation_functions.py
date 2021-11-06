# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 00:09:53 2021

@author: pavan
"""
import numpy as np
import sympy as smp
import Configure_Model_functions as cnfg
from scipy import interpolate as intp

def compute_fixed_points(self):
    
    
    ## Compute the coarse grained solutions
    Rsam=np.arange(self.stR1.value(),self.stR2.value(),20)
    r1=[]
    rlam1=[]
    r2=[]
    rlam2=[]
    r3=[]
    rlam3=[]        
    
    
    for i in range(0,len(Rsam)):
        
        r=Rsam[i]
        
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

                    if rlam==0:
                        r1.append(r)
                        rlam1.append(rlam)
                    if rlam>0:
                        r2.append(r)
                        rlam2.append(rlam)
                    
                else:
     
                    r3.append(r)
                    rlam3.append(rlam)
                    
    ## Use bisection to compute the approximate saddle node bifurcation
    r1temp=np.asarray(r1)
    idx=np.argmin(np.abs(r1temp-r2[0]))
    
    r1b=r1[idx-1]
    r2b=r1[idx]
    delta=np.abs(r1b-r2b)
    eps=0.5
    
    while delta>eps:
        
        rbis=(r1b+r2b)/2
        
        r=rbis
        found_sol=False
        
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
                        found_sol=True
                        r2b=rbis
                        
        if found_sol:
            
            delta=np.abs(r1b-r2b)
        
        if not found_sol:
            
            r1b=rbis
            delta=np.abs(r1b-r2b)
            
            
    ## Approximate Saddle node point and its distance to laminar
    self.Rsn=(r1b+r2b)/2
    

    ## Find fine grained values close to Rsn
    Rsam=np.linspace(self.Rsn-2,self.Rsn+2,30)
    
    for i in range(0,len(Rsam)):
        
        r=Rsam[i]
        
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

                    if rlam==0:
                        r1.append(r)
                        rlam1.append(rlam)
                    if rlam>0:
                        r2.append(r)
                        rlam2.append(rlam)
                    
                else:
     
                    r3.append(r)
                    rlam3.append(rlam)
    
    
    #Plot the data
    r1=np.sort(np.asarray(r1))
    r2=np.sort(np.asarray(r2))
    r3=np.sort(np.asarray(r3))
    
    rlam1=np.sort(np.asarray(rlam1))
    rlam2=np.sort(np.asarray(rlam2))
    rlam3=-np.sort(-np.asarray(rlam3))
    
    self.st_axs[0].plot(r1,rlam1,color='b')
    self.st_axs[0].plot(r2,rlam2,color='b')
    self.st_axs[0].plot(r3,rlam3,color='r')
    
    self.st_axs[0].set_xlabel("R")
    self.st_axs[0].set_ylabel(r"$\Delta$")
    
    self.st_fig.tight_layout()
    self.st_canvas.draw()
    
    
####

def stability_criteria_plot(self):
    
    ksam=np.linspace(self.stk1.value(),self.stk2.value(),100)
    
    
    ## Construct the criteria    
    r=self.stR.value()
    
    ## Get the parameters
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
                    
                    pik=self.dm.value()*self.dw.value()*ksam**4 - ksam**2*( (self.dm.value()*gw) + (self.dw.value()*fm) )
                    hk = detJ + pik
                    B = (self.dm.value()+self.dw.value())*ksam**2 - trJ
                    D = B**2 - 4*hk
                    
                    sigma1=np.zeros(len(ksam),dtype=complex)
                    
                    for p in range(0,len(ksam)):
                        
                        if D[p]<0:
                            temp=complex(0,np.sqrt(np.abs(D[p])))
                            sigma1[p]=0.5*(-B[p]+temp)
                                                        
                        if D[p]>0:
                            sigma1[p]=0.5*(-B[p]+np.sqrt(D[p]))
                                                
                    
    self.st_axs[1].plot(ksam,np.real(sigma1),color='orange')
    self.st_axs[1].plot(ksam,np.zeros(len(ksam)),ls='--',color='k')
    self.st_axs[1].set_xlabel("k")
    self.st_axs[1].set_ylabel(r"$Re(\sigma_1)$")
    self.st_fig.tight_layout()
    self.st_canvas.draw()
    
    
######

def find_turing(self):
    
    
    ## construct hmin(R)
    Rsam=np.linspace(self.Rsn+2,self.stR2.value(),30)
    ksam=np.linspace(self.stk1.value(),self.stk2.value(),100)
    hmin=np.zeros(len(Rsam))
    
    for i in range(0,len(Rsam)):
        
        r=Rsam[i]
        
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
                        
                        pik=self.dm.value()*self.dw.value()*ksam**4 - ksam**2*( (self.dm.value()*gw) + (self.dw.value()*fm) )
                        hk = detJ + pik
                        hmin[i]=hk.min()
        
        
        ## Interpolate to find the zero cross over
        intpf=intp.interp1d(hmin,Rsam)
        Rc=float(intpf(0))
        
        ## Set the slider and plot
        self.stR.setValue(Rc)
       
        
            