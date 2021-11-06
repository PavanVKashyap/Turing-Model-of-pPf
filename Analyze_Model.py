# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 11:48:40 2021

@author: pavan
"""

import numpy as np
import sympy as smp
from scipy.interpolate import interp1d as intp

## Function for getting parameters
def prm(R,f) :
    
    #Control Parameter
    R=R
    
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
    f1=f[0]#5.0
    f2=f[1]#5.0
    f3=f[2]#0.2
    f4=f[3]#0.2
    f5=f[4]#0.2
    f6=f[5]#5.0
    
    
    # Model parameters
    alpha_m = f1 * lm / R
    beta_m = f2 * bm * R**3
    sigma_m = f3 * sm
    
    alpha_w = f4 * lw / R
    beta_w = f5 * bw * R**2
    gamma_w = f6 * gw * R
    
    Mbar = M1 / R

    return np.asarray([alpha_m,beta_m,sigma_m,alpha_w,beta_w,gamma_w,Mbar,M0])



def plot_nullcline(self):
    
    r=self.cR.value()
    f=np.asarray([self.cf1.value(),self.cf2.value(),self.cf3.value(),
                  self.cf4.value(),self.cf5.value(),self.cf6.value()])
    
    ## Get the model parameters
    alpha_m=prm(r,f)[0]
    beta_m=prm(r,f)[1]
    sigma_m=prm(r,f)[2]
    alpha_w=prm(r,f)[3]
    beta_w=prm(r,f)[4]
    gamma_w=prm(r,f)[5]
    Mbar=prm(r,f)[6]
    M0=prm(r,f)[7]
    
    ## Sample space for nullclines
    wsam_max=1
    nsam=10**5
    wsam_m=np.linspace(0,wsam_max,nsam)
    wsam_w=np.linspace(1.01*np.sqrt(sigma_m/beta_w),wsam_max,nsam)
    
    ## Nullclines coordinates
    m_mnull=( (sigma_m*wsam_m**2) + (beta_m * Mbar * wsam_m**4) + (alpha_m * M0) ) / (alpha_m + (beta_m*wsam_m**4) ) 
    m_wnull=( alpha_w + (gamma_w*wsam_w**2) + (beta_w*Mbar*wsam_w**2)  ) / ( (beta_w*wsam_w**2) - sigma_m) 
    
    #Plot the null clines
    self.anphxs.cla()
    self.anphxs.plot(m_wnull,wsam_w,color='r',label=r'$W_n$')
    self.anphxs.plot(m_mnull,wsam_m,color='b',label=r'$M_n$')
    self.anphxs.set_title(f"Nullclines for R={r}")
    self.anphxs.set_xlabel("M")
    self.anphxs.set_ylabel("W")
    self.anphxs.set_xlim([self.cMmin.value(),self.cMmax.value()])
    self.anphxs.set_ylim([self.cWmin.value(),self.cWmax.value()])
    self.anphxs.legend(loc=1)

    self.anfig.tight_layout()
    self.ancanvas.draw()


### Categorize the stability conditions
def bifurcation_criteria(self):
    
    ## Parameters
    fp=np.asarray([self.cf1.value(),self.cf2.value(),self.cf3.value(),
                  self.cf4.value(),self.cf5.value(),self.cf6.value()])
    
    d=self.cdw.value()/self.cdm.value()
    
    self.anxs[0,1].cla()
    self.anxs[0,2].cla()
    self.anxs[1,1].cla()
    self.anxs[1,2].cla()
    
    ## Find the saddle node
    r1=50
    r2=900
    eps=0.5
    delta=np.abs(r2-r1)
    
    ## Bisection
    while delta>eps:
        
        rbis=(r1+r2)/2
        
        r=rbis
        
        ## Get the parameters
        alpha_m=prm(r,fp)[0]
        beta_m=prm(r,fp)[1]
        sigma_m=prm(r,fp)[2]
        alpha_w=prm(r,fp)[3]
        beta_w=prm(r,fp)[4]
        gamma_w=prm(r,fp)[5]
        Mbar=prm(r,fp)[6]
        M0=prm(r,fp)[7]
        
        #declare the symbols
        m,w=smp.symbols('m w', negative=False)
        
        ## Write the equations
        f = alpha_m*(M0 - m) - beta_m*(m-Mbar)*w**4 + sigma_m*w**2
        g = -alpha_w*w + beta_w*(m-Mbar)*w**3 - gamma_w*w**3 - sigma_m*m*w
        
        # Compute the solutions
        sol=smp.nonlinsolve([f,g], [m,w])
        
        ## Find number of solutions
        count=0
        for k in range(0,len(sol)):
            
            if complex(sol.args[k][0]).imag==0 and complex(sol.args[k][0]).real >= 0 and complex(sol.args[k][1]).imag==0 and complex(sol.args[k][1]).real >= 0:
                
                ms=complex(sol.args[k][0]).real
                ws=complex(sol.args[k][1]).real
                
                count+=1    
        
        
        if count==3:    
            r2=rbis
            delta=np.abs(r1-r2)
            
        if count==1:
            r1=rbis
            delta=np.abs(r1-r2)
            
    
    ## Saddle node bifurcation
    rsn=(r1+r2)/2
    self.summary.clear()
     
    ## Print on screen
    self.summary.append("Saddle node Bifurcation found")
    self.summary.append(f"Rsn={rsn:0.2f}")

    
    ## Sample the upper branch
    nsamp=40
    rsam=np.linspace(rsn+5,rsn+500,nsamp)
    hmin=np.zeros(len(rsam))
    
    ## Sort the upper branch solution and determine stability criteria
    s1=np.zeros([3,nsamp])
    s2=np.zeros([3,nsamp])
    t1=np.zeros([3,nsamp])
    t2=np.zeros([3,nsamp])
    
    for i in range(0,len(rsam)):
        
        r=rsam[i]
        
        ## Get the model parameters
        alpha_m=prm(r,fp)[0]
        beta_m=prm(r,fp)[1]
        sigma_m=prm(r,fp)[2]
        alpha_w=prm(r,fp)[3]
        beta_w=prm(r,fp)[4]
        gamma_w=prm(r,fp)[5]
        Mbar=prm(r,fp)[6]
        M0=prm(r,fp)[7]
        
        #declare the symbols
        m,w=smp.symbols('m w', negative=False)
        
        ## Write the equations
        f = alpha_m*(M0 - m) - beta_m*(m-Mbar)*w**4 + sigma_m*w**2
        g = -alpha_w*w + beta_w*(m-Mbar)*w**3 - gamma_w*w**3 - sigma_m*m*w
        
        # Compute the solutions
        sol=smp.nonlinsolve([f,g], [m,w])
        
        ## Variables to record the derivatives
        fm=np.zeros(3)
        fw=np.zeros(3)
        gm=np.zeros(3)
        gw=np.zeros(3)
        rlam=np.zeros(3)
        count=0
        
        ## Record the derivatives at all the solutions
        for k in range(0,len(sol)):
            
            if complex(sol.args[k][0]).imag==0 and complex(sol.args[k][0]).real >= 0 and complex(sol.args[k][1]).imag==0 and complex(sol.args[k][1]).real >= 0:
                
    
                ms=complex(sol.args[k][0]).real
                ws=complex(sol.args[k][1]).real
                
                ## Construct the Jacobian J
                fm[count] = -alpha_m - beta_m*ws**4
                fw[count] = -4*beta_m*(ms-Mbar)*ws**3 + 2*sigma_m*ws
                gm[count] = beta_w*ws**3 - sigma_m*ws
                gw[count] = -alpha_w + 3*beta_w*(ms-Mbar)*ws**2 - 3*gamma_w*ws**2 - sigma_m*ms
                
                rlam[count]=np.sqrt((ms-1)**2+(ws-0)**2)
                
                count+=1

                              
        ## Sort through the solutions
        temp=np.argsort(rlam)
        fm=fm[temp]
        fw=fw[temp]
        gm=gm[temp]
        gw=gw[temp]
                
        for k in range(0,3):
            
            s1[k,i]=fm[k]+gw[k]
            s2[k,i]=fm[k]*gw[k]-fw[k]*gm[k]
            
            t1[k,i]=d*fm[k]+gw[k]
            t2[k,i]=(d*fm[k]+gw[k])**2 - 4*d*(fm[k]*gw[k]-fw[k]*gm[k])
            
            if k==2:
                
                hmin[i]=(fm[k]*gw[k]-fw[k]*gm[k]) - (d*fm[k]+gw[k])**2/(4*d)
            
            
    ### Plot the stability conditions
    lbs=["L","U","T"]
    lst=['-','--','-']
    for k in range(0,3):
        
        self.anxs[0,1].plot(rsam,s1[k,:],ls=lst[k],label=lbs[k])
        self.anxs[0,2].plot(rsam,s2[k,:],ls=lst[k],label=lbs[k])
    
    self.anxs[0,1].plot(rsam,np.zeros(len(rsam)),ls='dashdot',color='k')
    self.anxs[0,1].set_xlabel("R")
    self.anxs[0,1].set_title("S1<0")
    #self.anxs[0,1].legend(loc="best",ncol=3)
    
    self.anxs[0,2].plot(rsam,np.zeros(len(rsam)),ls='dashdot',color='k')
    self.anxs[0,2].set_xlabel("R")
    self.anxs[0,2].set_title("S2>0")
    #self.anxs[0,2].legend(loc="best",ncol=3)
    
    ## Plot the turing criteria
    for k in range(0,3):
        
        self.anxs[1,1].plot(rsam,t1[k,:],ls=lst[k],label=lbs[k])
        self.anxs[1,2].plot(rsam,t2[k,:],ls=lst[k],label=lbs[k])
    
    self.anxs[1,1].plot(rsam,np.zeros(len(rsam)),ls='dashdot',color='k')
    self.anxs[1,1].set_xlabel("R")
    self.anxs[1,1].set_title("T1>0")
    #self.anxs[1,1].legend(loc="best",ncol=3)
    
    self.anxs[1,2].plot(rsam,np.zeros(len(rsam)),ls='dashdot',color='k')
    self.anxs[1,2].set_xlabel("R")
    self.anxs[1,2].set_title("T2>0")
    #self.anxs[1,2].legend(loc="best",ncol=3)
    
    self.anfig.tight_layout()
    self.ancanvas.draw()
    
    
    ## Test for Hopf instability
    hp1=s1[2,:][s1[2,:]>0]
    hp2=s2[2,:][s2[2,:]<0]
    if (len(hp1)==0) and (len(hp2)==0):
        
        ## Hopf bifurcation not found
        self.summary.append("\nNo Hopf bifurcation found")
        
    if (len(hp1)!=0) or (len(hp2)!=0):
        
        rh=[]

        if (len(hp1)!=0):
            
            intpf=intp(s1[2,:],rsam,kind='cubic')
            rh.append(intpf(0))
                
                
        if (len(hp2)!=0):
            
            intpf=intp(s2[2,:],rsam,kind='cubic')
            rh.append(intpf(0))
           
        
        ## Hopf bifurcation found
        self.summary.append("\nHopf bifurcation found")
        self.summary.append(f"Rh={np.min(np.asarray(rh)):0.2f}")
    
    
    #Test for Turing
    Tst1=np.sort(np.where(np.sign(t1[2,:])==1)[0])
    Tst2=np.sort(np.where(np.sign(t2[2,:])==1)[0])
    
    l1=np.min([len(Tst1),len(Tst2)])
    
    if l1 > 0:
            
        test=np.abs(Tst1[:l1]-Tst2[:l1]).min()
        
        if test==0:
            
            ## Find the critical point
            intpf=intp(hmin,rsam,kind='cubic')
            rc=intpf(0)
            
            ## Draw the nullclines at the critical point
            self.cR.setValue(rc)
            plot_nullcline(self)
            
            ## Add the details in the summary
            self.summary.append("\nTuring instability found")
            self.summary.append(f"Rt={rc:0.2f}")
            
            ## check for subcriticality
            eta,g,kc=weakly_nonlinear(self,rc)
            if g>0:
                self.summary.append(f"kc={kc:0.2f}")    
                self.summary.append("\nTuring instability : Subcritical")
                self.summary.append(f"eta={eta:0.2e}")
                self.summary.append(f"g={g:0.2e}")
            if g<0:
                self.summary.append(f"kc={kc:0.2f}")
                self.summary.append("\nTuring instability : Supercritical")
                self.summary.append(f"eta={eta:0.2e}")
                self.summary.append(f"g={g:0.2e}")
                
        else:
            self.summary.append("No Turing instability found")
                
    else:
        self.summary.append("No Turing instability found")
        
        
    ### Define the line segments
    self.summary.append("\nPlot legend :")
    self.summary.append("Blue : Laminar")
    self.summary.append("Green : Turbulent")
    self.summary.append("Orange : Unstable")
    


## Computet the dierviative of the Jacobian at Rc
def DJ(self,rc):
    
    dr=1.0
    
    rsam=[rc+dr,rc-dr]
    
    ## Parameters
    fp=np.asarray([self.cf1.value(),self.cf2.value(),self.cf3.value(),
                  self.cf4.value(),self.cf5.value(),self.cf6.value()])
    
    ## Variable to hold the values
    J=np.zeros([2,2,len(rsam)])

    for i in range(0,len(rsam)):
        
        r=rsam[i]
        
        ## Get the model parameters
        alpha_m=prm(r,fp)[0]
        beta_m=prm(r,fp)[1]
        sigma_m=prm(r,fp)[2]
        alpha_w=prm(r,fp)[3]
        beta_w=prm(r,fp)[4]
        gamma_w=prm(r,fp)[5]
        Mbar=prm(r,fp)[6]
        M0=prm(r,fp)[7]
        
        #declare the symbols
        m,w=smp.symbols('m w', negative=False)
        
        ## Write the equations
        f = alpha_m*(M0 - m) - beta_m*(m-Mbar)*w**4 + sigma_m*w**2
        g = -alpha_w*w + beta_w*(m-Mbar)*w**3 - gamma_w*w**3 - sigma_m*m*w
        
        # Compute the solutions
        sol=smp.nonlinsolve([f,g], [m,w])
        
        ## Record the derivatives at all the solutions
        for k in range(0,len(sol)):
            
            ## sort the real solutions
            if complex(sol.args[k][0]).imag==0 and complex(sol.args[k][0]).real >= 0 and complex(sol.args[k][1]).imag==0 and complex(sol.args[k][1]).real >= 0:
                
                ## Read the solution
                ms=complex(sol.args[k][0]).real
                ws=complex(sol.args[k][1]).real
                
                ## computet the distance to laminar
                rlam=np.sqrt((ms-1)**2+(ws-0)**2)
                
                ## Compute the Jacobian
                fm = -alpha_m - beta_m*ws**4
                fw = -4*beta_m*(ms-Mbar)*ws**3 + 2*sigma_m*ws
                gm = beta_w*ws**3 - sigma_m*ws
                gw = -alpha_w + 3*beta_w*(ms-Mbar)*ws**2 - 3*gamma_w*ws**2 - sigma_m*ms
                
                ## Jacobian trace and determinant
                trJ=fm+gw
                detJ=(fm*gw)-(fw*gm)
                
                ## Isolate the turbulent fixed point    
                if trJ<0 and detJ>0:
                    
                    if rlam > 0:
                        
                        J[0,0,i]=fm
                        J[0,1,i]=fw
                        J[1,0,i]=gm
                        J[1,1,i]=gw
                
                
    ## Compute th derivative at Rc
    dJdR=(J[:,:,0]-J[:,:,1])/(2*dr)
    
    return dJdR
        
    
        
### Check the weakly nonlinear analysis
def weakly_nonlinear(self,rc):
    
    ## Parameters
    fp=np.asarray([self.cf1.value(),self.cf2.value(),self.cf3.value(),
                  self.cf4.value(),self.cf5.value(),self.cf6.value()])
    
    dm=self.cdm.value()
    dw=self.cdw.value()
    
    r=rc
    
    ## Get the model parameters
    alpha_m=prm(r,fp)[0]
    beta_m=prm(r,fp)[1]
    sigma_m=prm(r,fp)[2]
    alpha_w=prm(r,fp)[3]
    beta_w=prm(r,fp)[4]
    gamma_w=prm(r,fp)[5]
    Mbar=prm(r,fp)[6]
    M0=prm(r,fp)[7]
    
    #declare the symbols
    m,w=smp.symbols('m w', negative=False)
    
    ## Write the equations
    f = alpha_m*(M0 - m) - beta_m*(m-Mbar)*w**4 + sigma_m*w**2
    g = -alpha_w*w + beta_w*(m-Mbar)*w**3 - gamma_w*w**3 - sigma_m*m*w
    
    # Compute the solutions
    sol=smp.nonlinsolve([f,g], [m,w])
    
    ## Record the derivatives at all the solutions
    for k in range(0,len(sol)):
        
        ## sort the real solutions
        if complex(sol.args[k][0]).imag==0 and complex(sol.args[k][0]).real >= 0 and complex(sol.args[k][1]).imag==0 and complex(sol.args[k][1]).real >= 0:
            
            ## Read the solution
            ms=complex(sol.args[k][0]).real
            ws=complex(sol.args[k][1]).real
            
            ## computet the distance to laminar
            rlam=np.sqrt((ms-1)**2+(ws-0)**2)
            
            ## Compute the Jacobian
            fm = -alpha_m - beta_m*ws**4
            fw = -4*beta_m*(ms-Mbar)*ws**3 + 2*sigma_m*ws
            gm = beta_w*ws**3 - sigma_m*ws
            gw = -alpha_w + 3*beta_w*(ms-Mbar)*ws**2 - 3*gamma_w*ws**2 - sigma_m*ms
            
            ## Jacobian trace and determinant
            trJ=fm+gw
            detJ=(fm*gw)-(fw*gm)
            
            ## Isolate the turbulent fixed point    
            if trJ<0 and detJ>0:
                
                if rlam > 0:
                    
                    ## Define the derivatives of f(M,W)
                    ## 1st order derivatives
                    fm = -alpha_m - beta_m*ws**4
                    fw = -4*beta_m*(ms-Mbar)*ws**3 + 2*sigma_m*ws
                    
                    ## 2nd order
                    fmm=0
                    fmw=-4*beta_m*ws**3
                    fww=-12*beta_m*(ms-Mbar)*ws**2 + 2*sigma_m
                    
                    ## 3rd order
                    fmmm=0
                    fmmw=0
                    fmww=-12*beta_m*ws**2
                    fwww=-24*beta_m*(ms-Mbar)*ws
                    
                    ###### Define the derivatives of g(M,W)
                    ## 1st order
                    gm = beta_w*ws**3 - sigma_m*ws
                    gw = -alpha_w + 3*beta_w*(ms-Mbar)*ws**2 - 3*gamma_w*ws**2 - sigma_m*ms
                    
                    ## 2nd order
                    gmm=0
                    gmw=3*beta_w*ws**2-sigma_m
                    gww=6*beta_w*(ms-Mbar)*ws - 6*gamma_w*ws
                    
                    ## 3rd order
                    gmmm=0
                    gmmw=0
                    gmww=6*beta_w*ws
                    gwww=6*beta_w*(ms-Mbar)-6*gamma_w
                            
                    ## Critical wave number
                    kc=np.sqrt((dm*gw+dw*fm)/(2*dm*dw))
                    
                    ## step 1 : alpha and beta
                    alpha = - (fm-kc**2*dm) / fw
                    #alpha = - gm / (gw - kc**2*dw)
                    beta =  - (fm-kc**2*dm) / gm
                    #beta = - fw / (gw - kc**2*dw)
                    
                    ## step 2 : |phi|
                    detphi = (fm-4*kc**2*dm)*(gw-4*kc**2*dw) - fw*gm
                    
                    ## setp 3 : G1h, G2h
                    G1h = fmm + 2*alpha*fmw + alpha**2*fww
                    G2h = gmm + 2*alpha*gmw + alpha**2*gww
                    
                    ## Setp 4 : G1, G2
                    G1 = 0.5*fmm + alpha*fmw + 0.5*alpha**2*fww
                    G2 = 0.5*gmm + alpha*gmw + 0.5*alpha**2*gww
                    
                    ## step 5 : a0h, b0h
                    a0h = -1/(fm*gw-gm*fw) * (G1h*gw - G2h*fw)
                    b0h = -1/(fm*gw-gm*fw) * (G2h*fm - G1h*gm)
                    
                    ## step 6 : a2h, b2h
                    a2h = -1/detphi * (G1*(gw-4*kc**2*dw) - G2*fw)
                    b2h = -1/detphi * (G2*(fm-4*kc**2*dm) - G1*gm)
                    
                    ## step 7 : g1, g2
                    g1a = (a0h+a2h)*( fmm + alpha*fmw + beta*gmm + alpha*beta*gmw)
                    g1b = (b0h+b2h)*( fmw + alpha*fww + beta*gmw + alpha*beta*gww)
                    g1 = g1a+g1b
                    
                    g2 = fmmm + beta*gmmm + 3*alpha*(fmmw + beta*gmmw) + 3*alpha**2*(fmww + beta*gmww) + alpha**3*(fwww + beta*gwww)
                    
                    ## Final step 
                    g=1/(1+alpha*beta) * (g1 + 0.5*g2)
                    
                    
                    ### Compute the first coefficient eta
                    J1=(-rc)*DJ(self,rc)
                    eta = (1/(1+alpha*beta)) * ( J1[0,0] + alpha*J1[0,1] + beta*J1[1,0] + alpha*beta*J1[1,1])
                    
                    
                
                    
                    
    
    return eta,g,kc             
                    
                    
                    
                

    
    
        
    
    
    
            
        
        
        
        