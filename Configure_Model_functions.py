import numpy as np
import Phase_space_sim_functions as phsim
import Phase_plane_analysis_functions as ppf


def prm(self,R) :
    
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

    return np.asarray([alpha_m,beta_m,sigma_m,alpha_w,beta_w,gamma_w,Mbar,M0])

def plot_preview(self,M,W):
    
    r=self.R_sim_E.value()
    npts=self.Np_pp_E.value()
    
    #Plot the null clines
    m_wnull, wsam_w, m_mnull, wsam_m = ppf.plot_parameters(self,r,npts)[5:]
    self.preview_window.prp_phaxs.plot(m_wnull,wsam_w,color='r',label=r'$W_n$')
    self.preview_window.prp_phaxs.plot(m_mnull,wsam_m,color='b',label=r'$M_n$')
    self.preview_window.prp_phaxs.set_xlabel("M")
    self.preview_window.prp_phaxs.set_ylabel("W")
    self.preview_window.prp_phaxs.set_xlim([0,1])
    self.preview_window.prp_phaxs.set_ylim([0,0.1])
    
    ## Geometry
    x=np.linspace(0,self.Lx_sim_E.value(),self.Nx_sim_E.value())
    
    ## Plot
    ## Plot M(x)
    self.preview_window.prp_axs[0,0].plot(x,M,color='blue')
    self.preview_window.prp_axs[0,0].set_xlabel("x")
    self.preview_window.prp_axs[0,0].set_ylabel("M")
    
    ## Plot W(x)
    self.preview_window.prp_axs[1,0].plot(x,W,color='orange')
    self.preview_window.prp_axs[1,0].set_xlabel("x")
    self.preview_window.prp_axs[1,0].set_ylabel("W")
    
    ## Plot the phase space
    self.preview_window.prp_phaxs.plot(M,W,color='k')
    self.preview_window.prp_fig.tight_layout()
    
    self.preview_window.prp_canvas.draw()
    
    
    

def noise_plot(self):
    
    ## Noise location
    txt=self.ns_loc.currentText()

    if txt=="Fixed Point":
        
        ## Compute the fixed points
        if self.R_sim_E.value()>self.Rsn:
            m0,w0=phsim.compute_fixed_points(self,float(self.R_sim_E.value()))
        else:
            m0=1.0
            w0=0.0
            
    if txt=="Custom":
        
        m0,w0=self.ns_M.value(),self.ns_W.value()

    ## Generate gaussian noise
    M=np.random.normal(loc=m0,scale=self.ns_Mamp.value(),size=self.Nx_sim_E.value())
    W=np.random.normal(loc=w0,scale=self.ns_Wamp.value(),size=self.Nx_sim_E.value())

    ## Make it periodic
    M[-1]=M[0]
    W[-1]=W[0]
    
    return M,W

    
    
    
def pulse_plot(self):
    
    ## Pulse location in phase space
    txt=self.pulse_loc.currentText()
    
    if txt=="Fixed Point":
        
        ## Compute the fixed points
        if self.R_sim_E.value()>self.Rsn:
            m0,w0=phsim.compute_fixed_points(self,float(self.R_sim_E.value()))
        else:
            m0=1.0
            w0=0.0
            
    if txt=="Custom":
        
        m0,w0=self.pulse_M.value(),self.pulse_W.value()
        
    ## Geometry
    x=np.linspace(0,self.Lx_sim_E.value(),self.Nx_sim_E.value())
        
    ### Add / Create a single pulse   
    if self.npulse.value()==1:
                    
        ## Pulse location in space
        xtxt=self.pulse_xloc.currentText()
        
        if xtxt=="Midpoint":
            x0=self.Lx_sim_E.value()/2.0
        if xtxt=="Custom":
            x0=self.pulse_x.value()
    
        ## Generate gaussian pulse
        M=m0 + (self.pulse_Mamp.value()-m0)*np.exp(-(x-x0)**2 / self.pulse_Mstd.value())
        W=w0 + (self.pulse_Wamp.value()-w0)*np.exp(-(x-x0)**2 / self.pulse_Wstd.value())
        
        
        
    ### Add / Create multiple pulses
    if self.npulse.value()>1:
        
        ## equidistant locations of pulses
        x0=np.linspace(0,self.Lx_sim_E.value(),self.npulse.value()+2)[1:-1]
        
        M=np.zeros(len(x))
        W=np.zeros(len(x))
        ## Generate gaussian pulses
        for i in range(0,len(x0)):
            
            M+=m0 + (self.pulse_Mamp.value()-m0)*np.exp(-(x-x0[i])**2 / self.pulse_Mstd.value())
            W+=w0 + (self.pulse_Wamp.value()-w0)*np.exp(-(x-x0[i])**2 / self.pulse_Wstd.value())
            
        ## Rescale M and W
        M-=np.ones(len(x))*m0*(self.npulse.value()-1)
        W-=np.ones(len(x))*w0*(self.npulse.value()-1)
        
    ## Make it periodic
    M[-1]=M[0]
    W[-1]=W[0]    
    
    return M,W

def pattern_plot(self):
        
    ## Pulse location in phase space
    txt=self.pat_loc.currentText()
    
    if txt=="Fixed Point":
        
        ## Compute the fixed points
        if self.R_sim_E.value()>self.Rsn:
            m0,w0=phsim.compute_fixed_points(self,float(self.R_sim_E.value()))
        else:
            m0=1.0
            w0=0.0
            
    if txt=="Custom":
        
        m0,w0=self.pat_M.value(),self.pat_W.value()
        
    
    ## Geometry
    x=np.linspace(0,self.Lx_sim_E.value(),self.Nx_sim_E.value())
    
    ## Generate sin and cosine signals
    k=2*np.pi*self.patnpks.value()/self.Lx_sim_E.value()
    M=m0 + self.pat_Mamp.value()*np.sin(k*x)
    W=w0 + self.pat_Wamp.value()*np.cos(k*x)
    
    ## Make it periodic
    M[-1]=M[0]
    W[-1]=W[0]
    
    return M,W


### Chimera Initial condition
def chimera_plot(self):
    
    ## initialize
    M=np.zeros(self.Nx_sim_E.value())
    W=np.zeros(self.Nx_sim_E.value())
    
    ## Geometry
    x=np.linspace(0,self.Lx_sim_E.value(),self.Nx_sim_E.value())
    
    ## Get the proportionalities
    prp=np.zeros(self.ch_IC_list.count()+1)
    for i in range(0,len(prp)-1):
        prp[i+1]=prp[i]+self.chimera_IC[i]["Proportion"]
        
    
    ## compile the IC
    for i in range(0,self.ch_IC_list.count()):
                
        x1=x[(x>=prp[i]*self.Lx_sim_E.value()) & (x<prp[i+1]*self.Lx_sim_E.value())]
        idx1=np.argmin(np.abs(x1[0]-x))
        idx2=np.argmin(np.abs(x1[-1]-x))
        
        xp=x1-x1[0]
        
        
        ## Add a pulse
        if self.chimera_IC[i]["type"] == "Pulse" :
            
            ## Pulse location in phase space
            txt=self.chimera_IC[i]["Location"]
            
            if txt=="Fixed Point":
                
                ## Compute the fixed points
                if self.R_sim_E.value()>self.Rsn:
                    m0,w0=phsim.compute_fixed_points(self,float(self.R_sim_E.value()))
                else:
                    m0=1.0
                    w0=0.0
                    
            if txt=="Custom":
                
                m0,w0=self.chimera_IC[i]["M"],self.chimera_IC[i]["W"]
                    
            ### Add / Create a single pulse   
            if self.chimera_IC[i]["No of pulses"]==1:
                            
                ## Pulse location in space
                xtxt=self.chimera_IC[i]["x location"]
                
                if xtxt=="Midpoint":
                    x0=xp.max()/2.0
                if xtxt=="Custom":
                    x0=self.chimera_IC[i]["x"]
            
                ## Generate gaussian pulse
                M[idx1:idx2+1]=m0 + (self.chimera_IC[i]["Mamp"]-m0)*np.exp(-(xp-x0)**2 / self.chimera_IC[i]["Mstd"])
                W[idx1:idx2+1]=w0 + (self.chimera_IC[i]["Wamp"]-w0)*np.exp(-(xp-x0)**2 / self.chimera_IC[i]["Wstd"])
                
                
                
            ### Add / Create multiple pulses
            if self.chimera_IC[i]["No of pulses"]>1:
                
                ## equidistant locations of pulses
                x0=np.linspace(0,xp.max(),self.chimera_IC[i]["No of pulses"]+2)[1:-1]
                
                ## Generate gaussian pulses
                for k in range(0,len(x0)):
                    
                    M[idx1:idx2+1]+=m0 + (self.chimera_IC[i]["Mamp"]-m0)*np.exp(-(xp-x0[k])**2 / self.chimera_IC[i]["Mstd"])
                    W[idx1:idx2+1]+=w0 + (self.chimera_IC[i]["Wamp"]-w0)*np.exp(-(xp-x0[k])**2 / self.chimera_IC[i]["Wstd"])
                
                ## Rescale M and W
                M[idx1:idx2+1]-=np.ones(len(xp))*m0*(self.chimera_IC[i]["No of pulses"]-1)
                W[idx1:idx2+1]-=np.ones(len(xp))*w0*(self.chimera_IC[i]["No of pulses"]-1)
                
                
        ## Add a pattern
        if self.chimera_IC[i]["type"] == "Pattern" :
            
            ## Pulse location in phase space
            txt=self.chimera_IC[i]["Location"]
            
            if txt=="Fixed Point":
                
                ## Compute the fixed points
                if self.R_sim_E.value()>self.Rsn:
                    m0,w0=phsim.compute_fixed_points(self,float(self.R_sim_E.value()))
                else:
                    m0=1.0
                    w0=0.0
                    
            if txt=="Custom":
                
                m0,w0=self.chimera_IC[i]["M"],self.chimera_IC[i]["W"]
                
        
            ## Generate sin and cosine signals
            k=2*np.pi*self.chimera_IC[i]["Npks"]/xp.max()
            M[idx1:idx2+1]=m0 + self.pat_Mamp.value()*np.sin(k*xp)
            W[idx1:idx2+1]=w0 + self.pat_Wamp.value()*np.sin(k*xp)
            
            
            
        ## Add Noise
        if self.chimera_IC[i]["type"] == "Noise" :
            
            txt=self.chimera_IC[i]["Location"]
        
            if txt=="Fixed Point":
                
                ## Compute the fixed points
                if self.R_sim_E.value()>self.Rsn:
                    m0,w0=phsim.compute_fixed_points(self,float(self.R_sim_E.value()))
                else:
                    m0=1.0
                    w0=0.0
                    
            if txt=="Custom":
                
                m0,w0=self.chimera_IC[i]["M"],self.chimera_IC[i]["W"]
        
            ## Generate gaussian noise
            M[idx1:idx2+1]=np.random.normal(loc=m0,scale=self.chimera_IC[i]["Mamp"],size=len(xp))
            W[idx1:idx2+1]=np.random.normal(loc=w0,scale=self.chimera_IC[i]["Wamp"],size=len(xp))
                                
                        
    ## Make the signal periodic
    M[-1]=M[0]
    W[-1]=W[0]    
    
    return M,W                
                    
    
        
        
        
        
    
    

    