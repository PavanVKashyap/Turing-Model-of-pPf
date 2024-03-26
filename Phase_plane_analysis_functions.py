import numpy as np
import Phase_space_sim_functions as phsim
import Configure_Model_functions as cnfg

### Get plot parameters
def plot_parameters(self,r,npts):
    
    ## get the parameters
    alpha_m=cnfg.prm(self,r)[0]
    beta_m=cnfg.prm(self,r)[1]
    sigma_m=cnfg.prm(self,r)[2]
    alpha_w=cnfg.prm(self,r)[3]
    beta_w=cnfg.prm(self,r)[4]
    gamma_w=cnfg.prm(self,r)[5]
    Mbar=cnfg.prm(self,r)[6]
    M0=cnfg.prm(self,r)[7]
    
    ## Sample space for nullclines
    wsam_max=1
    nsam=10**5
    wsam_m=np.linspace(0,wsam_max,nsam)
    wsam_w=np.linspace(1.01*np.sqrt(sigma_m/beta_w),wsam_max,nsam)
    
    ## Nullclines coordinates
    m_mnull=( (sigma_m*wsam_m**2) + (beta_m * Mbar * wsam_m**4) + (alpha_m * M0) ) / (alpha_m + (beta_m*wsam_m**4) ) 
    m_wnull=( alpha_w + (gamma_w*wsam_w**2) + (beta_w*Mbar*wsam_w**2)  ) / ( (beta_w*wsam_w**2) - sigma_m) 
    
    ## Sampling for quivers
    msam=np.linspace(0.00001,self.PP_Mmax.value(),npts)
    wsam=np.linspace(0.00001,self.PP_Wmax.value(),npts)
    
    mv=np.reshape(msam,[npts,1])*np.ones([npts,npts])
    wv=np.reshape(wsam,[1,npts])*np.ones([npts,npts])
    
    dmdt=alpha_m*(M0 - mv) - beta_m*(mv-Mbar)*wv**4 + sigma_m*wv**2
    dwdt=-alpha_w*wv + beta_w*(mv-Mbar)*wv**3 - gamma_w*wv**3 - sigma_m*mv*wv
    
    ## Normalize quivers
    vel=np.sqrt(dmdt**2+dwdt**2)

    
    return msam, wsam, dmdt, dwdt, vel, m_wnull, wsam_w, m_mnull, wsam_m


### Phase plot
def Plot_phase_plots(self):
    
    r=float(self.R_pp_sl.value())    
    npts=int(self.Np_pp_E.value())

    msam, wsam, dmdt, dwdt, vel = plot_parameters(self,r,75)[:5]
    
    ## Clear the plots
    for i in range(0,3):
        if self.cb_present[i]:
            self.cb[i].remove()
            self.cb_present[i]=False
        self.PP_axs[i].cla()
    
    ## Plot the dMdt phase space
    ctemp=np.min([np.abs(dmdt.min()),np.abs(dmdt.max())])
    cmin=-ctemp/3
    cmax=ctemp/3
    pt1=self.PP_axs[0].pcolormesh(msam,wsam,dmdt.T,cmap='seismic',shading='gouraud',vmin=cmin,vmax=cmax)
    self.PP_axs[0].set_xlabel("M")
    self.PP_axs[0].set_ylabel("W")
    self.PP_axs[0].set_title(r"$dMdt$")
    self.cb[0]=self.PP_fig.colorbar(pt1,ax=self.PP_axs[0])
    self.cb_present[0]=True

    
    ## Plot the dWdt phase space
    ctemp=np.min([np.abs(dwdt.min()),np.abs(dwdt.max())])
    cmin=-ctemp/1.5
    cmax=ctemp/1.5
    pt2=self.PP_axs[1].pcolormesh(msam,wsam,dwdt.T,cmap='seismic',shading='gouraud',vmin=cmin,vmax=cmax)
    self.PP_axs[1].set_xlabel("M")
    self.PP_axs[1].set_ylabel("W")
    self.PP_axs[1].set_title(r"$dWdt$")
    self.cb[1]=self.PP_fig.colorbar(pt2,ax=self.PP_axs[1])
    self.cb_present[1]=True
    
    #Plot the null clines
    m_wnull, wsam_w, m_mnull, wsam_m = plot_parameters(self,r,npts)[5:]
    self.PP_axs[2].plot(m_wnull,wsam_w,color='r',label=r'$W_n$')
    self.PP_axs[2].plot(m_mnull,wsam_m,color='b',label=r'$M_n$')
    self.PP_axs[2].set_title(f"Phase plot for R={r}")
    self.PP_axs[2].set_xlabel("M")
    self.PP_axs[2].set_ylabel("W")
    self.PP_axs[2].set_xlim([0,self.PP_Mmax.value()])
    self.PP_axs[2].set_ylim([0,self.PP_Wmax.value()])
    self.PP_axs[2].legend(loc=1)
    
    if self.Qv_pp_btn.isChecked():
            
        ## Plot the quivers
        msam, wsam, dmdt, dwdt, vel = plot_parameters(self,r,npts)[:5]
        pt3=self.PP_axs[2].quiver(msam,wsam,(dmdt/vel).T,(dwdt/vel).T,vel,cmap='jet')
        self.cb[2]=self.PP_fig.colorbar(pt3,ax=self.PP_axs[2])
        self.cb_present[2]=True
        
    self.PP_fig.tight_layout()
    self.PP_canvas.draw()
    
    
def quiver_change(self):
    
    ## Claer the plot
    self.PP_axs[2].cla()
    if self.cb_present[2]:
        self.cb[2].remove()
        self.cb_present[2]=False
    
    ## Get the inputs
    r=float(self.R_pp_sl.value())    
    npts=int(self.Np_pp_E.value())

    
    #Plot the null clines
    m_wnull, wsam_w, m_mnull, wsam_m = plot_parameters(self,r,npts)[5:]
    self.PP_axs[2].plot(m_wnull,wsam_w,color='r',label=r'$W_n$')
    self.PP_axs[2].plot(m_mnull,wsam_m,color='b',label=r'$M_n$')
    self.PP_axs[2].set_title(f"Phase plot for R={r}")
    self.PP_axs[2].set_xlabel("M")
    self.PP_axs[2].set_ylabel("W")
    self.PP_axs[2].set_xlim([0,self.PP_Mmax.value()])
    self.PP_axs[2].set_ylim([0,self.PP_Wmax.value()])
    self.PP_axs[2].legend(loc=1)
    
    if self.Qv_pp_btn.isChecked():
            
        ## Plot the quivers
        msam, wsam, dmdt, dwdt, vel = plot_parameters(self,r,npts)[:5]
        pt3=self.PP_axs[2].quiver(msam,wsam,(dmdt/vel).T,(dwdt/vel).T,vel,cmap='jet')
        self.cb[2]=self.PP_fig.colorbar(pt3,ax=self.PP_axs[2])
        self.cb_present[2]=True
        
    self.PP_fig.tight_layout()
    self.PP_canvas.draw()






    

    