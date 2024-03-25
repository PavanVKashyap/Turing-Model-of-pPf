from PyQt5.uic import loadUiType
from PyQt5 import QtGui
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.animation import FuncAnimation
import glob
import os


import Phase_plane_analysis_functions as ppf
import Phase_space_sim_functions as phsim
import Configure_Model_functions as cnfg
import simulation as sims
import Stability_calculation_functions as stb
import Analyze_Model as anmd

Ui_MainWindow, QMainWindow = loadUiType('Full_App.ui')
preview_MainWindow, QPreviewMainWindow = loadUiType('Preview_plot.ui') 
spt_MainWindow, QsptMainWindow = loadUiType('SPT.ui') 

## Create the class for the preview plot window
class Preplot_Window(preview_MainWindow, QPreviewMainWindow):
    def __init__(self,):
        super(Preplot_Window, self).__init__()
        self.setupUi(self)
        
        ## Preview plots   
        self.prp_fig=Figure(figsize=(15,5))
        self.prp_axs=self.prp_fig.subplots(nrows=2,ncols=2)
        self.prp_gs1=self.prp_axs[0,1].get_gridspec()
        self.prp_axs[0,1].remove()
        self.prp_axs[1,1].remove()
        self.prp_phaxs=self.prp_fig.add_subplot(self.prp_gs1[:,1])
        self.prp_canvas = FigureCanvas(self.prp_fig)
        self.preview_plot.addWidget(self.prp_canvas)
        self.prp_fig.tight_layout()
        self.prp_canvas.draw()
        
        ## add the toolbar
        prp_toolbar = NavigationToolbar(self.prp_canvas, self)
        self.preview_plot.addWidget(prp_toolbar)

## Create the class for the Space time plot window
class SPT_Window(spt_MainWindow, QsptMainWindow):
    def __init__(self,):
        super(SPT_Window, self).__init__()
        self.setupUi(self)
        
        ## Space time plots   
        self.sptfig=Figure(figsize=(8,15))
        self.sptaxs=self.sptfig.subplots(nrows=1,ncols=2)
        self.spt_canvas = FigureCanvas(self.sptfig)
        self.spt.addWidget(self.spt_canvas)
        self.sptfig.tight_layout()
        self.spt_canvas.draw()
        
        ## add the toolbar
        spt_toolbar = NavigationToolbar(self.spt_canvas, self)
        self.spt.addWidget(spt_toolbar)

##### Create the main class for the window
class Main(QMainWindow, Ui_MainWindow):
    def __init__(self, ):
        super(Main, self).__init__()
        self.setupUi(self)


###############--------------- Configure Model Tab ------------------------------############
        
        ## Plot Phase portrait   
        self.anfig=Figure(figsize=(15,5))
        self.anxs=self.anfig.subplots(nrows=2,ncols=3)
        self.angs1=self.anxs[0,0].get_gridspec()
        self.anxs[0,0].remove()
        self.anxs[1,0].remove()
        self.anphxs=self.anfig.add_subplot(self.angs1[:,0])
        self.ancanvas = FigureCanvas(self.anfig)
        self.condPlot.addWidget(self.ancanvas)
        self.anfig.tight_layout()
        self.ancanvas.draw()
        
        ## add the toolbar
        antoolbar = NavigationToolbar(self.ancanvas, self)
        self.condPlot.addWidget(antoolbar)
        
        ## Nullcline plot
        anmd.plot_nullcline(self)
        self.cR.editingFinished.connect(self.anplot_nullcline)
        self.cf1.editingFinished.connect(self.anplot_nullcline)
        self.cf2.editingFinished.connect(self.anplot_nullcline)
        self.cf3.editingFinished.connect(self.anplot_nullcline)
        self.cf4.editingFinished.connect(self.anplot_nullcline)
        self.cf5.editingFinished.connect(self.anplot_nullcline)
        self.cf6.editingFinished.connect(self.anplot_nullcline)
        
        self.cMmin.editingFinished.connect(self.anplot_nullcline)
        self.cMmax.editingFinished.connect(self.anplot_nullcline)
        self.cWmin.editingFinished.connect(self.anplot_nullcline)
        self.cWmax.editingFinished.connect(self.anplot_nullcline)
        
        ## Check the stability criteria
        self.cAnalyze.clicked.connect(self.analyze_bifurcation)
        

###############--------------- Phase Portrait Tab -------------------------------############

        ## Plot Phase portrait        
        self.PP_fig=Figure(figsize=(15,5))
        self.PP_axs=self.PP_fig.subplots(nrows=1,ncols=3)
        self.PP_canvas = FigureCanvas(self.PP_fig)
        self.PPplot.addWidget(self.PP_canvas)
        self.PP_fig.tight_layout()
        self.PP_canvas.draw()
        
        ## add the toolbar
        PP_toolbar = NavigationToolbar(self.PP_canvas, self)
        self.PPplot.addWidget(PP_toolbar)
        
        ## Initial plot
        self.cb=[[],[],[]]
        self.cb_present=[False,False,False]
        ppf.Plot_phase_plots(self)
        
        ## Link slider for Phase Portrait
        self.R_pp_sl.valueChanged.connect(self.change_R_pp)
        self.R_pp_sl.setMinimum(200)
        self.R_pp_sl.setMaximum(500)
        self.R_pp_sl.setValue(200)
        
        
        ## Remove / Plot quiver if button unselected
        self.Qv_pp_btn.toggled.connect(self.change_quiver_pp)
        
        ## Change quiver plot
        self.Np_pp_E.valueChanged.connect(self.change_Np)
               
        
#############---------------- Stability Analysis Tab -----------------------------#############

        ## Stability analysis plots   
        self.st_fig=Figure(figsize=(15,5))
        self.st_axs=self.st_fig.subplots(nrows=1,ncols=2)
        self.st_canvas = FigureCanvas(self.st_fig)
        self.stability_plots.addWidget(self.st_canvas)
        self.st_fig.tight_layout()
        self.st_canvas.draw()
        
        ## add the toolbar
        st_toolbar = NavigationToolbar(self.st_canvas, self)
        self.stability_plots.addWidget(st_toolbar)
        
        ### default value of Rsn
        self.Rsn=311
        
        ## Set the function for change in R value
        self.stR.valueChanged.connect(self.change_stRdisp)
        
        ## Change the Slider values
        self.stR1.valueChanged.connect(self.change_slider_range)
        self.stR2.valueChanged.connect(self.change_slider_range)
        
        ## Find the fixed points
        self.stFix.clicked.connect(self.find_fixed_points)
        
        ## Find Rc
        self.stRc.clicked.connect(self.find_Rc)
        
        
##############--------------- Phase Space Analysis Tab ---------------------------#############
        

        ## Phase Space Simulation plots   
        self.Ph_fig=Figure(figsize=(15,5))
        self.Ph_axs=self.Ph_fig.subplots(nrows=2,ncols=2)
        self.Ph_gs1=self.Ph_axs[0,1].get_gridspec()
        self.Ph_axs[0,1].remove()
        self.Ph_axs[1,1].remove()
        self.Ph_phaxs=self.Ph_fig.add_subplot(self.Ph_gs1[:,1])
        self.Ph_canvas = FigureCanvas(self.Ph_fig)
        self.Phplot.addWidget(self.Ph_canvas)
        self.Ph_fig.tight_layout()
        self.Ph_canvas.draw()
        
        ## add the toolbar
        Ph_toolbar = NavigationToolbar(self.Ph_canvas, self)
        self.Phplot.addWidget(Ph_toolbar)
        
        ## Draw phase plot
        self.ph_cb=[]
        self.ph_cb_present=False
        
        ## Start simulations
        self.start_ph_btn.clicked.connect(self.Phsim_start)
        self.stop_ph_btn.clicked.connect(self.Phsim_stop)
        self.pause_ph_btn.clicked.connect(self.Phsim_pause)
        
        ## Plot the full trajectory
        self.traj_ph_btn.clicked.connect(self.Phsim_pltTraj)
     
        ## Create a solution vector for the Phase plane analysis
        self.Mph=[]
        self.Wph=[]
        self.Tph=[]
        self.start_Phsim_ani=None
        self.paused_count=0
        self.start_count=0
        
        
############---------------------- Simulation Tab --------------------------------###############

        ## Simulation plots   
        self.sim_fig=Figure(figsize=(15,5))
        self.sim_axs=self.sim_fig.subplots(nrows=2,ncols=2)
        self.sim_gs1=self.sim_axs[0,1].get_gridspec()
        self.sim_axs[0,1].remove()
        self.sim_axs[1,1].remove()
        self.sim_phaxs=self.sim_fig.add_subplot(self.sim_gs1[:,1])
        self.sim_canvas = FigureCanvas(self.sim_fig)
        self.sim_plot.addWidget(self.sim_canvas)
        self.sim_fig.tight_layout()
        self.sim_canvas.draw()
        
        ## add the toolbar
        sim_toolbar = NavigationToolbar(self.sim_canvas, self)
        self.sim_plot.addWidget(sim_toolbar)
        
        # Start simulations
        self.start_sim_btn.clicked.connect(self.start_simulation)
        
        # Pause simulations
        self.pause_sim_btn.clicked.connect(self.pause_simulation)
        self.sim_paused_count=0
        
        # Stop simulations
        self.stop_sim_btn.clicked.connect(self.stop_simulation)
        
        ## Add button
        self.add_sim_btn.clicked.connect(self.modify_state)
        
        ## Initialize the space time plot window
        self.spt_window=SPT_Window()
        ## View SPT button
        self.sim_spt.clicked.connect(self.view_spt)
        
        ## Declare the solution variables
        self.Msim=None
        self.Wsim=None
        
        ## Declare teh animation variable
        self.sim_animation=None
        

#############-------------------- Model Parameters Tab ----------------------------###############
        
        ## Initialize the preview window
        self.preview_window=Preplot_Window()
        
        ## Noise preview button
        self.nspr_btn.clicked.connect(self.noise_show_preview)
        
        ## Pulse preview button
        self.pulsepr_btn.clicked.connect(self.pulse_show_preview)
        
        ## Pattern preview button
        self.patpr_btn.clicked.connect(self.pattern_show_preview)
        
        ## Set defaults
        self.Model_default.clicked.connect(self.MD_set_default)
        self.NS_default.clicked.connect(self.ns_set_default)
        self.pulse_default.clicked.connect(self.pulse_set_default)
        self.pat_default.clicked.connect(self.pat_set_default)
        
        ### Chimera IC 
        self.chimera_IC=[]
        self.chimera_prp=1
        self.ch_add.clicked.connect(self.chICadd)
        self.ch_remove_all.clicked.connect(self.chICremoveAll)
        self.ch_remove.clicked.connect(self.chRemove)
        self.ch_up.clicked.connect(self.MoveUp)
        self.ch_down.clicked.connect(self.MoveDown)
        self.ch_preview.clicked.connect(self.chimera_preview)
        
        self.ch_IC_list.itemSelectionChanged.connect(self.move_button_control)
        
           
###############################################################################################
###################                                                     #######################
###################             Functions for interaction               #######################
###################                                                     #######################
###############################################################################################
    

#############---------------- confingure Model Tab -------------------------------############

    
    ## Plot the null clines
    def anplot_nullcline(self):
        
        ## Plot the nullcline
        anmd.plot_nullcline(self)
        
    ## Analyze the bifurcations
    def analyze_bifurcation(self):
        
        ## Analyze and plot the criteria
        anmd.bifurcation_criteria(self)


##############---------------- Phase Portrait Tab --------------------------------############
    def change_R_pp(self):
        
        ## Plot 
        ppf.Plot_phase_plots(self)
        
        ## Display the R value
        self.R_pp_disp.display(self.R_pp_sl.value())
        
    def change_quiver_pp(self):
        
        ## Change quiver in Phase portrait
        ppf.quiver_change(self)
        
    def change_Np(self):
        
        ## change the quiver plot
        ppf.quiver_change(self)
        
##############------------------- Stability Analysis Tab -------------------------###########
    def change_slider_range(self):
        
        ## Set the slider range and position in the phase portrait tab
        self.R_pp_sl.setMinimum(int(self.stR1.value()))
        self.R_pp_sl.setMaximum(int(self.stR2.value()))
        self.R_pp_sl.setValue(int(self.stR1.value()))

        
    def change_stRdisp(self):
        
        ## Plot 
        if self.stR.value()>self.Rsn:
            self.st_axs[1].cla()
            stb.stability_criteria_plot(self)
                
        
    def find_fixed_points(self):
        
        ## Disable all entry options
        self.stR1.setEnabled(False)
        self.stR2.setEnabled(False)
        self.stk1.setEnabled(False)
        self.stk2.setEnabled(False)
        self.stR.setEnabled(False)
        self.stFix.setEnabled(False)
        self.stRc.setEnabled(False)
        
        ## Compute fixed points and plot
        self.st_axs[0].cla()
        stb.compute_fixed_points(self)
        
        ## Enable all entry options
        self.stR1.setEnabled(True)
        self.stR2.setEnabled(True)
        self.stk1.setEnabled(True)
        self.stk2.setEnabled(True)
        self.stR.setEnabled(True)
        self.stFix.setEnabled(True)
        self.stRc.setEnabled(True)
        
    def find_Rc(self):
        
        ## Disable all entry options
        self.stR1.setEnabled(False)
        self.stR2.setEnabled(False)
        self.stk1.setEnabled(False)
        self.stk2.setEnabled(False)
        self.stR.setEnabled(False)
        self.stFix.setEnabled(False)
        self.stRc.setEnabled(False)
        
        ## compute the turing bifurcation value
        self.st_axs[1].cla()
        stb.find_turing(self)
        
        ## Enable all entry options
        self.stR1.setEnabled(True)
        self.stR2.setEnabled(True)
        self.stk1.setEnabled(True)
        self.stk2.setEnabled(True)
        self.stR.setEnabled(True)
        self.stFix.setEnabled(True)
        self.stRc.setEnabled(True)
        
    
                
                
        
        

        
##############------------------- Phase plane analysis Tab -----------------------##########

    def Plot_phase_plots_Phtab(self):
        
        r=float(self.R_ph_E.value())    
        npts=int(self.Np_pp_E.value())
        
        ## Plot M(t)
        self.Ph_axs[0,0].plot(0,self.M0_ph_E.value(),animated=True,color='b',label="M(t)")
        self.Ph_axs[0,0].plot(np.arange(0,1000),np.ones(1000)*self.ms,animated=True,color='k',label="Fixed Point")
        self.Ph_axs[0,0].set_xlabel("Time")
        self.Ph_axs[0,0].set_ylabel("M(t)")
        # self.Ph_axs[0,0].legend(loc='best')
        self.Ph_axs[0,0].set_xlim([0,1000])
        self.Ph_axs[0,0].set_ylim([0,self.PP_Mmax.value()])
        
        ## Plot W(t)
        self.Ph_axs[1,0].plot(0,self.W0_ph_E.value(),animated=True,color='orange',label="W(t)")
        self.Ph_axs[1,0].plot(np.arange(0,1000),np.ones(1000)*self.ws,animated=True,color='k',label="Fixed Point")
        self.Ph_axs[1,0].set_xlabel("Time")
        self.Ph_axs[1,0].set_ylabel("W(t)")
        # self.Ph_axs[1,0].legend(loc='best')
        self.Ph_axs[1,0].set_xlim([0,1000])
        self.Ph_axs[1,0].set_ylim([0,self.PP_Wmax.value()])
        
        ## Plot phase plot
        m_wnull, wsam_w, m_mnull, wsam_m = ppf.plot_parameters(self,r,npts)[5:]
        self.Ph_phaxs.plot(m_wnull,wsam_w,color='r',animated=True,label=r'$W_n$')
        self.Ph_phaxs.plot(m_mnull,wsam_m,color='b',animated=True,label=r'$M_n$')
        self.Ph_scat=self.Ph_phaxs.scatter(self.M0_ph_E.value(),self.W0_ph_E.value(),c=[0],s=[5],cmap='jet',marker='o',zorder=0)
        self.Ph_phaxs.set_xlabel("M")
        self.Ph_phaxs.set_ylabel("W")
        self.Ph_phaxs.set_xlim([0,self.PP_Mmax.value()])
        self.Ph_phaxs.set_ylim([0,self.PP_Wmax.value()])
        # self.Ph_phaxs.legend(loc='best')
        
        self.Ph_fig.tight_layout()
        self.Ph_canvas.draw()
        
        return (self.Ph_axs[0,0].lines[0],self.Ph_axs[0,0].lines[1],
                self.Ph_axs[1,0].lines[0],self.Ph_axs[1,0].lines[1],
                self.Ph_phaxs.lines[0], self.Ph_phaxs.lines[1],
                self.Ph_scat)

    def update_Ph_sim_plot(self,i):
    
        ## Call noise free Runge Kutta or With Noise Euler Maruyama
        if self.NS_ph_check.isChecked():
            phsim.Euler_Maruyama(self,i)
        else:
            phsim.Runge_Kutta(self,i)
            
        ## Compute the time vector
        self.Tph.append(i*self.dt_ph_E.value())
        
        ## Compute the fixed points
        if self.R_ph_E.value()>self.Rsn:
            self.ms,self.ws=phsim.compute_fixed_points(self,float(self.R_ph_E.value()))
        else:
            self.ms=1.0
            self.ws=0.0
    

        ## M(t)
        self.Ph_axs[0,0].lines[0].set_data(np.asarray(self.Tph),np.asarray(self.Mph))
        self.Ph_axs[0,0].lines[1].set_ydata(np.ones(1000)*self.ms)
        
        ## W(t)
        self.Ph_axs[1,0].lines[0].set_data(np.asarray(self.Tph),np.asarray(self.Wph))
        self.Ph_axs[1,0].lines[1].set_ydata(np.ones(1000)*self.ws)
        
        ## Update trajectory in phase space
        m_wnull, wsam_w, m_mnull, wsam_m = ppf.plot_parameters(self,self.R_ph_E.value(),self.Np_pp_E.value())[5:]
        self.Ph_phaxs.lines[0].set_data(m_wnull,wsam_w)
        self.Ph_phaxs.lines[1].set_data(m_mnull,wsam_m)
        
        temp=np.vstack((self.Mph,self.Wph)).T
        self.Ph_scat.set_offsets(temp)
        
        carray=np.arange(len(self.Tph))
        self.Ph_scat.set_array(np.asarray(carray))
        self.Ph_scat.set_clim([carray[0],carray[-1]])
        
        sz=np.ones(len(self.Tph))*5
        sz[-1]=120
        self.Ph_scat.set_sizes(sz)
        
        return (self.Ph_axs[0,0].lines[0],self.Ph_axs[0,0].lines[1],
                self.Ph_axs[1,0].lines[0],self.Ph_axs[1,0].lines[1],
                self.Ph_phaxs.lines[0], self.Ph_phaxs.lines[1],
                self.Ph_scat)
        
        
    def Phsim_start(self):
        
        self.start_count+=1
        
        if self.start_count%2!=0:
        
            if self.start_count>1:
                self.start_Phsim_ani._stop()
                
            self.start_ph_btn.setText("Re-start")
                
        
        if self.start_count%2==0:
            self.start_Phsim_ani._stop()
            self.start_ph_btn.setText("Start")
            
        ## Clear All plots
        self.Ph_axs[0,0].cla()
        self.Ph_axs[1,0].cla()
        self.Ph_phaxs.cla()
        
        ##Reset the paused state
        self.paused_count=0
        self.pause_ph_btn.setText("Pause")
            
            
        ## Set the initial condition
        self.Mph=[self.M0_ph_E.value()]
        self.Wph=[self.W0_ph_E.value()]
        self.Tph=[self.dt_ph_E.value()]

        ## Compute the fixed points
        if self.R_ph_E.value()>self.Rsn:
            self.ms,self.ws=phsim.compute_fixed_points(self,float(self.R_ph_E.value()))
        else:
            self.ms=1.0
            self.ws=0.0
        
        ## Start animation
        self.start_Phsim_ani = FuncAnimation(self.Ph_canvas.figure, self.update_Ph_sim_plot, init_func=self.Plot_phase_plots_Phtab, blit=True, interval=10)
        
    def Phsim_pause(self):
        
        self.paused_count+=1
        
        if self.paused_count%2!=0:
            self.start_Phsim_ani.pause()
            self.pause_ph_btn.setText("Resume")
            
        if self.paused_count%2==0:
            self.start_Phsim_ani.resume()
            self.pause_ph_btn.setText("Pause")
            
    def Phsim_stop(self):
        
        if self.start_count>0:
            self.start_Phsim_ani._stop()
            self.start_count=0
            self.paused_count=0
            self.pause_ph_btn.setText("Pause")
            self.start_ph_btn.setText("Start")
            
    def Phsim_pltTraj(self):
        
        
        ## Stop the animation if in progress
        if self.start_count>0:
            self.start_Phsim_ani._stop()
            self.start_count=0
            self.paused_count=0
            self.pause_ph_btn.setText("Pause")
            self.start_ph_btn.setText("Start")
        
        ## Clear the axes
        self.Ph_axs[0,0].cla()
        self.Ph_axs[1,0].cla()
        self.Ph_phaxs.cla()
        
        ## Disable all inputs and buttons
        self.R_ph_E.setEnabled(False)
        self.dt_ph_E.setEnabled(False)
        self.NS_ph_E.setEnabled(False)
        self.M0_ph_E.setEnabled(False)
        self.W0_ph_E.setEnabled(False)
        self.iter_ph_E.setEnabled(False)
        
        self.start_ph_btn.setEnabled(False)
        self.stop_ph_btn.setEnabled(False)
        self.pause_ph_btn.setEnabled(False)
        
        ## Call the function to compute and plot the trajectory
        if self.NS_ph_check.isChecked():
            phsim.Euler_Maruyama_full_traj(self)
        else:
            phsim.Runge_Kutta_full_traj(self)
            
        ## Enable all inputs and buttons
        self.R_ph_E.setEnabled(True)
        self.dt_ph_E.setEnabled(True)
        self.NS_ph_E.setEnabled(True)
        self.M0_ph_E.setEnabled(True)
        self.W0_ph_E.setEnabled(True)
        self.iter_ph_E.setEnabled(True)
        
        self.start_ph_btn.setEnabled(True)
        self.stop_ph_btn.setEnabled(True)
        self.pause_ph_btn.setEnabled(True)
        
        
        
############ ----------------------------- Simulation Tab ----------------------------###########
    

    def sim_init_plot(self):
        
        r=float(self.R_sim_E.value())    
        npts=int(self.Np_pp_E.value())
        ## Geometry
        x=np.linspace(0,self.Lx_sim_E.value(),self.Nx_sim_E.value())

        ## Plot M(x)
        self.sim_axs[0,0].plot(x,self.Msim,animated=True,color='b')
        self.sim_axs[0,0].set_xlabel("x")
        self.sim_axs[0,0].set_ylabel("M(x)")
        self.sim_axs[0,0].set_ylim([0,self.PP_Mmax.value()])
        
        ## Plot W(x)
        self.sim_axs[1,0].plot(x,self.Wsim,animated=True,color='orange')
        self.sim_axs[1,0].set_xlabel("x")
        self.sim_axs[1,0].set_ylabel("W(x)")
        self.sim_axs[1,0].set_ylim([0,self.PP_Wmax.value()])

        ## Plot phase plot
        m_wnull, wsam_w, m_mnull, wsam_m = ppf.plot_parameters(self,r,npts)[5:]
        self.sim_phaxs.plot(m_wnull,wsam_w,color='r',animated=False,label=r'$W_n$')
        self.sim_phaxs.plot(m_mnull,wsam_m,color='b',animated=False,label=r'$M_n$')
        self.sim_phaxs.plot(self.Msim,self.Wsim,color='k',animated=True)
        self.sim_phaxs.plot(self.Msim[0],self.Wsim[0],color='yellow',ls='None',marker='o',ms=8,animated=True)
        self.sim_phaxs.set_xlabel("M")
        self.sim_phaxs.set_ylabel("W")
        self.sim_phaxs.set_xlim([0,self.PP_Mmax.value()])
        self.sim_phaxs.set_ylim([0,self.PP_Wmax.value()])
        
        
        ## Add the time of simulation
        self.simtime= self.sim_phaxs.text(0.65, 0.9, '', transform=self.sim_phaxs.transAxes)
        
        self.sim_fig.tight_layout()
        self.sim_canvas.draw()
        
        return (self.sim_axs[0,0].lines[0],
                self.sim_axs[1,0].lines[0],
                self.sim_phaxs.lines[2],self.sim_phaxs.lines[3],self.simtime)

    def update_sim_plot(self,i):
        
        sims.simulation_time_step(self)
        
        self.sim_axs[0,0].lines[0].set_ydata(self.Msim)
        self.sim_axs[1,0].lines[0].set_ydata(self.Wsim)
        self.sim_phaxs.lines[2].set_data(self.Msim,self.Wsim)
        self.sim_phaxs.lines[3].set_data(self.Msim[0],self.Wsim[0])
        
        ##Update time
        self.simtime.set_text(f"Time={self.solver.sim_time:0.2f}")
        
        
        return (self.sim_axs[0,0].lines[0],
                self.sim_axs[1,0].lines[0],
                self.sim_phaxs.lines[2],self.sim_phaxs.lines[3],self.simtime)

    def start_simulation(self):
        
        ## Clear all the plots
        self.sim_axs[0,0].cla()
        self.sim_axs[1,0].cla()
        self.sim_phaxs.cla()
        
        ## Remove the old data
        tcheck=glob.glob("case")
        if len(tcheck)>=1:
            os.system("rm -r case")
        
        ## Disable the editing options
        self.R_sim_E.setEnabled(False)
        self.Lx_sim_E.setEnabled(False)
        self.Nx_sim_E.setEnabled(False)
        self.NS_sim_E.setEnabled(False)
        self.NS_sim_check.setEnabled(False)
        self.adv_sim_check.setEnabled(False)
        self.start_sim_btn.setEnabled(False)
        self.spt_dt.setEnabled(False)
        
        ## Initialize the solver
        self.solver=sims.solver_initialize(self)
        
        ## Start animation
        self.sim_animation = FuncAnimation(self.sim_canvas.figure, self.update_sim_plot, init_func=self.sim_init_plot, blit=True, interval=10)
        
        
        
    def pause_simulation(self):
        
        self.sim_paused_count+=1
        
        if self.sim_paused_count%2!=0:
            self.sim_animation.pause()
            self.pause_sim_btn.setText("Resume")
            
        if self.sim_paused_count%2==0:
            self.sim_animation.resume()
            self.pause_sim_btn.setText("Pause")
            
    def stop_simulation(self):
        
        ## Start animation
        self.sim_animation._stop()
        
        ## Enable the editing options
        self.R_sim_E.setEnabled(True)
        self.Lx_sim_E.setEnabled(True)
        self.Nx_sim_E.setEnabled(True)
        self.NS_sim_E.setEnabled(True)
        self.NS_sim_check.setEnabled(True)
        self.adv_sim_check.setEnabled(True)
        self.start_sim_btn.setEnabled(True)
        self.spt_dt.setEnabled(True)
        
    def modify_state(self):
        
        #Condition
        txt=self.IC_select.currentText()
        
        if txt=='Noise':
            
            Mch,Wch=cnfg.noise_plot(self)
            
        if txt=='Pulse':
            
            Mch,Wch=cnfg.pulse_plot(self)
            
        if txt=='Pattern':
            
            Mch,Wch=cnfg.pattern_plot(self)
            
        if txt=='Chimera':
            
            Mch,Wch=cnfg.chimera_plot(self)
            
        if txt=="Self":
            
            Mch=np.roll(self.Msim,np.random.randint(0,self.Nx_sim_E.value()))
            Wch=np.roll(self.Wsim,np.random.randint(0,self.Nx_sim_E.value()))
        
            
        ## Modify the state
        self.m.set_scales(1)
        self.w.set_scales(1)
        self.m['g']+=Mch
        self.w['g']+=Wch
        
        
    ### View Space time
    def view_spt(self):
        
        tcheck=glob.glob("case")
        if len(tcheck)>=1:
            sims.spt_plot(self)
            self.spt_window.show()
        else:
            print("No data available")
            self.spt_window.show()
        
        
            
############ ----------------------------- Model Parameters Tab -----------------------###########

    def noise_show_preview(self):
        
        self.preview_window.prp_axs[0,0].cla()
        self.preview_window.prp_axs[1,0].cla()
        self.preview_window.prp_phaxs.cla()
        # Get the data
        M,W=cnfg.noise_plot(self)
        ## Plot
        cnfg.plot_preview(self, M, W)
        self.preview_window.show()
    
    def pulse_show_preview(self):
        
        self.preview_window.prp_axs[0,0].cla()
        self.preview_window.prp_axs[1,0].cla()
        self.preview_window.prp_phaxs.cla()
        
        M,W=cnfg.pulse_plot(self)
        cnfg.plot_preview(self, M, W)
        
        self.preview_window.show()
        
    def pattern_show_preview(self):
        
        self.preview_window.prp_axs[0,0].cla()
        self.preview_window.prp_axs[1,0].cla()
        self.preview_window.prp_phaxs.cla()
        
        M,W=cnfg.pattern_plot(self)
        cnfg.plot_preview(self, M, W)
        
        self.preview_window.show()
        
    def MD_set_default(self):
        
        self.f1.setValue(5.0)
        self.f2.setValue(5.0)
        self.f3.setValue(0.2)
        self.f4.setValue(0.2)
        self.f5.setValue(0.2)
        self.f6.setValue(5.0)
        self.dm.setValue(10.0)
        self.dw.setValue(0.2)
        
    def ns_set_default(self):
        
        self.ns_M.setValue(1.0)
        self.ns_W.setValue(0.0)
        self.ns_Mamp.setValue(0.04)
        self.ns_Wamp.setValue(0.01)
        
    def pulse_set_default(self):
        
        self.pulse_M.setValue(1.0)
        self.pulse_W.setValue(0.0)
        self.pulse_x.setValue(25)
        self.pulse_Mamp.setValue(0.5)
        self.pulse_Mstd.setValue(1.5)
        self.pulse_Wamp.setValue(0.08)
        self.pulse_Wstd.setValue(1.5)
        
    def pat_set_default(self):
        
        self.patnpks.setValue(10)
        self.pat_M.setValue(1.0)
        self.pat_W.setValue(0.0)
        self.pat_Mamp.setValue(0.05)
        self.pat_Wamp.setValue(0.01)
        
    #### chimera IC
    def chICadd(self):
        
        
        ## Add the item to list
        if self.ch_IC_select.currentText() == "Pulse" :
            
            pulse={"type" : "Pulse",
                   "Proportion" : self.ch_prp.value(),
                   "No of pulses" : self.npulse.value(),
                   "Location" : self.pulse_loc.currentText(),
                   "M" : self.pulse_M.value(),
                   "W" : self.pulse_W.value(),
                   "x location" : self.pulse_xloc.currentText(),
                   "x" : self.pulse_x.value(),
                   "Mamp" : self.pulse_Mamp.value(),
                   "Mstd" : self.pulse_Mstd.value(),
                   "Wamp" : self.pulse_Wamp.value(),
                   "Wstd" : self.pulse_Wstd.value()}
            
            self.ch_IC_list.addItem(f"Pulse-{self.ch_prp.value()}")
            self.chimera_IC.append(pulse)
            self.chimera_prp-=self.ch_prp.value()
            self.ch_prp.setMaximum(self.chimera_prp)
            
        if self.ch_IC_select.currentText() == "Pattern" :
            
            pattern={"type":"Pattern",
                     "Proportion" : self.ch_prp.value(),
                     "Npks": self.patnpks.value(),
                     "Location":self.pat_loc.currentText(),
                     "M" : self.pat_M.value(),
                     "W" : self.pat_W.value(),
                     "Mamp" : self.pat_Mamp.value(),
                     "Wamp" : self.pat_Wamp.value()}
            
            self.ch_IC_list.addItem(f"Pattern-{self.ch_prp.value()}")
            self.chimera_IC.append(pattern)
            self.chimera_prp-=self.ch_prp.value()
            self.ch_prp.setMaximum(self.chimera_prp)
            
        if self.ch_IC_select.currentText() == "Noise" :
            
            noise={"type" : "Noise",
                   "Proportion" : self.ch_prp.value(),
                   "Location" : self.ns_loc.currentText(),
                   "M" : self.ns_M.value(),
                   "W" : self.ns_W.value(),
                   "Mamp" : self.ns_Mamp.value(),
                   "Wamp" : self.ns_Wamp.value()}
            
            self.ch_IC_list.addItem(f"Noise-{self.ch_prp.value()}")
            self.chimera_IC.append(noise)
            self.chimera_prp-=self.ch_prp.value()
            self.ch_prp.setMaximum(self.chimera_prp)
            
        
        ## Restrict adding more IC
        if self.chimera_prp==0:
            
            self.ch_add.setEnabled(False)
            
        if self.chimera_prp>0:
            
            self.ch_add.setEnabled(True)
            
    ## Remove single item from list
    def chRemove(self):
        
        ## Remove selected item
        row=self.ch_IC_list.currentRow()
        self.ch_IC_list.takeItem(row)
        
        ## Update the proportionality
        self.chimera_prp+=self.chimera_IC[row]["Proportion"]
        self.ch_prp.setMaximum(self.chimera_prp)
        
        ## Remove the item from list
        self.chimera_IC.remove(self.chimera_IC[row])
        
        ## Reset the add button
        if self.chimera_prp>0:
            
            self.ch_add.setEnabled(True)
        
        
    ## Clear the entire list
    def chICremoveAll(self):
        
        # Clear the list
        self.chimera_IC=[]
        
        ## Clear the list widget
        self.ch_IC_list.clear()
        
        ## Reset the maxima of proportionality
        self.chimera_prp=1
        self.ch_prp.setMaximum(self.chimera_prp)
        
        ## Reset the Add button
        self.ch_add.setEnabled(True)
        
    
    ## Move the items up
    def MoveUp(self):
        
        ## Get the row index
        rowindex=self.ch_IC_list.currentRow()
        
        ## Move the IC in the list
        temp1=self.chimera_IC[rowindex]
        temp2=self.chimera_IC[rowindex-1]
        
        self.chimera_IC[rowindex-1]=temp1
        self.chimera_IC[rowindex]=temp2
        
        ## Get the item from the list
        currentitem=self.ch_IC_list.takeItem(rowindex)
        
        ## Insert it above
        self.ch_IC_list.insertItem(rowindex-1,currentitem)
        
        ## Change the selection
        self.ch_IC_list.setCurrentRow(rowindex-1)
        
    ## Move the items Down
    def MoveDown(self):
        
        ## Get the row index
        rowindex=self.ch_IC_list.currentRow()
        
        ## Move the IC in the list
        temp1=self.chimera_IC[rowindex]
        temp2=self.chimera_IC[rowindex+1]
        
        self.chimera_IC[rowindex+1]=temp1
        self.chimera_IC[rowindex]=temp2
        
        ## Get the item from the list
        currentitem=self.ch_IC_list.takeItem(rowindex)
        
        ## Insert it below
        self.ch_IC_list.insertItem(rowindex+1,currentitem)
        
        ## Change the selection
        self.ch_IC_list.setCurrentRow(rowindex+1)        
        
        
    ## Restrict the move buttons
    def move_button_control(self):
        
        ## Condition for up button
        self.ch_up.setEnabled(not bool(self.ch_IC_list.currentRow()==0))
        
        ## Condition for down button
        self.ch_down.setEnabled(not bool(self.ch_IC_list.currentRow()==self.ch_IC_list.count()-1))
        
    ## Get the preview
    def chimera_preview(self):
        
        self.preview_window.prp_axs[0,0].cla()
        self.preview_window.prp_axs[1,0].cla()
        self.preview_window.prp_phaxs.cla()
        
        M,W=cnfg.chimera_plot(self)
        cnfg.plot_preview(self, M, W)
        
        self.preview_window.show()
        
    
        
        

if __name__ == '__main__':
    import sys
    from PyQt5 import QtWidgets
    
        
    app = QtWidgets.QApplication(sys.argv)
    main = Main()
    main.show()
    sys.exit(app.exec_())

        