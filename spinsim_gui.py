#Spinsim_GSL GUI
#GUI interface for Spinsim_GSL command-line program
#Written I. Yulaev (iyulaev@ucsd.edu) 2010-02-21
#University of California, San Diego ECE Department
#Magnet icon taken from openoffice.org website

import time
import re
import thread

#Stuff for PyQt
import sys
from PyQt4 import QtCore, QtGui, Qt
from qt_ui_v2 import Ui_Dialog
from TimeDep_Handler import TimeDep_Handler

#stuff for command line calls
from subprocess import *
import os

#stuff for creating the plots
from numpy import *
from pylab import *
import matplotlib.pyplot as plt

#TODO: Add to tree diagram: handleMakePlotT, external process call for that, handleSetLimits, button for that

print_command_only = 0

#Parameters for switching SI -> CGS (the names are defined backwards, sorry)
AM_PER_OE = 1.0/79.53
K_PER_ERGCC = 10.0
MS_PER_EMUCC = 1.0/1000.0

class MyForm(QtGui.QMainWindow):

	cgs_units = 0

	"""This class initializes the GUI and also contains code for reading in the parameters and running the spinsim C code"""
	def __init__(self, parent=None):
		""" This function initializes the GUI and sets up all of the default parameters """
		QtGui.QWidget.__init__(self, parent)
		self.ui = Ui_Dialog()
		self.ui.setupUi(self)
		
		#make sure functions know we're in the init method right now
		self.initializing = 1
		
		#set up all of the variables (with some default values) 
		#that we will eventually extract from the gui
		self.sv_thetaoff = 0.05
		self.sv_percK2 = 0.5
		self.sv_Alpha = [0.01,0.01,0.01,0.01]
		self.sv_K1 = [3e5, 3e6, 0.0, 0.0]
		self.sv_Ms = [650000, 600000, 0, 0]
		self.sv_lposz = [0.0, 7.0, 0.0, 0.0]
		self.sv_simtype = 0
		self.sv_jmin = -1e11
		self.sv_jmax = 1e11
		self.sv_n = 30
		self.sv_hzmin = -1e6
		self.sv_hzmax = 1e6
		self.sv_nfield = 30
		self.sv_tf = 2.5e-8
		self.sv_j_timdep_funs = 1
		self.sv_hz_timdep_funs = 1
		self.sv_dx = 50
		self.sv_dy = 50
		self.sv_temperature = 0
		self.sv_initmag_x = [0.0,0.0,0.0,0.0]
		self.sv_initmag_y = [0.1, 0.0, 0.0, 0.0]
		self.sv_initmag_z = [0.9, 1, 0, 0]
		self.sv_dz = [3,3,0,0]
		self.sv_pol = [0.3,0.3,0.3,0.3]
		self.sv_layerEn = [1,1,0,0]
		self.sv_layerFixed = [0,1,0,0]
		self.sv_numthreads = 1
		self.sv_numFrames = 50
		self.sv_fieldInPlane = 0
		
		#relaxation time to not apply field nor current
		self.sv_relax = 0.0
		
		#af coupling terms: coupling b/w 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
		self.sv_af = [0.0,0.0,0.0,0.0,0.0,0.0]
		
		#constant variable; list of currently implemented simtype's
		self.VALID_SIMTYPES = [6, 9, 10, 26, 41, 6+(1<<31), 9+(1<<31), 10+(1<<31), 26+(1<<31), 41+(1<<31)]
		
		#handle switching to CGS units
		self.switch_units_to_cgs(self.cgs_units)
		#[/end if cgs_units]

		#set up all of the GUI fields to the "default" values
		#for non-layer dependent things
		self.variablesToDialogNonLayer()
		
		#set layer 1 as the first layer that's selected
		self.handleLayer1Button()
		
		#setup progress bar
		self.ui.progressBar_Progress.setVisible(False)
		self.ui.label_Progress.setVisible(False)
		
		#setup threading lock
		self.spinsim_thread_lock = thread.allocate_lock()
		self.progress_var_lock = thread.allocate_lock()
		self.progress_var = 0
		
		#setup threading lock for plot limit poll-er
		self.spinsim_plot_lock = thread.allocate_lock()
		self.plot_lims = [0,0,0,0]
		self.plot_lims_updated = 0
		
		#setup timer for updating progress bar when the spinsim c program is runnings
		self.progress_timer = QtCore.QTimer()
		QtCore.QObject.connect(self.progress_timer, QtCore.SIGNAL("timeout()"), self.handleUpdateProgress)
		
		#setup the comboBox (and its label) for plotting
		self.ui.comboBox_PlotType.setEditable(False)
		self.ui.comboBox_PlotType.setVisible(True)
		#self.ui.label_GeneratePlot.setVisible(False)
		#self.ui.pushButton_MakePlot.setVisible(False)
		
		#make sure number of frames is not visible (until we select the plot type such that it should be)
		self.ui.label_numFrames.setVisible(False)
		self.ui.lineEdit_numFrames.setVisible(False)
		
		#hide the layer 4 button for now
		self.ui.pushButton_layer4.setVisible(False)
		
		#make sure functions know we're no longer in the init method
		self.initializing = 0
		
	#[/__init__()]
	
	def switch_units_to_cgs(self, cgs_units):
		"""This function translates the SI units entered into the GUI into CGS, and also re-sets the GUI
		labels."""
		if(cgs_units != 0):
			count = 0
			
			#Switch all labels
			self.ui.label_5.setText("Saturation Magnetization, in emu/cc")
			self.ui.label_18.setText("Layer Anisotropy, erg/cc")
			self.ui.label_11.setText("Min/Max Field & # of Points (Oe)")
			
			#Switch entered values
			for element in self.sv_K1:
				self.sv_K1[count] = element * K_PER_ERGCC
				count = count + 1
				
			count = 0
				
			for element in self.sv_Ms:
				self.sv_Ms[count] = element * MS_PER_EMUCC
				count = count+1
			
			self.sv_hzmin = self.sv_hzmin * AM_PER_OE 
			self.sv_hzmax = self.sv_hzmax * AM_PER_OE
			
	#[/switch_units_to_cgs()]
	
	def switch_units_to_si(self, cgs_units):
		"""This function translates the SI units entered into the GUI into CGS, and also re-sets the GUI
		labels."""
		if(cgs_units != 0):
			count = 0
			
			#Switch entered values
			for element in self.sv_K1:
				self.sv_K1[count] = element / K_PER_ERGCC
				count = count + 1
				
			count = 0
				
			for element in self.sv_Ms:
				self.sv_Ms[count] = element / MS_PER_EMUCC
				count = count+1
			
			self.sv_hzmin = self.sv_hzmin / AM_PER_OE 
			self.sv_hzmax = self.sv_hzmax / AM_PER_OE
			
	#[/switch_units_to_si()]
	
	def printDebugVars(self):
		"""This function prints out all of the currently stored internal variables"""
		
		self.dialogToVariablesNonLayer()
		
		print("sv_thetaoff = %f" % self.sv_thetaoff)
		print("sv_Alpha = %f" % self.sv_Alpha)
		print("sv_K1 = [%f, %f, %f, %f]" % (self.sv_K1[0],self.sv_K1[1],self.sv_K1[2],self.sv_K1[3]))
		print("sv_Ms = [%f, %f, %f, %f]" % (self.sv_Ms[0],self.sv_Ms[1],self.sv_Ms[2],self.sv_Ms[3]))
		print("sv_lposz = [%f, %f, %f, %f" % (self.sv_lposz[0],self.sv_lposz[1],self.sv_lposz[2],self.sv_lposz[3]))
		print("sv_jmin/sv_jmax/sv_n = %f / %f / %f" % (self.sv_jmin, self.sv_jmax, self.sv_n))
		
	def dialogToVariablesNonLayer(self):
		"""This function grabs all of the text in the non-layer-dependent fields
		and puts it into the GUIs internal variable storage"""
		#self.sv_Alpha = float(self.ui.lineEdit_Alpha.text())
		self.sv_jmin = float(self.ui.lineEdit_jmin.text())
		self.sv_jmax = float(self.ui.lineEdit_jmax.text())
		self.sv_n = int(self.ui.lineEdit_n.text())
		self.sv_hzmin = float(self.ui.lineEdit_hzmin.text())
		self.sv_hzmax = float(self.ui.lineEdit_hzmax.text())
		self.sv_nfield = int(self.ui.lineEdit_nfield.text())
		self.sv_tf = float(self.ui.lineEdit_IntegrationTime.text())
		self.sv_temperature = float(self.ui.lineEdit_SampTemp.text())
		self.sv_percK2 = float(self.ui.lineEdit_perK2.text())
		self.sv_thetaoff = float(self.ui.lineEdit_thetaoff.text())
		self.sv_numthreads = int(self.ui.lineEdit_threadcount.text())
		self.sv_dx = float(self.ui.lineEdit_dx.text())
		self.sv_dy = float(self.ui.lineEdit_dy.text())
		
		if(self.ui.checkBox_fieldInPlane.checkState() == QtCore.Qt.Checked):
			self.sv_fieldInPlane = 1
		else:
			self.sv_fieldInPlane = 0
		
		#af coupling terms
		self.sv_af[0] = float(self.ui.lineEdit_af_12.text());
		self.sv_af[1] = float(self.ui.lineEdit_af_13.text());
		self.sv_af[2] = float(self.ui.lineEdit_af_14.text());
		self.sv_af[3] = float(self.ui.lineEdit_af_23.text());
		self.sv_af[4] = float(self.ui.lineEdit_af_24.text());
		self.sv_af[5] = float(self.ui.lineEdit_af_34.text());
		
		#relaxation time
		self.sv_relax = float(self.ui.lineEdit_relax.text());
	#[/dialogToVariablesNonLayer]
	
	def variablesToDialogNonLayer(self):
		"""This function takes all of the non-layer-dependent variables from memory and writes
		them into the dialog text boxes"""
		#self.ui.lineEdit_Alpha.setText(str(self.sv_Alpha))
		self.ui.lineEdit_jmin.setText(str(self.sv_jmin))
		self.ui.lineEdit_jmax.setText(str(self.sv_jmax))
		self.ui.lineEdit_n.setText(str(self.sv_n))
		self.ui.lineEdit_hzmin.setText(str(self.sv_hzmin))
		self.ui.lineEdit_hzmax.setText(str(self.sv_hzmax))
		self.ui.lineEdit_nfield.setText(str(self.sv_nfield))
		self.ui.lineEdit_IntegrationTime.setText(str(self.sv_tf))
		self.ui.lineEdit_SampTemp.setText(str(self.sv_temperature))
		self.ui.lineEdit_perK2.setText(str(self.sv_percK2))
		self.ui.lineEdit_thetaoff.setText(str(self.sv_thetaoff))
		self.ui.lineEdit_threadcount.setText(str(self.sv_numthreads))
		self.ui.lineEdit_dx.setText(str(self.sv_dx))
		self.ui.lineEdit_dy.setText(str(self.sv_dy))
		self.ui.lineEdit_numFrames.setText(str(self.sv_numFrames))
				
		if(self.sv_fieldInPlane == 1):
			self.ui.checkBox_layerFixed.setCheckState(QtCore.Qt.Checked)
		else:
			self.ui.checkBox_layerFixed.setCheckState(QtCore.Qt.Unchecked)
		
		#af coupling terms
		self.ui.lineEdit_af_12.setText(str(self.sv_af[0]));
		self.ui.lineEdit_af_13.setText(str(self.sv_af[1]));
		self.ui.lineEdit_af_14.setText(str(self.sv_af[2]));
		self.ui.lineEdit_af_23.setText(str(self.sv_af[3]));
		self.ui.lineEdit_af_24.setText(str(self.sv_af[4]));
		self.ui.lineEdit_af_34.setText(str(self.sv_af[5]));
		
		#relaxation time
		self.ui.lineEdit_relax.setText(str(self.sv_relax));
		
	#[/variablesToDialogNonLayer()]
		
	def handleLayer1Button(self):
		"""Handler for layer 1 button - stores the old layer's fields to internal variables,
		sets up the push button colors and GUI labels, and writes from internal variables
		to the new layer's GUI fields"""
		#store the old variable's fields, if the program is NOT initializing
		if(self.initializing != 1):
			self.storeLayern(self.ui.curr_layer)
		#set up the push buttons, lables, etc
		self.ui.pushButton_layer1.setStyleSheet("* { color: red }");
		self.ui.pushButton_layer2.setStyleSheet("* { color: black }");
		self.ui.pushButton_layer3.setStyleSheet("* { color: black }");
		self.ui.pushButton_layer4.setStyleSheet("* { color: black }");
		self.ui.label_layerParams.setText("Layer 1 Parameters");
		self.ui.curr_layer = 1
		#write from the internal variables to the GUI fields for the new layer
		self.handleLayernButton(self.ui.curr_layer)
		
	def handleLayer2Button(self):
		"""Handler for layer 1 button - stores the old layer's fields to internal variables,
		sets up the push button colors and GUI labels, and writes from internal variables
		to the new layer's GUI fields"""
		#store the old variable's fields, if the program is NOT initializing
		if(self.initializing != 1):
			self.storeLayern(self.ui.curr_layer)
		#set up the push buttons, lables, etc
		self.ui.pushButton_layer1.setStyleSheet("* { color: black }");
		self.ui.pushButton_layer2.setStyleSheet("* { color: red }");
		self.ui.pushButton_layer3.setStyleSheet("* { color: black }");
		self.ui.pushButton_layer4.setStyleSheet("* { color: black }");
		self.ui.label_layerParams.setText("Layer 2 Parameters");
		self.ui.curr_layer = 2
		#write from the internal variables to the GUI fields for the new layer
		self.handleLayernButton(self.ui.curr_layer)
		
	def handleLayer3Button(self):
		"""Handler for layer 3 button - stores the old layer's fields to internal variables,
		sets up the push button colors and GUI labels, and writes from internal variables
		to the new layer's GUI fields"""
		#store the old variable's fields, if the program is NOT initializing
		if(self.initializing != 1):
			self.storeLayern(self.ui.curr_layer)
		#set up the push buttons, lables, etc
		self.ui.pushButton_layer1.setStyleSheet("* { color: black }");
		self.ui.pushButton_layer2.setStyleSheet("* { color: black }");
		self.ui.pushButton_layer3.setStyleSheet("* { color: red }");
		self.ui.pushButton_layer4.setStyleSheet("* { color: black }");
		self.ui.label_layerParams.setText("Layer 3 Parameters");
		self.ui.curr_layer = 3
		#write from the internal variables to the GUI fields for the new layer
		self.handleLayernButton(self.ui.curr_layer)
		
	def handleLayer4Button(self):
		"""Handler for layer 4 button - stores the old layer's fields to internal variables,
		sets up the push button colors and GUI labels, and writes from internal variables
		to the new layer's GUI fields"""
		#store the old variable's fields, if the program is NOT initializing
		if(self.initializing != 1):
			self.storeLayern(self.ui.curr_layer)
		#set up the push buttons, lables, etc
		self.ui.pushButton_layer1.setStyleSheet("* { color: black }");
		self.ui.pushButton_layer2.setStyleSheet("* { color: black }");
		self.ui.pushButton_layer3.setStyleSheet("* { color: black }");
		self.ui.pushButton_layer4.setStyleSheet("* { color: red }");
		self.ui.label_layerParams.setText("Layer 4 Parameters");
		self.ui.curr_layer = 4
		#write from the internal variables to the GUI fields for the new layer
		self.handleLayernButton(self.ui.curr_layer)
		
	def storeLayern(self, n):
		"""This method stores the current GUI settings into the variables for layer n"""
		self.sv_Ms[n-1] = float(self.ui.lineEdit_SatMag.text())
		self.sv_K1[n-1] = float(self.ui.lineEdit_K1.text())
		self.sv_lposz[n-1] = float(self.ui.lineEdit_layerpos_z.text())
		self.sv_dz[n-1] = float(self.ui.lineEdit_dz.text())
		self.sv_initmag_x[n-1] = float(self.ui.lineEdit_InitMagX.text())
		self.sv_initmag_y[n-1] = float(self.ui.lineEdit_InitMagY.text())
		self.sv_initmag_z[n-1] = float(self.ui.lineEdit_InitMagZ.text())
		self.sv_Alpha[n-1] = float(self.ui.lineEdit_Alpha.text())
		self.sv_pol[n-1] = float(self.ui.lineEdit_Pol.text())
		
		if(self.ui.checkBox_layerEn.checkState() == QtCore.Qt.Checked):
			self.sv_layerEn[n-1] = 1
		else:
			self.sv_layerEn[n-1] = 0
			
		if(self.ui.checkBox_layerFixed.checkState() == QtCore.Qt.Checked):
			self.sv_layerFixed[n-1] = 1
		else:
			self.sv_layerFixed[n-1] = 0
	#[/storeLayern()]
	
	def handleLayernButton(self, n):
		"""This method sets up the GUI fields based on the current variables stored for
		layer n. Typically gets called when the <layer n> button gets pressed"""
		self.ui.lineEdit_SatMag.setText(str(self.sv_Ms[n-1]))
		self.ui.lineEdit_K1.setText(str(self.sv_K1[n-1]))
		self.ui.lineEdit_layerpos_z.setText(str(self.sv_lposz[n-1]))
		self.ui.lineEdit_dz.setText(str(self.sv_dz[n-1]))
		self.ui.lineEdit_InitMagX.setText(str(self.sv_initmag_x[n-1]))
		self.ui.lineEdit_InitMagY.setText(str(self.sv_initmag_y[n-1]))
		self.ui.lineEdit_InitMagZ.setText(str(self.sv_initmag_z[n-1]))
		self.ui.lineEdit_Alpha.setText(str(self.sv_Alpha[n-1]))
		self.ui.lineEdit_Pol.setText(str(self.sv_pol[n-1]))
		
		#set up checkboxes
		if(self.sv_layerEn[n-1] == 1):
			self.ui.checkBox_layerEn.setCheckState(QtCore.Qt.Checked)
		else:
			self.ui.checkBox_layerEn.setCheckState(QtCore.Qt.Unchecked)
		
		if(self.sv_layerFixed[n-1] == 1):
			self.ui.checkBox_layerFixed.setCheckState(QtCore.Qt.Checked)
		else:
			self.ui.checkBox_layerFixed.setCheckState(QtCore.Qt.Unchecked)
	#[/handleLayernButton()]
	
	def handleUpdateProgress(self):
		"""This function gets called by a timer to update the progress bar"""
		#print("Updating progress...")
		#acquire lock and set progress bar
		self.progress_var_lock.acquire()
		self.ui.progressBar_Progress.setValue(self.progress_var)
		
		#if we get to 100%, we're done, stop the timer
		if(self.progress_var == 100):
			self.progress_timer.stop()
			self.ui.comboBox_PlotType.setVisible(True)
			self.ui.label_GeneratePlot.setVisible(True)
			self.ui.pushButton_MakePlot.setVisible(True)
		
		self.progress_var_lock.release()
		
	def handleSetLimits(self):
		"""This function gets called when the user hits the "Limits from Plot" button in the GUI. It updates the
			j and Hz limits based on the last zoom-box limit settings from the last figure call. The settings themselves
			get passed to spinsim_gui through the stdout pipe from the spinsim C engine. """
			
		#check to see if our plot limits were updated. If so, post the updates to the GUI
		self.spinsim_plot_lock.acquire()
		
		if(self.plot_lims_updated == 1):
			self.ui.lineEdit_jmax.setText(str(self.plot_lims[3]));
			self.ui.lineEdit_hzmax.setText(str(self.plot_lims[1]));
			self.ui.lineEdit_jmin.setText(str(self.plot_lims[2]));
			self.ui.lineEdit_hzmin.setText(str(self.plot_lims[0]));
			self.plot_lims_updated = 0
					
		self.spinsim_plot_lock.release()
		
	def handleTimeDepDialog(self):
		"""This function gets called when the Setup Time Dependent Functions button is pressed. It
		initializes the time dependent function setup dialog and displays it. The dialog itself passes
		the new settings back to spinsim_gui via callbacks"""
		tdd = TimeDep_Handler(self, self.sv_hz_timdep_funs, self.sv_j_timdep_funs)
		
	def makePlotSlotT(self,plotToMake,numFrames):
		"""This function calls spinsim_make_plots.py via the shell to create the figure window and draw our plots. It also creates
			*.png file with the plots in the file"""
		
		#if we're doing an animation, we must also give spinsim_make_plots the # of frames we want
		if(plotToMake == 6):	
			cmdstr = "python spinsim_make_plots.py %d %d" % ((plotToMake+1), (numFrames+1))
		#default to '2' layers for plotting an x-y time domain plot, since as of
		#2011-09-25 I handle this well even if one of the layers isn't generated
		elif(plotToMake == 5):
			cmdstr = "python spinsim_make_plots.py %d %d" % ((plotToMake+1), 2)
		else:
			cmdstr = "python spinsim_make_plots.py %d" % (plotToMake+1)
			
		#create plots by calling a new process
		plot_process = Popen(cmdstr, shell=True, stdout=PIPE, cwd=os.getcwd())
		
		#loop to poll status of plot and update limits
		while(plot_process.poll() == None):
			time.sleep(1)
			
			#grab stdout
			plot_stdout = plot_process.stdout
			#Ivan: Why does this line only return when we exit the program?
			#because I wasn't flushing stdout in the other process...duh...
			line = plot_stdout.readline()
			#print("DEBUG: Got line \"%s\"" % line)
			line = line.rstrip('\n') #chomp
			
			#get the progress from stdout
			lims = re.match('\S+ = (\S+) (\S+) (\S+) (\S+)',line)
			if(lims != None):
				print("Spinsim_GUI: New limits are %f %f %f %f" % (float(lims.group(1)), float(lims.group(2)), float(lims.group(3)), float(lims.group(4))))
				#set the progress var
				self.spinsim_plot_lock.acquire()
				
				self.plot_lims[0] = float(lims.group(1))
				self.plot_lims[1] = float(lims.group(2))
				self.plot_lims[2] = float(lims.group(3))
				self.plot_lims[3] = float(lims.group(4))
				self.plot_lims_updated = 1
					
				self.spinsim_plot_lock.release()
			#[/if(lims != None)
		#[/while plot_process.poll() == None]
		
	def handlePlotTypeSelect(self, selected):
		"""This handler is called when a new plot type is selected in the plot type drop-down menu.
		If the user selects the Animation plot type, make visible the text entry field for number of
		frames"""
		if(selected == 6):
			self.ui.label_numFrames.setVisible(True)
			self.ui.lineEdit_numFrames.setVisible(True)
		else:
			self.ui.label_numFrames.setVisible(False)
			self.ui.lineEdit_numFrames.setVisible(False)
			
		#print("Selected plot type %d\n" % selected)
		
	def makePlotSlot(self):
		"""This function calls makePlotSlotT in a new thread. Used for calling
		the python program that actually creates the plots using matplotlib."""
		plotToMake = self.ui.comboBox_PlotType.currentIndex()
		
		#get the # of frames from the UI, if the user has selected the animation plot type
		if(plotToMake == 6):
			self.sv_numFrames = int(self.ui.lineEdit_numFrames.text())

		thread.start_new_thread(self.makePlotSlotT, (plotToMake,self.sv_numFrames))
	
	def runHsweepSlot(self):
		"""This function is called when "OK" is pressed. It reads all of the parameters and calls the runHSweepFun() function, in another thread"""
		#Extract all values from text boxes
		self.dialogToVariablesNonLayer()
		self.storeLayern(self.ui.curr_layer)		
		
		#check that number of threads, current density point count, field point count
		#are all valid values
		if(self.sv_numthreads < 1):
			self.sv_numthreads = 1
			print("Thread Count must be set to something greater than 1 - defaulting to 1")
			Qt.QMessageBox.warning(self, "Error", "Thread Count must be set to something greater than 1 - defaulting to 1")
		if(self.sv_n < 1):
			self.sv_n = 1
			print("Current conditions count must be set to something greater than 1 - defaulting to 1")
			Qt.QMessageBox.warning(self, "Error", "Current conditions count must be set to something greater than 1 - defaulting to 1")
		if(self.sv_nfield < 1):
			self.sv_nfield = 1
			print("External Magnetic Field count must be set to something greater than 1 - defaulting to 1")
			Qt.QMessageBox.warning(self, "Error", "External Magnetic Field count must be set to something greater than 1 - defaulting to 1")
		
		#check that	we've implemented the simtype they want to do
		if self.generateSimType() in self.VALID_SIMTYPES:
			#do nothing
			print("DEBUG: Simtype OK")
		else:
			print("ERROR: Got unsupported simtype %d" % self.generateSimType())
			Qt.QMessageBox.warning(self, "Error", "Invalid Sim Type")
			return
		
		if(print_command_only == 0):
		
			#launch new thread to run the spinsim code and update progress
			if(self.spinsim_thread_lock.acquire(0)):
				#make progress bar visible
				self.ui.progressBar_Progress.setVisible(True)
				self.ui.label_Progress.setVisible(True)
				#reset progress bar
				self.progress_var_lock.acquire()
				self.ui.progressBar_Progress.setValue(0)
				self.progress_var = 0
				self.progress_var_lock.release()
				#start the progress timer
				self.progress_timer.start(1000)
				#spawn the thread
				thread.start_new_thread(self.runHsweepFun, (0,))
					
			#can't run more than one instance of spinsim, sorry
			else:
				print("Sorry, process already running. Wait for current process to finish before restarting.\n")
				Qt.QMessageBox.information(self, "Info", "Sorry, process already running. Wait for current process to finish before restarting.")
				
		else:
			#spawn the thread
			thread.start_new_thread(self.runHsweepFun, (0,))
			
		#print("Returned from runHsweepSlot")
	
	def runHsweepFun(self, garbage):
		"""This function gets called by runHSweepSlot() which is the button press handler. It runs in another thread and periodically updates
			the progress via a locked variable"""
			
		#Ivan: This is the syntax to use
		cmd_str = "%s/spinsim" % (os.getcwd())
		
		#put together simtype based on the current layerEn/layerFixed parameters
		self.sv_simtype = self.generateSimType()
		
		#If we're in 'CGS Mode' switch the units back to SI for just the duration of the cmd_str call
		self.switch_units_to_si(self.cgs_units)
		
		#Ivan : This is the syntax to use
		#Todo: Put in functionality for selecting something other than time-independent current and external field (last 2 params)
		cmd_args = " %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d %e %e %d %e %e %d %e %d %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d %d %e %e %e %e %e %e %e %e %e %e" % (
			self.sv_thetaoff, self.sv_percK2, self.sv_Alpha[0], self.sv_Ms[0], self.sv_Ms[1], self.sv_Ms[2], self.sv_Ms[3], self.sv_K1[0], self.sv_K1[1], self.sv_K1[2], self.sv_K1[3],
			 self.sv_lposz[0], self.sv_lposz[1], self.sv_lposz[2], self.sv_lposz[3], self.sv_simtype, self.sv_jmin, self.sv_jmax, self.sv_n, self.sv_hzmin, self.sv_hzmax,
			 self.sv_nfield, self.sv_tf, self.sv_j_timdep_funs, self.sv_hz_timdep_funs, self.sv_dx, self.sv_dy, self.sv_temperature, self.sv_initmag_x[0], self.sv_initmag_y[0],
			 self.sv_initmag_z[0], self.sv_initmag_x[1], self.sv_initmag_y[1], self.sv_initmag_z[1], self.sv_initmag_x[2], self.sv_initmag_y[2], self.sv_initmag_z[2],
			 self.sv_initmag_x[3], self.sv_initmag_y[3], self.sv_initmag_z[3], self.sv_dz[0], self.sv_dz[1], self.sv_dz[2], self.sv_dz[3], self.sv_pol[0], self.sv_pol[1],
			 self.sv_pol[2], self.sv_pol[3], 0, self.sv_numthreads, self.sv_af[0], self.sv_af[1], self.sv_af[2], self.sv_af[3], self.sv_af[4], self.sv_af[5], self.sv_Alpha[1], self.sv_Alpha[2], self.sv_Alpha[3], self.sv_relax)
			
		print ("Will run the spinsim C code with command:\n%s%s\n" % (cmd_str, cmd_args))	
		
		#switch units back to CGS, if necessary
		self.switch_units_to_cgs(self.cgs_units)
		
		if(print_command_only == 0):
			time_start = time.time()	
			#spawn the spinsim process as a child process and DO NOT BLOCK
			#spinsim_pid = os.spawnv(os.P_NOWAIT, cmd_str, cmd_args)
			spinsim_process = Popen(cmd_str + cmd_args, shell=True, stdout=PIPE, cwd=os.getcwd())
		
			spinsim_pid = spinsim_process.pid
			print("Running spinsim with PID %d" % spinsim_pid)
			
			#loop to poll status of spinsim and update progress bar
			spinsim_running = True
			while(spinsim_running):
				time.sleep(1)
				#print("DEBUG: Polling loop woke up")
				
				#grab stdout
				spinsim_stdout = spinsim_process.stdout
				line = spinsim_stdout.readline()
				line = line.rstrip('\n') #chomp
				#print("DEBUG: Read line %s" % line)
				
				#get the progress from stdout
				progress = re.search('\d+',line)
				if(progress != None):
					#print("Current progress is " + progress.group())
					#set the progress var
					self.progress_var_lock.acquire()
					
					if(int(progress.group()) >= 100):
						self.progress_var = 100
					else:
						self.progress_var = int(progress.group())
						
					self.progress_var_lock.release()
					
					if(int(progress.group()) >= 100):
						spinsim_running = False

				#[/if(progress != None)
				
			#[/while spinsim_running]
			
			time_stop = time.time()
			print("Process took %.1f seconds to run" % (time_stop-time_start))
			
			#done running, release lock
			self.spinsim_thread_lock.release()
			
	#[/runHsweepSlot()]
	
	def handleLoadParamsFromFile(self):
		"""Handler for "load parameters from file" button"""
		#get filename from dialog
		#get filename from dialog
		filename = Qt.QFileDialog.getOpenFileName(self, "", "", "*.ini")
		#if it's a valid length, open the file!
		if(len(filename) > 1):
			self.loadParamsFromFile(str(filename))
	
	def loadParamsFromFile(self, filename):
		"""This function loads the parameters from file 'filename' and updates the GUI with the loaded params
		Intended for either hand-coding parameter file or saving your runs using saveParamsToFile"""
		
		#open the file for reading
		try:
			my_file = open(filename, 'r')
		except IOError:
			print "Couldn't open file! Sorry."
			Qt.QMessageBox.warning(self, "Info", "Couldn't open INI file! Sorry.")
			return 0
		
		#Parse each line based on the first token (tokens split by '=' character)
		#If it's a valid token, read in the second token to the correct variable. If it's a comment, skip.
		#Otherwise, print "unrecognized line" error.
		for line in my_file:
			#split the line on the '=' character
			string_split = line.partition('=')
			
			if(line[0] == '#'):
				print("Skipping comment")
			#match line like "first_token=3.14159"
			elif(string_split[0] == "thetaoff"):
				self.sv_thetaoff = float(string_split[2])
			elif(string_split[0] == "percK2"):
				self.sv_percK2 = float(string_split[2])
			elif(string_split[0] == "Alpha"):
				self.sv_Alpha[0] = float(string_split[2])
				self.sv_Alpha[1] = float(string_split[2])
				self.sv_Alpha[2] = float(string_split[2])
				self.sv_Alpha[3] = float(string_split[2])
				
			elif(string_split[0] == "Alpha_0"):
				self.sv_Alpha[0] = float(string_split[2])
			elif(string_split[0] == "Alpha_1"):
				self.sv_Alpha[1] = float(string_split[2])
			elif(string_split[0] == "Alpha_2"):
				self.sv_Alpha[2] = float(string_split[2])
			elif(string_split[0] == "Alpha_3"):
				self.sv_Alpha[3] = float(string_split[2])
			
			elif(string_split[0] == "K1_1"):
				self.sv_K1[0] = float(string_split[2])
			elif(string_split[0] == "K1_2"):
				self.sv_K1[1] = float(string_split[2])
			elif(string_split[0] == "K1_3"):
				self.sv_K1[2] = float(string_split[2])
			elif(string_split[0] == "K1_4"):
				self.sv_K1[3] = float(string_split[2])
				
			elif(string_split[0] == "Ms_1"):
				self.sv_Ms[0] = float(string_split[2])
			elif(string_split[0] == "Ms_2"):
				self.sv_Ms[1] = float(string_split[2])
			elif(string_split[0] == "Ms_3"):
				self.sv_Ms[2] = float(string_split[2])
			elif(string_split[0] == "Ms_4"):
				self.sv_Ms[3] = float(string_split[2])
				
			elif(string_split[0] == "lposz_1"):
				self.sv_lposz[0] = float(string_split[2])
			elif(string_split[0] == "lposz_2"):
				self.sv_lposz[1] = float(string_split[2])
			elif(string_split[0] == "lposz_3"):
				self.sv_lposz[2] = float(string_split[2])
			elif(string_split[0] == "lposz_4"):
				self.sv_lposz[3] = float(string_split[2])
				
			
			
			elif(string_split[0] == "jmin"):
				self.sv_jmin = float(string_split[2])
			elif(string_split[0] == "jmax"):
				self.sv_jmax = float(string_split[2])
			elif(string_split[0] == "n"):
				self.sv_n = int(string_split[2])
			elif(string_split[0] == "hzmin"):
				self.sv_hzmin = float(string_split[2])
			elif(string_split[0] == "hzmax"):
				self.sv_hzmax = float(string_split[2])
			elif(string_split[0] == "nfield"):
				self.sv_nfield = int(string_split[2])
				
			elif(string_split[0] == "t_f"):
				self.sv_tf = float(string_split[2])
			elif(string_split[0] == "j_timdep_funs"):
				self.sv_j_timdep_funs = int(string_split[2])
			elif(string_split[0] == "hz_timdep_funs"):
				self.sv_hz_timdep_funs = int(string_split[2])
			elif(string_split[0] == "dx"):
				self.sv_dx = float(string_split[2])
			elif(string_split[0] == "dy"):
				self.sv_dy = float(string_split[2])
			elif(string_split[0] == "temperature"):
				self.sv_temperature = float(string_split[2])
				
			elif(string_split[0] == "initmag_x_1"):
				self.sv_initmag_x[0] = float(string_split[2])
			elif(string_split[0] == "initmag_x_2"):
				self.sv_initmag_x[1] = float(string_split[2])
			elif(string_split[0] == "initmag_x_3"):
				self.sv_initmag_x[2] = float(string_split[2])
			elif(string_split[0] == "initmag_x_4"):
				self.sv_initmag_x[3] = float(string_split[2])
				
			elif(string_split[0] == "initmag_y_1"):
				self.sv_initmag_y[0] = float(string_split[2])
			elif(string_split[0] == "initmag_y_2"):
				self.sv_initmag_y[1] = float(string_split[2])
			elif(string_split[0] == "initmag_y_3"):
				self.sv_initmag_y[2] = float(string_split[2])
			elif(string_split[0] == "initmag_y_4"):
				self.sv_initmag_y[3] = float(string_split[2])
				
			elif(string_split[0] == "initmag_z_1"):
				self.sv_initmag_z[0] = float(string_split[2])
			elif(string_split[0] == "initmag_z_2"):
				self.sv_initmag_z[1] = float(string_split[2])
			elif(string_split[0] == "initmag_z_3"):
				self.sv_initmag_z[2] = float(string_split[2])
			elif(string_split[0] == "initmag_z_4"):
				self.sv_initmag_z[3] = float(string_split[2])
				
			elif(string_split[0] == "dz_1"):
				self.sv_dz[0] = float(string_split[2])
			elif(string_split[0] == "dz_2"):
				self.sv_dz[1] = float(string_split[2])
			elif(string_split[0] == "dz_3"):
				self.sv_dz[2] = float(string_split[2])
			elif(string_split[0] == "dz_4"):
				self.sv_dz[3] = float(string_split[2])
				
			elif(string_split[0] == "polarisation_1"):
				self.sv_pol[0] = float(string_split[2])
			elif(string_split[0] == "polarisation_2"):
				self.sv_pol[1] = float(string_split[2])
			elif(string_split[0] == "polarisation_3"):
				self.sv_pol[2] = float(string_split[2])
			elif(string_split[0] == "polarisation_4"):
				self.sv_pol[3] = float(string_split[2])
				
			elif(string_split[0] == "t_relax"):
				self.sv_relax = float(string_split[2])
				
			elif(string_split[0] == "layeren_1"):
				self.sv_layerEn[0] = int(string_split[2])
			elif(string_split[0] == "layeren_2"):
				self.sv_layerEn[1] = int(string_split[2])
			elif(string_split[0] == "layeren_3"):
				self.sv_layerEn[2] = int(string_split[2])
			elif(string_split[0] == "layeren_4"):
				self.sv_layerEn[3] = int(string_split[2])
				
			elif(string_split[0] == "layerfixed_1"):
				self.sv_layerFixed[0] = int(string_split[2])
			elif(string_split[0] == "layerfixed_2"):
				self.sv_layerFixed[1] = int(string_split[2])
			elif(string_split[0] == "layerfixed_3"):
				self.sv_layerFixed[2] = int(string_split[2])
			elif(string_split[0] == "layerfixed_4"):
				self.sv_layerFixed[3] = int(string_split[2])
				
			elif(string_split[0] == "numthreads"):
				self.sv_numthreads = int(string_split[2])
			elif(string_split[0] == "numFrames"):
				self.sv_numFrames = int(string_split[2])
				
			elif(string_split[0] == "af_bw_12"):
				self.sv_af[0] = float(string_split[2])
			elif(string_split[0] == "af_bw_13"):
				self.sv_af[1] = float(string_split[2])
			elif(string_split[0] == "af_bw_14"):
				self.sv_af[2] = float(string_split[2])
			elif(string_split[0] == "af_bw_23"):
				self.sv_af[3] = float(string_split[2])
			elif(string_split[0] == "af_bw_24"):
				self.sv_af[4] = float(string_split[2])
			elif(string_split[0] == "af_bw_34"):
				self.sv_af[5] = float(string_split[2])
				
			else:
				print("Couldn't parse line %s in file %s, are you sure the format is 'token=123.0' with no extra spaces?" % (line, filename))
		
		#close the file
		my_file.close()
		
		#load new parameters into GUI fields
		self.variablesToDialogNonLayer()
		self.handleLayernButton(self.ui.curr_layer)
		
		return 0
	#[/loadParamsFromFile]
	
	def handleSaveParamsToFile(self):
		"""Handler for "save parameters to file" button"""
		#get filename from dialog
		filename = Qt.QFileDialog.getSaveFileName(self, "", "", "*.ini")
		#if it's a valid length, open the file!
		if(len(filename) > 1):
			self.saveParamsToFile(str(filename))

	def saveParamsToFile(self, filename):
		"""This function saves the parameters to file 'filename'
		Intended for saving your runs, to be used later with loadParamsFromFile"""
		
		#grab latest parameters from GUI fields
		self.dialogToVariablesNonLayer()
		self.storeLayern(self.ui.curr_layer)
		
		print("Debug: Saving params to file %s" % filename)
		
		#Open file for writing
		try:
			my_file = open(filename, 'w')
		except IOError:
			print ("Couldn't open file %s! Sorry." % my_file)
			Qt.QMessageBox.warning(self, "Info", "Couldn't open INI file! Sorry.")
			return 0
		
		#write out all of the parameters, with the correct token string so that
		#loadParamsFromFile() can read and parse the file back in.
		my_file.write("thetaoff=%e\n" % self.sv_thetaoff)
		my_file.write("percK2=%e\n" % self.sv_percK2)
		my_file.write("Alpha_0=%e\n" % self.sv_Alpha[0])
		my_file.write("Alpha_1=%e\n" % self.sv_Alpha[1])
		my_file.write("Alpha_2=%e\n" % self.sv_Alpha[2])
		my_file.write("Alpha_3=%e\n" % self.sv_Alpha[3])
		
		my_file.write("K1_1=%e\n" % self.sv_K1[0])
		my_file.write("K1_2=%e\n" % self.sv_K1[1])
		my_file.write("K1_3=%e\n" % self.sv_K1[2])
		my_file.write("K1_4=%e\n" % self.sv_K1[3])
			
		my_file.write("Ms_1=%e\n" % self.sv_Ms[0])
		my_file.write("Ms_2=%e\n" % self.sv_Ms[1])
		my_file.write("Ms_3=%e\n" % self.sv_Ms[2])
		my_file.write("Ms_4=%e\n" % self.sv_Ms[3])
			
		my_file.write("lposz_1=%e\n" % self.sv_lposz[0])
		my_file.write("lposz_2=%e\n" % self.sv_lposz[1])
		my_file.write("lposz_3=%e\n" % self.sv_lposz[2])
		my_file.write("lposz_4=%e\n" % self.sv_lposz[3])
		
		my_file.write("jmin=%e\n" % self.sv_jmin)
		my_file.write("jmax=%e\n" % self.sv_jmax)
		my_file.write("n=%d\n" % self.sv_n)
		my_file.write("hzmin=%e\n" % self.sv_hzmin)
		my_file.write("hzmax=%e\n" % self.sv_hzmax)
		my_file.write("nfield=%d\n" % self.sv_nfield)
			
		my_file.write("t_f=%e\n" % self.sv_tf)
		my_file.write("j_timdep_funs=%d\n" % self.sv_j_timdep_funs)
		my_file.write("hz_timdep_funs=%d\n" % self.sv_hz_timdep_funs)
		my_file.write("dx=%e\n" % self.sv_dx)
		my_file.write("dy=%e\n" % self.sv_dy)
		my_file.write("temperature=%e\n" % self.sv_temperature)
			
		my_file.write("initmag_x_1=%e\n" % self.sv_initmag_x[0])
		my_file.write("initmag_x_2=%e\n" % self.sv_initmag_x[1])
		my_file.write("initmag_x_3=%e\n" % self.sv_initmag_x[2])
		my_file.write("initmag_x_4=%e\n" % self.sv_initmag_x[3])
			
		my_file.write("initmag_y_1=%e\n" % self.sv_initmag_y[0])
		my_file.write("initmag_y_2=%e\n" % self.sv_initmag_y[1])
		my_file.write("initmag_y_3=%e\n" % self.sv_initmag_y[2])
		my_file.write("initmag_y_4=%e\n" % self.sv_initmag_y[3])
			
		my_file.write("initmag_z_1=%e\n" % self.sv_initmag_z[0])
		my_file.write("initmag_z_2=%e\n" % self.sv_initmag_z[1])
		my_file.write("initmag_z_3=%e\n" % self.sv_initmag_z[2])
		my_file.write("initmag_z_4=%e\n" % self.sv_initmag_z[3])
			
		my_file.write("dz_1=%e\n" % self.sv_dz[0])
		my_file.write("dz_2=%e\n" % self.sv_dz[1])
		my_file.write("dz_3=%e\n" % self.sv_dz[2])
		my_file.write("dz_4=%e\n" % self.sv_dz[3])
			
		my_file.write("polarisation_1=%e\n" % self.sv_pol[0])
		my_file.write("polarisation_2=%e\n" % self.sv_pol[1])
		my_file.write("polarisation_3=%e\n" % self.sv_pol[2])
		my_file.write("polarisation_4=%e\n" % self.sv_pol[3])
			
		my_file.write("layeren_1=%d\n" % self.sv_layerEn[0])
		my_file.write("layeren_2=%d\n" % self.sv_layerEn[1])
		my_file.write("layeren_3=%d\n" % self.sv_layerEn[2])
		my_file.write("layeren_4=%d\n" % self.sv_layerEn[3])
			
		my_file.write("layerfixed_1=%d\n" % self.sv_layerFixed[0])
		my_file.write("layerfixed_2=%d\n" % self.sv_layerFixed[1])
		my_file.write("layerfixed_3=%d\n" % self.sv_layerFixed[2])
		my_file.write("layerfixed_4=%d\n" % self.sv_layerFixed[3])
			
		my_file.write("numthreads=%d\n" % self.sv_numthreads)
		my_file.write("numFrames=%d\n" % self.sv_numFrames)

		my_file.write("af_bw_12=%e\n" % self.sv_af[0]);
		my_file.write("af_bw_13=%e\n" % self.sv_af[1]);
		my_file.write("af_bw_14=%e\n" % self.sv_af[2]);
		my_file.write("af_bw_23=%e\n" % self.sv_af[3]);
		my_file.write("af_bw_24=%e\n" % self.sv_af[4]);
		my_file.write("af_bw_34=%e\n" % self.sv_af[5]);
		
		my_file.write("t_relax=%e\n" % self.sv_relax);
		
		my_file.close()
		return 0	
	#[/saveParamsToFile]
	
	def generateSimType(self):
		"""This function reads sv_layerEn and sv_layerFixed variables and packs them
		into the standard 'simtype' form. See the documentation for more information about
		what the "generic" models that we pack down into are. Mostly, it's just moving the layers
		so that the enabled layers are all adjacent to each other."""
		
		layerPos = [0,0,0,0]
		count = 0
		lp_count = 0
		#go through all enabled layers and count them
		while(count < 4):
			if(self.sv_layerEn[count] != 0):
				layerPos[lp_count] = count
				lp_count += 1
			count += 1
		
		#put together integer to return
		returned = 0
		count = 0
		while(count < lp_count):
			if(self.sv_layerFixed[layerPos[count]] == 1):
				returned += (1 << (2*count))
			else:
				returned += (2 << (2*count))
			count += 1
		#[/while]
		
		#set field-in-plane simulation function
		returned += self.sv_fieldInPlane << 31
		
		return(returned)
	
	def timedep_writeback(self, j_funs, hz_funs):
		"""This function gets called by TimeDep_Handler to pass the j_timdep_funs and
		hz_timdep_funs variables back to (this) parent dialog"""
		self.sv_j_timdep_funs = j_funs
		self.sv_hz_timdep_funs = hz_funs
		print("DEBUG: spinsim GUI parent dialog got %f, %f for j_funs and hz_funs" % (self.sv_j_timdep_funs, self.sv_hz_timdep_funs))

#Ivan: Removed spinsim_make_plots() 2010-02-16, code can be found in revision dated 2010-02-14

#Launch the GUI
if __name__ == "__main__":
	app = QtGui.QApplication(sys.argv)
	myapp = MyForm()
	myapp.show()
	app.exec_()

