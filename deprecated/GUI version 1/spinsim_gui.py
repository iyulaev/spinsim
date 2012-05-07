#Spinsim_GSL GUI
#GUI interface for Spinsim_GSL command-line program
#Written I. Yulaev (iyulaev@ucsd.edu) 2010-01-16
#University of California, San Diego ECE Department
#Magnet icon taken from openoffice.org website

import time
import re
import thread

#Stuff for PyQt
import sys
from PyQt4 import QtCore, QtGui
from qt_ui import Ui_Dialog

#stuff for command line calls
from subprocess import *
import os

#stuff for creating the plots
from numpy import *
from pylab import *
import matplotlib.pyplot as plt

#TODO: Add to tree diagram: handleMakePlotT, external process call for that, handleSetLimits, button for that

class MyForm(QtGui.QMainWindow):
	"""This class initializes the GUI and also contains code for reading in the parameters and running the spinsim C code"""
	def __init__(self, parent=None):
		""" This function initializes the GUI and sets up all of the default parameters """
		QtGui.QWidget.__init__(self, parent)
		self.ui = Ui_Dialog()
		self.ui.setupUi(self)

		#set up all of the "default" values
		self.ui.lineEdit_K1.setText("3e5")
		self.ui.lineEdit_Alpha.setText("0.01")
		self.ui.lineEdit_Ms.setText("650000");
		self.ui.lineEdit_Mrs.setText("600000");
		self.ui.lineEdit_dx.setText("50");
		self.ui.lineEdit_dy.setText("50");
		self.ui.lineEdit_dz.setText("3");
		self.ui.lineEdit_jmax.setText("1e11");
		self.ui.lineEdit_hzmax.setText("1e6");
		self.ui.lineEdit_jmin.setText("-1e11");
		self.ui.lineEdit_hzmin.setText("-1e6");
		self.ui.lineEdit_t_f.setText("2.5e-8");
		self.ui.lineEdit_perK2.setText("0.5");
		self.ui.lineEdit_thetaoff.setText("0.05");
		self.ui.lineEdit_n.setText("30");
		self.ui.lineEdit_nfield.setText("30");
		self.ui.lineEdit_threadcount.setText("1");
		
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
		
		#setup timer for updating progress bar
		self.progress_timer = QtCore.QTimer()
		QtCore.QObject.connect(self.progress_timer, QtCore.SIGNAL("timeout()"), self.handleUpdateProgress)
		
		#setup the comboBox (and its label) for plotting
		self.ui.comboBox_PlotType.setEditable(False)
		self.ui.comboBox_PlotType.setVisible(True)
		#self.ui.label_GeneratePlot.setVisible(False)
		#self.ui.pushButton_MakePlot.setVisible(False)
		self.plotcount = 0
		
		
	#[/__init__()]
	
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
			j and Hz limits based on the last zoom-box limit settings from the last figure call"""
			
		#check to see if our plot limits were updated. If so, post the updates to the GUI
		self.spinsim_plot_lock.acquire()
		
		if(self.plot_lims_updated == 1):
			self.ui.lineEdit_jmax.setText(str(self.plot_lims[3]));
			self.ui.lineEdit_hzmax.setText(str(self.plot_lims[1]));
			self.ui.lineEdit_jmin.setText(str(self.plot_lims[2]));
			self.ui.lineEdit_hzmin.setText(str(self.plot_lims[0]));
			self.plot_lims_updated = 0
					
		self.spinsim_plot_lock.release()
		
	def makePlotSlotT(self,plotToMake):
		"""This function calls spinsim_make_plots.py via the shell to create the figure window and draw our plots. It also creates
			*.png file with the plots in the file"""
			
		cmdstr = "python spinsim_make_plots.py %d" % (plotToMake+1)
		
		#create plots by calling a new process
		plot_process = Popen(cmdstr, shell=True, stdout=PIPE, cwd=os.getcwd())
		
		#loop to poll status of plot and update limits
		while(plot_process.poll() == None):
			time.sleep(1)
			
			#grab stdout
			plot_stdout = plot_process.stdout
			#Ivan: Why does this line only return when we exit the program?
			#because I wasn't flushing stdout in the other process...duh...s
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
		
	def makePlotSlot(self):
		"""This function calls makePlotSlotT in a new thread"""
		plotToMake = self.ui.comboBox_PlotType.currentIndex()

		thread.start_new_thread(self.makePlotSlotT, (plotToMake,))
		
		#spinsim_make_plots(plotToMake,self.plotcount)
		self.plotcount += 1
		#self.ui.pushButton_MakePlot.setVisible(False)
		
	def handleSinglePlot(self):
		"""This function handler calls a new thread to call a process to create a figure showing the single-point time domain solution"""
		#not sure why we NEED to pass a tuple for parameters
		thread.start_new_thread(self.handleSinglePlotT, (0,))
		
	def handleSinglePlotT(self, garbage):
		"""This function calls spinsim_make_plots.py via the shell to create the figure window and draw our single-point plot.
			It also creates *.png file with the plots in the file"""
		
		#parameter garbage isn't used
		
		cmdstr = "python spinsim_make_plots.py 6"
		
		#create plots by calling a new process
		plot_process = Popen(cmdstr, shell=True, cwd=os.getcwd())
	
	def runHsweepSlot(self):
		"""This function is called when "OK" is pressed. It reads all of the parameters and calls the runHSweepFul() function, in another thread"""
		#Extract all the values of out the text boxes
		k1_float = float(self.ui.lineEdit_K1.text())
		alpha_float = float(self.ui.lineEdit_Alpha.text())
		ms_float = float(self.ui.lineEdit_Ms.text())
		mrs_float = float(self.ui.lineEdit_Mrs.text())
		dx_float = float(self.ui.lineEdit_dx.text())
		dy_float = float(self.ui.lineEdit_dy.text())
		dz_float = float(self.ui.lineEdit_dz.text())
		jmax_float = float(self.ui.lineEdit_jmax.text())
		hzmax_float = float(self.ui.lineEdit_hzmax.text())
		jmin_float = float(self.ui.lineEdit_jmin.text())
		hzmin_float = float(self.ui.lineEdit_hzmin.text())
		t_f_float = float(self.ui.lineEdit_t_f.text())
		perK2_float = float(self.ui.lineEdit_perK2.text())
		thetaoff_float = float(self.ui.lineEdit_thetaoff.text())
		n_int = int(self.ui.lineEdit_n.text())
		nfield_int = int(self.ui.lineEdit_nfield.text())
		threadcount_int = int(self.ui.lineEdit_threadcount.text())
		
		if(threadcount_int < 1):
			threadcount_int = 1
			print("Thread Count must be set to something greater than 1 - defaulting to 1")
		if(n_int < 1):
			n_int = 1
			print("Current conditions count must be set to something greater than 1 - defaulting to 1")
		if(nfield_int < 1):
			nfield_int = 1
			print("External Magnetic Field count must be set to something greater than 1 - defaulting to 1")
		
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
			thread.start_new_thread(self.runHsweepFun, (k1_float, alpha_float, ms_float, mrs_float, dx_float, dy_float, dz_float,
				jmin_float, jmax_float, hzmin_float, hzmax_float, t_f_float, perK2_float, thetaoff_float, n_int, nfield_int, threadcount_int))
				
		#can't run more than one instance of spinsim, sorry
		else:
			print("Sorry, process already running. Wait for current process to finish before restarting.\n")
			
		#print("Returned from runHsweepSlot")
	
	def runHsweepFun(self, k1_float, alpha_float, ms_float, mrs_float, dx_float, dy_float, dz_float,
		jmin_float, jmax_float, hzmin_float, hzmax_float, t_f_float, perK2_float, thetaoff_float, n_int, nfield_int, threadcount_int):
		"""This function gets called by runHSweepSlot() which is the button press handler. It runs in another thread and periodically updates
			the progress via a locked variable"""
			
		#debug print
		#cmd_str = "%s/spinsim %e %f %f %f %f %f %f %e %e %e %f %f %d %d %d > output.spinsim" % (os.getcwd(), k1_float,
		#	alpha_float, ms_float, mrs_float, dx_float, dy_float, dz_float, jmax_float, hzmax_float, t_f_float, perK2_float,
		#	thetaoff_float, n_int, nfield_int, threadcount_int)
		#Ivan: This is the syntax to use
		cmd_str = "%s/spinsim" % (os.getcwd())
		
		#cmd_args = ["spinsim", str(k1_float), str(alpha_float), str(ms_float), str(mrs_float), str(dx_float),
		#	str(dy_float), str(dz_float), str(jmax_float), str(hzmax_float), str(t_f_float), str(perK2_float),
		#	str(thetaoff_float), str(n_int), str(nfield_int), str(threadcount_int)]
		#Ivan : This is the syntax to use
		#Todo: Put in functionality for selecting something other than time-independent current and external field (last 2 params)
		cmd_args = " %e %f %f %f %f %f %f %e %e %e %e %e %f %f %d %d %d 1 1" % (k1_float,
			alpha_float, ms_float, mrs_float, dx_float, dy_float, dz_float, jmin_float, jmax_float, hzmin_float, hzmax_float,
			t_f_float, perK2_float,	thetaoff_float, n_int, nfield_int, threadcount_int)
			
		#print ("Will run the spinsim C code with command:\n%s\n" % cmd_str)	
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
	
	def loadParamsFromFile(self, filename):
		"""This function loads the parameters from file 'filename' and updates the GUI with the loaded params
		Intended for either hand-coding parameter file or saving your runs using saveParamsToFile"""
		try:
			my_file = open(filename, 'r')
			("Couldn't open file %s! Sorry." % my_file)
		except IOError:
			print "Couldn't open file! Sorry."
			return 0
		
		for line in my_file:
			#split the line on the '=' character
			string_split = line.partition('=')
			
			if(line[0] == '#'):
				print("Skipping comment")
			#match line like "first_token=3.14159"
			elif(string_split[0] == "first_token"):
				self.first_token = float(string_split[2])
			elif(string_split[0] == "second_token"):
				self.second_token = float(string_split[2])
			else:
				print("Couldn't parse line %s in file %s, are you sure the format is 'token=123.0' with no extra spaces?" % (line, filename))
		
		my_file.close()
		return 0
	#[/loadParamsFromFile]

	def saveParamsToFile(self, filename):
		"""This function saves the parameters to file 'filename'
		Intended for saving your runs, to be used later with loadParamsFromFile"""
		try:
			my_file = open(filename, 'w')
		except IOError:
			print ("Couldn't open file %s! Sorry." % my_file)
			return 0
		
		my_file.write('first_token=%e\n' % self.first_token)
		my_file.write('second_token=%e\n' % self.second_token)
		
		my_file.close()
		return 0	
	#[/saveParamsToFile]

#Ivan: Removed spinsim_make_plots() 2010-02-16, code can be found in revision dated 2010-02-14

#Launch the GUI
if __name__ == "__main__":
	app = QtGui.QApplication(sys.argv)
	myapp = MyForm()
	myapp.show()
	app.exec_()

