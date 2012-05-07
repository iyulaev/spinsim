#plotting libraries
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#system library for getting command line args
import sys

import time

#stuff for command line calls
from subprocess import *
import os

old_lims = []

def handle_axes_manip(event):
	#print "Got a axes_enter_event!"
	
	if(event.inaxes != None):
		if((event.inaxes.get_xlim()[0] == 0) & (event.inaxes.get_xlim()[1] == 1)):
			return
		if((old_lims[0],old_lims[1],old_lims[2],old_lims[3]) != (event.inaxes.get_xlim()[0], event.inaxes.get_xlim()[1], event.inaxes.get_ylim()[0], event.inaxes.get_ylim()[1])):
			print("N-xlim-ylim = %f %f %f %f\n" % (event.inaxes.get_xlim()[0], event.inaxes.get_xlim()[1], event.inaxes.get_ylim()[0], event.inaxes.get_ylim()[1]))
			sys.stdout.flush()
			old_lims[0] = event.inaxes.get_xlim()[0]
			old_lims[1] = event.inaxes.get_xlim()[1]
			old_lims[2] = event.inaxes.get_ylim()[0]
			old_lims[3] = event.inaxes.get_ylim()[1]

def spinsim_make_plots(plottype):
	"""This function generates the magsweep plots from data at a pre-determined location. On plottype==5, generates all four plots: 
	spin up, spin down, sum of two, difference of two. Alternately, it does spinup, spindown, sum, difference for each of the plottypes 1-4.
	It also makes a PNG (150dpi) for each plot that gets generated."""
	
	#this variable isn't really used anymore
	plotcount = 0
	
	#print("Generating plot of type %d" % plottype)
	
	#This closes the current figure, but doesn't seem to help with opening another one
	#close('all')
	
	#define all of the names of the data files
	filename_up = "magsweepup_output"
	filename_down = "magsweepdown_output"
	jfile = 'magsweep_output_j'
	Hzfile = 'magsweep_output_hz'
	
	#f = open(filename, 'r')
	#load in all of the data files
	try:
		data_up = loadtxt(filename_up, unpack=True)
		data_down = loadtxt(filename_down, unpack=True)
		j = loadtxt(jfile, unpack=True)
		hz = loadtxt(Hzfile,unpack=True)
	except IOError as e:
		print("ERROR: Got an IO Error trying to read the files. Are you sure you've already run the spinsim tool?\nError was %s" % e)

	#added 2010-01-30
	#sets up a new figure, with four sub-plots
	fig = plt.figure(plotcount, figsize=(16,10), facecolor='w')
	fig.subplots_adjust(left=0.1, wspace=0.2)
	#link axes_enter_event with the handler for that event
	fig.canvas.mpl_connect('axes_enter_event', handle_axes_manip)
	fig.canvas.mpl_connect('axes_leave_event', handle_axes_manip)
	
	#remember the old axes limits, cause we're going to have to see if we need to change them or not
	old_lims.append(hz[0])
	old_lims.append(hz[len(hz)-1])
	old_lims.append(j[0])
	old_lims.append(j[len(j)-1])
	
	#Set up subplot for sweep up
	if((plottype == 1) | (plottype == 5)):	
		ax1 = fig.add_subplot(221 if (plottype==5) else 111)
		ax1.set_title('Magnetization - Sweep Up')
		ax1.set_xlabel('Applied Field (A/m)')
		ax1.set_ylabel('Applied Current (A/m**2)')
		#ax1.pcolor(X,Y,data_up.T)
		#imshow creates the actual plot
		#note that flipud() is used to flip the data upside down. For some reason imshow likes to plot the data upside down...
		cax1 = ax1.imshow(flipud(data_up.T), aspect='auto', interpolation='nearest',extent=[hz[0], hz[len(hz)-1], j[0], j[len(j)-1]])
		#colorbar creates the color bar (legend)
		cbar = fig.colorbar(cax1)
	
	if((plottype == 2) | (plottype == 5)):	
		#Set up subplot for sweep down
		ax2 = fig.add_subplot(222 if (plottype==5) else 111)
		ax2.set_title('Magnetization - Sweep Down')
		ax2.set_xlabel('Applied Field (A/m)')
		ax2.set_ylabel('Applied Current (A/m**2)')
		#ax2.pcolor(X,Y,data_down.T)
		cax2 = ax2.imshow(flipud(data_down.T), aspect='auto', interpolation='nearest',extent=[hz[0], hz[len(hz)-1], j[0], j[len(j)-1]])
		cbar = fig.colorbar(cax2)
	
	#Set up subplot for sum of the two
	if((plottype == 4) | (plottype == 5)):	
		ax3 = fig.add_subplot(223 if (plottype==5) else 111)
		ax3.set_title('Magnetization - (Up + Down)')
		ax3.set_xlabel('Applied Field (A/m)')
		ax3.set_ylabel('Applied Current (A/m**2)')
		#ax3.pcolor(X,Y,(data_up.T + data_down.T))
		cax3 = ax3.imshow(flipud(data_up.T + data_down.T), aspect='auto', interpolation='nearest',extent=[hz[0], hz[len(hz)-1], j[0], j[len(j)-1]])
		cbar = fig.colorbar(cax3)
	
	if((plottype == 3) | (plottype == 5)):	
		#Set up subplot for difference of the two
		ax4 = fig.add_subplot(224 if (plottype==5) else 111)
		ax4.set_title('Magnetization - (Up - Down)')
		ax4.set_xlabel('Applied Field (A/m)')
		ax4.set_ylabel('Applied Current (A/m**2)')
		#ax4.pcolor(X,Y,(data_up.T - data_down.T))
		cax4 = ax4.imshow(flipud(data_up.T - data_down.T), aspect='auto', interpolation='nearest',extent=[hz[0], hz[len(hz)-1], j[0], j[len(j)-1]])
		cbar = fig.colorbar(cax4)
	
	#X, Y = meshgrid(hz,j)
	
	#make up a filename to set the pltos to
	filestring = "magsweep_plots.png"
	if(plottype == 1):
		filestring = "magsweep_up_plot.png"
	if(plottype == 2):
		filestring = "magsweep_down_plot.png"
	if(plottype == 3):
		filestring = "magsweep_sum_plot.png"
	if(plottype == 4):
		filestring = "magsweep_diff_plot.png"
	
	#save the figure to a file
	plt.savefig(filestring, facecolor='w', edgecolor='w', orientation='landscape', dpi=150)
	
	#plt.draw()
	plt.show()
	
	#plt.close('all')
#[/spinsim_make_plots]
	
def spinsim_make_singleplot(plottype, numlays):
	"""This function generates several plots for a single point solution. It generates x,y,z magnetization vs
	time (on separate plots) and also draws a 3D plot of the magnetization in an x,y,z co-ordinate space"""
	#define file names
	filename_t = "magsingle_t"
	filename_x = "magsingle_x"
	filename_y = 'magsingle_y'
	filename_z = 'magsingle_z'
	
	if(numlays == 2):
		filename_t2 = "magsingle_t2"
		filename_x2 = "magsingle_x2"
		filename_y2 = 'magsingle_y2'
		filename_z2 = 'magsingle_z2'
	
	#try to load the files in
	try:
		data_t = loadtxt(filename_t, unpack=True)
		data_x = loadtxt(filename_x, unpack=True)
		data_y = loadtxt(filename_y, unpack=True)
		data_z = loadtxt(filename_z, unpack=True)
	except IOError as e:
		print("ERROR: Got an IO Error trying to read the files. Are you sure you've already run the spinsim tool for single-point?\nError was %s" % e)
		data_t = None;
		
	if(numlays == 2):
		try:
			data_t2 = loadtxt(filename_t2, unpack=True)
			data_x2 = loadtxt(filename_x2, unpack=True)
			data_y2 = loadtxt(filename_y2, unpack=True)
			data_z2 = loadtxt(filename_z2, unpack=True)
		except IOError as e:
			print("ERROR: Got an IO Error trying to read the files. Are you sure you've already run the spinsim tool for single-point?\nError was %s" % e)
			data_t2 = None;
		
	#probably need to decimate the data if it's too big -> figure that 50,000 points is too big
	if(data_t != None):
		while (len(data_t) > 25000):
			data_t = data_t[1:len(data_t):3]	
			data_x = data_x[1:len(data_x):3]
			data_y = data_y[1:len(data_y):3]
			data_z = data_z[1:len(data_z):3]
	#[/while]
	if(numlays > 1 and data_t2 != None):
		while (len(data_t2) > 25000):
			data_t2 = data_t2[1:len(data_t):3]	
			data_x2 = data_x2[1:len(data_x):3]
			data_y2 = data_y2[1:len(data_y):3]
			data_z2 = data_z2[1:len(data_z):3]
		#[/while]
	
	#check for all-same in x,y,z data
	#don't worry about this check if we're doing two layers because holy fuck
	#i'm lazy
	if(data_t != None):
		all_same = 1
		first_val = data_x[0]
		for point in data_x:
			if(point != first_val):
				all_same = 0
		if(all_same == 1):
			data_x[0] = data_x[0] - 0.01
			
		all_same = 1
		first_val = data_y[0]
		for point in data_y:
			if(point != first_val):
				all_same = 0
		if(all_same == 1):
			data_y[0] = data_y[0] - 0.01
			
		all_same = 1
		first_val = data_z[0]
		for point in data_z:
			if(point != first_val):
				all_same = 0
		if(all_same == 1):
			data_z[0] = data_z[0] - 0.01
	
		
	#sets up a new figure, with four sub-plots
	fig = plt.figure(0, figsize=(16,9), facecolor='w')
	fig.subplots_adjust(left=0.1, wspace=0.2, hspace=0.5)
	#link axes_enter_event with the handler for that event
	#we're not supporting re-plot by zoom so whatever
	#fig.canvas.mpl_connect('axes_enter_event', handle_axes_manip)
	#fig.canvas.mpl_connect('axes_leave_event', handle_axes_manip)
	
	#put together the X versus T plot
	sp1 = fig.add_subplot(311)
	if(data_t!=None):
		sp1.plot(data_t, data_x)
	if(numlays==2 and data_t2!=None):
		sp1.plot(data_t2, data_x2, 'r--')
	sp1.grid(True)
	sp1.set_title('X versus T')
	sp1.set_ylabel('X component of magnetization')
	sp1.set_xlabel('Time (s)')
	
	#put together the Y versus T plot
	sp2 = fig.add_subplot(312)
	if(data_t!=None):
		sp2.plot(data_t, data_y)
	if(numlays==2 and data_t2!=None):
		sp2.plot(data_t2, data_y2, 'r--')
	sp2.grid(True)
	sp2.set_title('Y versus T')
	sp2.set_ylabel('Y component of magnetization')
	sp2.set_xlabel('Time (s)')
	
	#put together the Z versus T plot
	sp3 = fig.add_subplot(313)
	if(data_t!=None):	
		sp3.plot(data_t, data_z)
	if(numlays==2 and data_t2!=None):
		sp3.plot(data_t2, data_z2, 'r--')
	sp3.grid(True)
	sp3.set_title('Z versus T')
	sp3.set_ylabel('Z component of magnetization')
	sp3.set_xlabel('Time (s)')
	
	#save the figure to a file
	plt.savefig('magsingle_xyz.png', facecolor='w', edgecolor='w', orientation='landscape', dpi=150)
	
	#put together the 3D plot
	#see http://matplotlib.sourceforge.net/examples/mplot3d/lines3d_demo.html for the code on which
	#this is based
	fig2 = plt.figure(facecolor='w')
	ax_3d = Axes3D(fig2)
	if(data_t != None):
		ax_3d.plot(data_x, data_y, data_z, label='3D visualization of L1 magnetization vector movement')
	if(data_t2 != None):
		ax_3d.plot(data_x2, data_y2, data_z2,  'r--',label='3D visualization of L2 magnetization vector movement')
	#ax_3d.plot([0.9,1],[0,1],[0.9,1], label='3D visualization of magnetization vector movement')
	ax_3d.legend()
	#save the 3D plot
	plt.savefig('magsingle_3d.png', facecolor='w', edgecolor='w', orientation='landscape', dpi=150)
	
	#plt.draw()
	plt.show()
#[/spinsim_make_singleplot]

def spinsim_make_xyplot(plottype):
	"""This function generates x-y plots for the magnetization, based on either an H(z) external field sweep or a J current density sweep. 
	Plottype 8 will give H(z) vs magnetization plot, plottype 9 will give J vs magnetization plot."""

	#define file names
	
	if(plottype == 8) :
		filename_x = "magsweep_output_hz"
		filename_y1 = "magsweepdown_output"
		filename_y2 = "magsweepup_output"
		filename_out = "magsweep_h-sweep.png"
	if(plottype == 9) :
		filename_x = "magsweep_output_j"
		filename_y1 = "magsweepdown_output"
		filename_y2 = "magsweepup_output"
		filename_out = "magsweep_j-sweep.png"
	
	#try to load the files in
	try:
		data_x = loadtxt(filename_x, unpack=True)
		data_y1 = loadtxt(filename_y1, unpack=True)
		data_y2 = loadtxt(filename_y2, unpack=True)
	except IOError as e:
		print("ERROR: Got an IO Error trying to read the files. Are you sure you've already run the spinsim tool for single-point?\nError was %s" % e)
		
	#if(plottype == 9): 
	#	data_y1 = data_y1.transpose()
	#	data_y2 = data_y2.transpose()
		
	#sets up a new figure, with four sub-plots
	fig = plt.figure(0, figsize=(16,9), facecolor='w')
	fig.subplots_adjust(left=0.1, wspace=0.2, hspace=0.5)
	#link axes_enter_event with the handler for that event
	#we're not supporting re-plot by zoom so whatever
	#fig.canvas.mpl_connect('axes_enter_event', handle_axes_manip)
	#fig.canvas.mpl_connect('axes_leave_event', handle_axes_manip)
	
	#put together the X versus T plot
	sp1 = fig.add_subplot(111)
	sp1.plot(data_x, data_y1)
	sp1.plot(data_x, data_y2, 'r')
	sp1.set_ylim(ymin=-1.1)
	sp1.set_ylim(ymax=1.1)
	sp1.grid(True)
	if(plottype == 8): sp1.set_title('Magnetization versus H(z)')
	if(plottype == 9): sp1.set_title('Magnetization versus J')
	sp1.set_ylabel('Magnetization - red is sweep-up, blue sweep-down')
	if(plottype == 8): sp1.set_xlabel('H(z)')
	if(plottype == 9): sp1.set_xlabel('J')
	
	#save the figure to a file
	plt.savefig(filename_out, facecolor='w', edgecolor='w', orientation='landscape', dpi=75)
	
	#plt.draw()
	plt.show()
#[/spinsim_make_xyplot]

def spinsim_make_multiplot(plottype, NUM_FRAMES):
	"""This function generates a series of plots with more and more data. If one were to
	take these plots and "play them" with a short time delay, it'd look like an animation. Parameter 
	plottype is discarded. Parameter NUM_FRAMES determines how many frames of the animation you'll 
	get. It should be a positive, non-zero integer."""
	
	#define file names
	filename_t = "magsingle_t"
	filename_x = "magsingle_x"
	filename_y = 'magsingle_y'
	filename_z = 'magsingle_z'
	
	#try to load the files in
	try:
		data_t = loadtxt(filename_t, unpack=True)
		data_x = loadtxt(filename_x, unpack=True)
		data_y = loadtxt(filename_y, unpack=True)
		data_z = loadtxt(filename_z, unpack=True)
	except IOError as e:
		print("ERROR: Got an IO Error trying to read the files. Are you sure you've already run the spinsim tool for single-point?\nError was %s" % e)
	
	count = 1
	
	#print("Length of data is %d" % (data_t.size))
	#print("Array data 1/10th is")
	#print(data_t[0:(data_t.size/10)])
	#print("length of that is %d" % ((data_t[0:(data_t.size/10)]).size))
	
	#change this to change the number of frames in each animation
	fig_list = []
	ax_list = []
	
	#put together the frames one at a time
	while(count < NUM_FRAMES):
		fig_list.append(plt.figure(count-1, figsize=(8,8), facecolor='w'))
		ax_list.append(Axes3D(fig_list[count-1]))
		
		ax_list[count-1].plot(data_x[0:(data_x.size / NUM_FRAMES * count)], \
			data_y[0:(data_y.size / NUM_FRAMES * count)],\
			data_z[0:(data_z.size / NUM_FRAMES * count)],\
			label='3D visualization of magnetization vector movement')
			
		ax_list[count-1].legend()
		
		#save as a unique file
		filename_str = "magsingle_3d_%d.png" % (count)
		#IY: on 2010-02-24, reduced DPI to 50 from 150 to make the images smaller.
		plt.savefig(filename_str, facecolor='w', edgecolor='w', orientation='landscape', dpi=50)
		
		count += 1
	#[/while]
	
	print ("DEBUG: About to run apngasm")
	
	#run apngasm
	#create plots by calling a new process
	cmdstr = "./apngasm magsingle_anim.png magsingle_3d_1.png"
	apngasm_process = Popen(cmdstr, shell=True, cwd=os.getcwd())
	
	#wait for apngasm to finish
	while(apngasm_process.poll() == None):
			time.sleep(1)
			
	print ("DEBUG: apngasm finished, deleting files!")
	
	#delete all of the intermediate old files
	count = 1
	while(count < NUM_FRAMES):
		os.unlink("magsingle_3d_%d.png" % count)
		count += 1
		
#[/spinsim_make_multiplot]
	
#main
if((int(sys.argv[1]) > 9) | (int(sys.argv[1]) < 1) ):
	print("ERROR: spinsim_make_plots.py got invalid command line argument")
	print("Got %s expected 1-9" % sys.argv[1])

#plot types 1-5 are the plots of the entire Hz/j space
#plot type 6 plots a single point solution, drawing the x,y,z,magnetization components vs time
#plot type 7 makes an animated solution, drawing the x,y,z, components in a 3d plane and making a new drawing for every time step
#in the case of plot type 7, the second command line argument tells us how many steps to make.
#plot type 8 makes an x-y plot of h-sweep vs applied field
#plot type 9 makes an x-y plot of j-sweep vs applied field
if(int(sys.argv[1]) <= 5):			
	spinsim_make_plots(int(sys.argv[1]))
if(int(sys.argv[1]) == 6):
	if(len(sys.argv) > 2):
		spinsim_make_singleplot(int(sys.argv[1]), int(sys.argv[2]))
	else:
		spinsim_make_singleplot(int(sys.argv[1]), 2)
if(int(sys.argv[1]) == 7):
	spinsim_make_multiplot(int(sys.argv[1]), int(sys.argv[2]))
if(int(sys.argv[1]) == 8):
	spinsim_make_xyplot(int(sys.argv[1]))
#argv==9 makes it plot both l0 and l1 on a single single-plot axes set
if(int(sys.argv[1]) == 9):
	spinsim_make_xyplot(int(sys.argv[1]))
	
	
	
