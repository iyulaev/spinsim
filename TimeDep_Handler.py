#Stuff for PyQt
from PyQt4 import QtCore, QtGui, Qt
from timedep_dialog import Ui_Dialog

class TimeDep_Handler(QtGui.QMainWindow):
	def __init__(self, parent, hz_timdep_funs, j_timdep_funs):
		"""This method sets up the dialog for the time-dependent function editor"""
		QtGui.QWidget.__init__(self, parent)
		self.ui = Ui_Dialog()
		self.ui.setupUi(self)
		
		#save the parameters we were given
		self.hz_funs = hz_timdep_funs
		self.j_funs = j_timdep_funs
		self.parent_dialog = parent
		
		#intialize variables for sin dependent/independent, timedep/timeindep
		self.j_sin_dep = [0,0,0]
		self.j_sin_indep = [0,0,0]
		self.j_user_dep = ""
		self.j_user_indep = ""
		
		self.hz_sin_dep = [0,0,0]
		self.hz_sin_indep = [0,0,0]
		self.hz_user_dep = ""
		self.hz_user_indep = ""
		
		#load the rest of the parameters from the file
		self.loadFromFile()
		self.variablesToDialog()
		
		self.show()
		
	def variablesToDialog(self):
		"""This function sets up all of the GUI fields according to the stored
		variables. It is called by the __init__ function of TimeDep_Handler"""
		#Set up all of the j and hz text fields
		self.ui.j_dep_ampl.setText(str(self.j_sin_dep[0]))
		self.ui.j_dep_freq.setText(str(self.j_sin_dep[1]))
		self.ui.j_dep_phas.setText(str(self.j_sin_dep[2]))
		self.ui.j_indep_ampl.setText(str(self.j_sin_indep[0]))
		self.ui.j_indep_freq.setText(str(self.j_sin_indep[1]))
		self.ui.j_indep_phas.setText(str(self.j_sin_indep[2]))
		self.ui.j_userdep_file.setText(self.j_user_dep)
		self.ui.j_userindep_file.setText(self.j_user_indep)
		
		self.ui.hz_dep_ampl.setText(str(self.hz_sin_dep[0]))
		self.ui.hz_dep_freq.setText(str(self.hz_sin_dep[1]))
		self.ui.hz_dep_phas.setText(str(self.hz_sin_dep[2]))
		self.ui.hz_indep_ampl.setText(str(self.hz_sin_indep[0]))
		self.ui.hz_indep_freq.setText(str(self.hz_sin_indep[1]))
		self.ui.hz_indep_phas.setText(str(self.hz_sin_indep[2]))
		self.ui.hz_userdep_file.setText(self.hz_user_dep)
		self.ui.hz_userindep_file.setText(self.hz_user_indep)
		
		#DC current enabled?
		if((self.j_funs % 2) == 1):
			self.ui.j_dc_en.setCheckState(QtCore.Qt.Checked)
		else:
			self.ui.j_dc_en.setCheckState(QtCore.Qt.Unchecked)
		#amplitude-dependent, sine wave
		if(((self.j_funs>>1) % 2) == 1):
			self.ui.j_sin_dep_en.setCheckState(QtCore.Qt.Checked)
		else:
			self.ui.j_sin_dep_en.setCheckState(QtCore.Qt.Unchecked)
		#amplitude-independent, sine wave
		if(((self.j_funs>>2) % 2) == 1):
			self.ui.j_sin_indep_en.setCheckState(QtCore.Qt.Checked)
		else:
			self.ui.j_sin_indep_en.setCheckState(QtCore.Qt.Unchecked)
		#user amplitude-dependent function
		if(((self.j_funs>>3) % 2) == 1):
			self.ui.j_userdep_en.setCheckState(QtCore.Qt.Checked)
		else:
			self.ui.j_userdep_en.setCheckState(QtCore.Qt.Unchecked)
		#user amplitude-independent function
		if(((self.j_funs>>4) % 2) == 1):
			self.ui.j_userindep_en.setCheckState(QtCore.Qt.Checked)
		else:
			self.ui.j_userindep_en.setCheckState(QtCore.Qt.Unchecked)
			
		#DC magnetic field enabled?
		if((self.hz_funs % 2) == 1):
			self.ui.hz_dc_en.setCheckState(QtCore.Qt.Checked)
		else:
			self.ui.hz_dc_en.setCheckState(QtCore.Qt.Unchecked)
		#amplitude-dependent, sine wave
		if(((self.hz_funs>>1) % 2) == 1):
			self.ui.hz_sin_dep_en.setCheckState(QtCore.Qt.Checked)
		else:
			self.ui.hz_sin_dep_en.setCheckState(QtCore.Qt.Unchecked)
		#amplitude-independent, sine wave
		if(((self.hz_funs>>2) % 2) == 1):
			self.ui.hz_sin_indep_en.setCheckState(QtCore.Qt.Checked)
		else:
			self.ui.hz_sin_indep_en.setCheckState(QtCore.Qt.Unchecked)
		#user amplitude-dependent function
		if(((self.hz_funs>>3) % 2) == 1):
			self.ui.hz_userdep_en.setCheckState(QtCore.Qt.Checked)
		else:
			self.ui.hz_userdep_en.setCheckState(QtCore.Qt.Unchecked)
		#user amplitude-independent function
		if(((self.hz_funs>>4) % 2) == 1):
			self.ui.hz_userindep_en.setCheckState(QtCore.Qt.Checked)
		else:
			self.ui.hz_userindep_en.setCheckState(QtCore.Qt.Unchecked)
		
	def dialogToVariables(self):
		"""This function saves all of the GUI fields to the internal class variables"""
		#extract all of the textbox fields from GUI
		self.j_sin_dep[0] = float(self.ui.j_dep_ampl.text())
		self.j_sin_dep[1] = float(self.ui.j_dep_freq.text())
		self.j_sin_dep[2] = float(self.ui.j_dep_phas.text())
		self.j_sin_indep[0] = float(self.ui.j_indep_ampl.text())
		self.j_sin_indep[1] = float(self.ui.j_indep_freq.text())
		self.j_sin_indep[2] = float(self.ui.j_indep_phas.text())
		self.j_user_dep = self.ui.j_userdep_file.text()
		self.j_user_indep = self.ui.j_userindep_file.text()
		
		self.hz_sin_dep[0] = float(self.ui.hz_dep_ampl.text())
		self.hz_sin_dep[1] = float(self.ui.hz_dep_freq.text())
		self.hz_sin_dep[2] = float(self.ui.hz_dep_phas.text())
		self.hz_sin_indep[0] = float(self.ui.hz_indep_ampl.text())
		self.hz_sin_indep[1] = float(self.ui.hz_indep_freq.text())
		self.hz_sin_indep[2] = float(self.ui.hz_indep_phas.text())
		self.hz_user_dep = self.ui.hz_userdep_file.text()
		self.hz_user_indep = self.ui.hz_userindep_file.text()
		
		self.j_funs = 0		
		#extract checkboxes from GUI; set up internal variables
		if(self.ui.j_dc_en.checkState() ==  QtCore.Qt.Checked):
			self.j_funs += 1
		if(self.ui.j_sin_dep_en.checkState() ==  QtCore.Qt.Checked):
			self.j_funs += 2
		if(self.ui.j_sin_indep_en.checkState() ==  QtCore.Qt.Checked):
			self.j_funs += 4
		if(self.ui.j_userdep_en.checkState() ==  QtCore.Qt.Checked):
			self.j_funs += 8
		if(self.ui.j_userindep_en.checkState() ==  QtCore.Qt.Checked):
			self.j_funs += 16
			
		self.hz_funs = 0
		#do the same thing, this time for Hz rather than j
		if(self.ui.hz_dc_en.checkState() ==  QtCore.Qt.Checked):
			self.hz_funs += 1
		if(self.ui.hz_sin_dep_en.checkState() ==  QtCore.Qt.Checked):
			self.hz_funs += 2
		if(self.ui.hz_sin_indep_en.checkState() ==  QtCore.Qt.Checked):
			self.hz_funs += 4
		if(self.ui.hz_userdep_en.checkState() ==  QtCore.Qt.Checked):
			self.hz_funs += 8
		if(self.ui.hz_userindep_en.checkState() ==  QtCore.Qt.Checked):
			self.hz_funs += 16
		
		
	def loadFromFile(self):
		"""This function loads the fun_params file and stores what it read from the file to the internal
		class variables. Typically called right when the dialog loads so that we might start from the last
		edits we made to the fun_params file"""
		try:
			my_file = open("fun_params", 'r')
		except IOError:
			print "Couldn't open file! Sorry."
			Qt.QMessageBox.warning(self, "Info", "Couldn't open fun_params file! Maybe it doesn't exist?")
			return 0
		
		for line in my_file:
			#split the line on whitespace
			string_split = line.split()
			
			if(line[0] == '#'):
				print("Skipping comment")
			#match line like "j_sin_dep	0.2	1000000000	1.571"
			elif(string_split[0] == "j_sin_dep"):
				self.j_sin_dep[0] = float(string_split[1])
				self.j_sin_dep[1] = float(string_split[2])
				self.j_sin_dep[2] = float(string_split[3])
			elif(string_split[0] == "j_sin_indep"):
				self.j_sin_indep[0] = float(string_split[1])
				self.j_sin_indep[1] = float(string_split[2])
				self.j_sin_indep[2] = float(string_split[3])
			elif(string_split[0] == "hz_sin_dep"):
				self.hz_sin_dep[0] = float(string_split[1])
				self.hz_sin_dep[1] = float(string_split[2])
				self.hz_sin_dep[2] = float(string_split[3])
			elif(string_split[0] == "hz_sin_indep"):
				self.hz_sin_indep[0] = float(string_split[1])
				self.hz_sin_indep[1] = float(string_split[2])
				self.hz_sin_indep[2] = float(string_split[3])
			elif(string_split[0] == "j_user_dep"):
				self.j_user_dep = string_split[1]
			elif(string_split[0] == "j_user_indep"):
				self.j_user_indep = string_split[1]
			elif(string_split[0] == "hz_user_dep"):
				self.hz_user_dep = string_split[1]
			elif(string_split[0] == "hz_user_indep"):
				self.hz_user_indep = string_split[1]
			else:
				print("Couldn't parse line %s" % line)
		
		my_file.close()
		return 0
		
	def saveToFile(self):
		"""This function saves the parameters to file 'filename'
		Intended for saving your runs, to be used later with loadParamsFromFile"""
		
		#grab latest parameters from GUI fields
		self.dialogToVariables()
		
		try:
			my_file = open("fun_params", 'w')
		except IOError:
			print ("Couldn't open file %s! Sorry." % my_file)
			Qt.QMessageBox.warning(self, "Info", "Couldn't open fun_params file for writing! Sorry.")
			return 0
			
		my_file.write("j_sin_dep\t%.5e\t%.5e\t%.5e\n" % (self.j_sin_dep[0], self.j_sin_dep[1], self.j_sin_dep[2]))
		my_file.write("j_sin_indep\t%.5e\t%.5e\t%.5e\n" % (self.j_sin_indep[0], self.j_sin_indep[1], self.j_sin_indep[2]))
		my_file.write("hz_sin_dep\t%.5e\t%.5e\t%.5e\n" % (self.hz_sin_dep[0], self.hz_sin_dep[1], self.hz_sin_dep[2]))
		my_file.write("hz_sin_indep\t%.5e\t%.5e\t%.5e\n" % (self.hz_sin_indep[0], self.hz_sin_indep[1], self.hz_sin_indep[2]))
		my_file.write("j_user_dep\t%s\n" % self.j_user_dep)
		my_file.write("j_user_indep\t%s\n" % self.j_user_indep)
		my_file.write("hz_user_dep\t%s\n" % self.hz_user_dep)
		my_file.write("hz_user_indep\t%s\n" % self.hz_user_indep)
		
		my_file.close()
		return 0
		
	def accept(self):
		"""This method should save the GUI fields to internal variables, save
		internal variables to file, and send the j_funs and hz_funs back to the
		parent dialog class"""
		self.dialogToVariables()
		self.saveToFile()
		self.parent_dialog.timedep_writeback(self.j_funs, self.hz_funs)
		print("User hit accept!")
		self.hide()
	
	
	#Ivan: I think hiding the dialog and not "closing" it or destorying it in any way
	#is bad practice. This should probably be addressed at some point
	def reject(self):
		"""This methods simply hides the dialog"""
		print("User hit reject!")
		self.hide()

