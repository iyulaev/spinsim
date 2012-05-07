# the default command is used the compile the C spinsim program
default: 
	gcc -pthread -lgsl -lgslcblas -lm -o spinsim spinsim.c

#pendantic does the same as default, but yells at you more
pedantic:
	gcc -W -Wall -ansi -pedantic -std=c99 -pthread -lgsl -lgslcblas -lm -o spinsim spinsim.c
	
# GUI compiles all of the python GUI files
gui:
	pyuic4 qt_ui_v2.ui > qt_ui_v2.py
	pyuic4 timedep_dialog.ui > timedep_dialog.py
	
#Compiles default, gui, and apngasm
all:
	gcc -pthread -lgsl -lgslcblas -lm -o spinsim spinsim.c
	pyuic4 qt_ui_v2.ui > qt_ui_v2.py
	pyuic4 timedep_dialog.ui > timedep_dialog.py
	gcc -O2 -o apngasm apngasm.c -lz -lpng 

#removes compiled files
clean:
	-rm *.pyc
	-rm spinsim
	
#removes compiled GUI files
cleangui:
	-rm *.pyc
	-rm qt_ui_v2.py
	-rm timedep_dialog.py

#removes everything that cleangui does and also apngasm
cleanall:
	-rm *.pyc
	-rm spinsim
	-rm qt_ui_v2.py
	-rm timedep_dialog.py
	-rm apngasm

#remove compiled files AND generated data runs
#be careful with this one!
cleandata:
	-rm *.pyc
	-rm spinsim
	-rm magsweep*
	-rm magsingle*
	
