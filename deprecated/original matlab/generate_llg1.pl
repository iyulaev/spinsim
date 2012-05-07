#!/usr/bin/perl

#This program generates llg1 based on the number of iterations given
#The output file "llg1_<m>.m" will have the Hz(k) parameter equal to (-1*magnitude) * (m/n)*(2*magnitude)
#Assumes the "source" file will be located in path in filename "llg1_src.m.src", else change line 14 (?)
#Written by Ivan Yulaev 2009-12-22, contact iyulaev@ucsd.edu

#Arguments will be 
#1: number of loop iterations, assumes number of inner and outer are the same (perhaps change this later)
#2: magnitude
#3: percent of K1 that K2 is
#4: thetaoffset

#Note: In the llg1_src.m.src file, the keyword "PINEAPPLE" will be replaced by the value of Hz(k) appropriate
#for that file. Thus, "PINEAPPLE" is a protected keyword.
#Also "APPLESEED" will be replaced by "n=<n>", so "APPLESEED" is also protected

my($source_name) = "<llg1_src.m.src";

#Check the number of arguments we got
if($#ARGV != 3) {
	my($temp) = $#ARGV + 1;
	print("Got wrong number of arguments. Expecting 4, got $temp\n");
	exit;
}

#set the number of iterations we will run
my($num_iterations) = $ARGV[0];
#set the magnitude we will use
my($magnitude) = $ARGV[1];
#set the percent of K1 that K2 is
my($perck2) = $ARGV[2];
#set the thetaoffset
my($tofset) = $ARGV[3];

#outer loop for generating files
my($count2) = 0;
while($count2 < $num_iterations) {
	#the loop that shits out $num_iterations files
	my($count) = 0;
	while($count < $num_iterations) {
		#open relevant file handles
		open(SRC_FILE, $source_name);
		open(DEST_FILE, (">llg1_" . ($count2+1) . "_" . ($count+1) . ".m"));
		
		#print("DEBUG: DEST file for iteration $count will be: ");
		#print (">llg1_" . $count . ".m");
		#print("\n");

		#go through the file and shit out an instance
		while(<SRC_FILE>) {
			my($line) = $_;
			chomp($line);
			
			#if we find a line that has keyword "PINEAPPLE" replace the line with "banana = #" (# is whatever we should have there, to match
			#the behavior of linspace(-$magnitude, $magnitude, $num_iterations) on the $count-th iteration)
			if($line =~ m/PINEAPPLE/) {
				#print("DEBUG: Found line pineapple.\n");
				my($curr_mag) = (-1.0 * $magnitude) + (2.0 * ($count/($num_iterations-1)) * $magnitude);
				#print("DEBUG: Will write current magnitude at $curr_mag for iteration $count\n");
				print DEST_FILE "\tbanana = $curr_mag;\n";
			}
			#replace the line with keyword APPLESEED with the value of n
			elsif($line =~ m/APPLESEED/) {
				print DEST_FILE "\tn = $num_iterations;\n";
			}
			#line with keyword IMPORT_THETAOFFSET will be replaced by given thetaoffset
			elsif($line =~ m/IMPORT_THETAOFFSET/) {
				print DEST_FILE "\tthetaoff = $tofset;\n";
			}
			#line with keyword IMPORT_PERCK2 will be replaced by given percent of k2
			elsif($line =~ m/IMPORT_PERCK2/) {
				print DEST_FILE "\tpercK2 = $perck2;\n";
			}
			#line with keyword IMPORT_EYE will be replaced by current outer loop iteration
			elsif($line =~ m/IMPORT_EYE/) {
				my($temp) = $count2 + 1;
				print DEST_FILE "\teye = $temp;\n";
			}
			#otherwise don't change the line!
			else {
				print DEST_FILE "$line\n";
			}
		}
		
		#be a good program and close relevant file handles
		close(DEST_FILE);
		close(SRC_FILE);
		$count = $count + 1;
	}
	$count2 = $count2 + 1;
}

#print("DEBUG: Done running $count iterations\n");
