/* 	timedepfuns.c
	This file implements the time-dependent functions for current and applied field for Spinsim
	It allows the use of sinusoidal and user-defined functions for Hz and j variables, used in the llg
	differential equation functions.
	Written and Engineered 2010-01-23 I. Yulaev (iyulaev@ucsd.edu)
	University of California San Diego, ECE Department, Fullerton Lab
	See the manual for information on how to use this module
*/

//Defines for first parameter of TimeDepFun
#define CALC_FIELD 0
#define CALC_CURRENT 1
//Size of the line buffer in LoadTDParams() for reading strings from a file
#define LINEBUF_SIZE 256
//What file to read the parameters from
#define TIME_DEP_PARAM_FILE "fun_params"
//Bitwise values for parsing the functions to use
#define DC_FUN 1
#define SIN_DEP_FUN 2
#define SIN_INDEP_FUN 4
#define USER_DEP_FUN 8
#define USER_INDEP_FUN 16
//Print out debug values if defined as 1
//#define DEBUG 1

#include <string.h>

//Amplitude (scalar), frequency, phase for DC-dependent sine function for j
double j_sin_dep_params[3];
//Amplitude (absolute, in A/m**2), frequency, phase for DC-independent sine function for j
double j_sin_indep_params[3];
//List of time-value (scalar) pairs for user-defined DC-dependent function for j
double * j_user_dep_params;
//List of time-value (absolute) pairs for user-defined DC-independent function for j
double * j_user_indep_params;

//Amplitude (scalar), frequency, phase for ext-field-dependent sine function for Hz
double hz_sin_dep_params[3];
//Amplitude (absolute, in T), frequency, phase for ext-field-independent sine function for Hz
double hz_sin_indep_params[3];
//List of time-value (scalar) pairs for user-defined ext-field-dependent function for Hz
double * hz_user_dep_params;
//List of time-value (absolute, in T) pairs for user-defined ext-field-independent function for Hz
double * hz_user_indep_params;

/*LoadTDUserParams() loads in a time-dependent parameters file from user input
User file format is "<time in seconds>\t<either current in A/m**2 or field in T>"
Params:
(1) User file path string is given as parameter
(2) Also given whether we're doing field or current
(3) Also given what function we're filling in (dependent or independent)
Returns 0 on success, else failure */
int LoadTDUserParams(char * file_name, int j_or_hz, int what_function) {
	//Variable for counting # of lines (to allocate the right amount of space)
	int number_of_lines = 0;
	//Variable for storing which variable we wrote to
	double* debug_data_read;
	//Temporary line buffer
	char linebuf[LINEBUF_SIZE];
	
	//First pass: Count the number of lines in the file
	//Open the parameter files
	FILE* tdp_user_file;
	tdp_user_file = fopen(file_name, "r");
//	Check that we could open the file
	if(tdp_user_file == NULL) {
		printf("ERROR: LoadTDUserParams() couldn't open user file %s\n", file_name);
		return(1);
	}
	while(fgets(linebuf, LINEBUF_SIZE, tdp_user_file)) {
		if(strlen(linebuf) > 1 && linebuf[0] != '#') {
			number_of_lines++;
		}
	}
	
	close(tdp_user_file);
	
	//Allocate enough memory, we need 2 doubles for each "line" in our user function file because it's (time, amplitude) pairs
	if(j_or_hz == CALC_FIELD && what_function == USER_DEP_FUN) {
		hz_user_dep_params = (double*) malloc(2 * sizeof(double) * (number_of_lines+1));
		hz_user_dep_params[0] = (double) number_of_lines;
	}
	else if(j_or_hz == CALC_FIELD && what_function == USER_INDEP_FUN) {
		hz_user_indep_params = (double*) malloc(2 * sizeof(double) * (number_of_lines+1));
		hz_user_indep_params[0] = (double) number_of_lines;
	}
	else if(j_or_hz == CALC_CURRENT && what_function == USER_DEP_FUN) {
		j_user_dep_params = (double*) malloc(2 * sizeof(double) * (number_of_lines+1));
		j_user_dep_params[0] = (double) number_of_lines;
	}
	else if(j_or_hz == CALC_CURRENT && what_function == USER_INDEP_FUN) {
		j_user_indep_params = (double*) malloc(2 * sizeof(double) * (number_of_lines+1));
		j_user_indep_params[0] = (double) number_of_lines;
	}
	else {
		printf("ERROR: LoadTDUserParams given unrecognized parameters, couldn't load user file %s.\n", file_name);
		return(1);
	}
	
	//Reload the file and go through it again, this time putting everything into the variable
	int count = 2;
	tdp_user_file = fopen(file_name, "r");
//	Check that we could open the file
	if(tdp_user_file == NULL) {
		printf("ERROR: LoadTDUserParams() couldn't open user file %s\n", file_name);
		return(1);
	}
	
	//Go through the file, again, and this time save the values to the appropriate global variable
	//These variables are stored and will be later called by llg1() functions to calculate the h field and
	//current tensity at time t using TimeDepFun()
	while(fgets(linebuf, LINEBUF_SIZE, tdp_user_file)) {
		if(strlen(linebuf) > 1 && linebuf[0] != '#') {
			char * curr_token;
			
			if(j_or_hz == CALC_FIELD && what_function == USER_DEP_FUN) {
				curr_token = strtok(linebuf, "\t\n");
				hz_user_dep_params[count++] = atof(curr_token);
				curr_token = strtok(NULL, "\t\n");
				hz_user_dep_params[count++] = atof(curr_token);
				debug_data_read = hz_user_dep_params;
			}
			else if(j_or_hz == CALC_FIELD && what_function == USER_INDEP_FUN) {
				curr_token = strtok(linebuf, "\t\n");
				hz_user_indep_params[count++] = atof(curr_token);
				curr_token = strtok(NULL, "\t\n");
				hz_user_indep_params[count++] = atof(curr_token);
				debug_data_read = hz_user_indep_params;
			}
			else if(j_or_hz == CALC_CURRENT && what_function == USER_DEP_FUN) {
				curr_token = strtok(linebuf, "\t\n");
				j_user_dep_params[count++] = atof(curr_token);
				curr_token = strtok(NULL, "\t\n");
				j_user_dep_params[count++] = atof(curr_token);
				debug_data_read = j_user_dep_params;
			}
			else if(j_or_hz == CALC_CURRENT && what_function == USER_INDEP_FUN) {
				curr_token = strtok(linebuf, "\t\n");
				j_user_indep_params[count++] = atof(curr_token);
				curr_token = strtok(NULL, "\t\n");
				j_user_indep_params[count++] = atof(curr_token);
				debug_data_read = j_user_indep_params;
			}
		}
	}
	
	close(tdp_user_file);
	
	#ifdef DEBUG
	printf("DEBUG: For file %s, read in data:\n", file_name);
	for(count = 0; count < number_of_lines+1; count++) {
		printf("%e\t%e\n", debug_data_read[2*count], debug_data_read[2*count+1]);
	}
	#endif
	
	return(0);	
}

/*LoadTDParams() loads the time dependent parameters from a time-dependent parameter file (fun_params)
For sin functions, it simply stores the amplitude/frequency/phase information in the global variable
For user-defined functions, it calls LoadTDUserParams() to load the user table into memory
For format of a time-dependent parameter file, please see documentation 
Returns 0 on success, else failure*/
int LoadTDParams() {
	//line buffer
	char linebuf[LINEBUF_SIZE];
	
	//Open the parameter files
	FILE* tdp_file;
	tdp_file = fopen(TIME_DEP_PARAM_FILE, "r");
	
	if(tdp_file == NULL) {
		printf("ERROR: Couldn't open time dependent parameter file. \n");
		return(1);
	}
	
//	Read in the parameter file one line at a time, and parse each line
	while(fgets(linebuf, LINEBUF_SIZE, tdp_file)) {
//		Extract first token from the line
		char temp_buf[256];
		char temp_char = linebuf[0];
		int count = 0;
		
		if(strlen(linebuf) > 1 && linebuf[0] != '#') {
			#ifdef DEBUG
			printf("DEBUG: Loaded line: %s", linebuf);
			#endif
			
			while(temp_char != '\t') {
				temp_buf[count++] = temp_char;
				temp_char = linebuf[count];
			}
			temp_buf[count] = '\0';
		
	//		Try to match line with "j_sin_dep"
			if(strcmp(temp_buf, "j_sin_dep") == 0 && strlen(temp_buf) != 0 && temp_buf[0] != '#') {
	//			We will tokenize the string. Throw out the first one, though
				char* curr_token;
				int count2 = 0;
				curr_token = strtok(linebuf, "\t\n");
				curr_token = strtok(NULL, "\t\n");
				//TODO: This loop is not safe!
				//Store the read in values to global variable
				while(curr_token != NULL) {
					j_sin_dep_params[count2++] = atof(curr_token);
					curr_token = strtok(NULL, "\t\n");
				}
			}
	//		Try to match line with "j_sin_indep"
			else if(strcmp(temp_buf, "j_sin_indep") == 0 && strlen(temp_buf) != 0 && temp_buf[0] != '#') {
	//			We will tokenize the string. Throw out the first one, though
				char* curr_token;
				int count2 = 0;
				curr_token = strtok(linebuf, "\t\n");
				curr_token = strtok(NULL, "\t\n");
				while(curr_token != NULL) {
					j_sin_indep_params[count2++] = atof(curr_token);
					curr_token = strtok(NULL, "\t\n");
				}
			}
	//		Try to match line with "hz_sin_dep"
			else if(strcmp(temp_buf, "hz_sin_dep") == 0 && strlen(temp_buf) != 0 && temp_buf[0] != '#') {
	//			We will tokenize the string. Throw out the first one, though
				char* curr_token;
				int count2 = 0;
				curr_token = strtok(linebuf, "\t\n");
				curr_token = strtok(NULL, "\t\n");
				while(curr_token != NULL) {
					hz_sin_dep_params[count2++] = atof(curr_token);
					curr_token = strtok(NULL, "\t\n");
				}
			}
	//		Try to match line with "hz_sin_indep"
			else if(strcmp(temp_buf, "hz_sin_indep") == 0 && strlen(temp_buf) != 0 && temp_buf[0] != '#') {
	//			We will tokenize the string. Throw out the first one, though
				char* curr_token;
				int count2 = 0;
				curr_token = strtok(linebuf, "\t\n");
				curr_token = strtok(NULL, "\t\n");
				while(curr_token != NULL) {
					hz_sin_indep_params[count2++] = atof(curr_token);
					curr_token = strtok(NULL, "\t\n");
				}
			}
	//		Try to match line with "hz_user_dep"
			else if(strcmp(temp_buf, "hz_user_dep") == 0 && strlen(temp_buf) != 0 && temp_buf[0] != '#') {
	//			We will tokenize the string. Throw out the first one, though
				char* curr_token;
				int count2 = 0;
				curr_token = strtok(linebuf, "\t\n");
				curr_token = strtok(NULL, "\t\n");
				
				#ifdef DEBUG
				printf("DEBUG: LoadTDParams about to try to load user dependent hz function from file %s\n", curr_token);
				#endif
				LoadTDUserParams(curr_token, CALC_FIELD, USER_DEP_FUN);
			}
	//		Try to match line with "hz_user_indep"
			else if(strcmp(temp_buf, "hz_user_indep") == 0 && strlen(temp_buf) != 0 && temp_buf[0] != '#') {
	//			We will tokenize the string. Throw out the first one, though
				char* curr_token;
				int count2 = 0;
				curr_token = strtok(linebuf, "\t\n");
				curr_token = strtok(NULL, "\t\n");
				
				#ifdef DEBUG
				printf("DEBUG: LoadTDParams about to try to load user independent hz function from file %s\n", curr_token);
				#endif
				LoadTDUserParams(curr_token, CALC_FIELD, USER_INDEP_FUN);
			}
	//		Try to match line with "j_user_dep"
			else if(strcmp(temp_buf, "j_user_dep") == 0 && strlen(temp_buf) != 0 && temp_buf[0] != '#') {
	//			We will tokenize the string. Throw out the first one, though
				char* curr_token;
				int count2 = 0;
				curr_token = strtok(linebuf, "\t\n");
				curr_token = strtok(NULL, "\t\n");
				
				#ifdef DEBUG
				printf("DEBUG: LoadTDParams about to try to load user dependent j function from file %s\n", curr_token);
				#endif
				LoadTDUserParams(curr_token, CALC_CURRENT, USER_DEP_FUN);
			}
	//		Try to match line with "j_user_indep"
			else if(strcmp(temp_buf, "j_user_indep") == 0 && strlen(temp_buf) != 0 && temp_buf[0] != '#') {
	//			We will tokenize the string. Throw out the first one, though
				char* curr_token;
				int count2 = 0;
				curr_token = strtok(linebuf, "\t\n");
				curr_token = strtok(NULL, "\t\n");
				
				#ifdef DEBUG
				printf("DEBUG: LoadTDParams about to try to load user independent j function from file %s\n", curr_token);
				#endif
				LoadTDUserParams(curr_token, CALC_CURRENT, USER_INDEP_FUN);
			}
			else {
				printf("ERROR: Unrecognized string parsing file, %s\n", temp_buf);
			}
		} //[/if strlen > 1]
	}
	
	close(tdp_file);
	
	#ifdef DEBUG
	printf("DEBUG: Got the following values from the file:\n");
	printf("j_sin_dep = %e\t%e\t%e\n", j_sin_dep_params[0], j_sin_dep_params[1], j_sin_dep_params[2]);
	printf("j_sin_indep = %e\t%e\t%e\n", j_sin_indep_params[0], j_sin_indep_params[1], j_sin_indep_params[2]);
	printf("hz_sin_dep = %e\t%e\t%e\n", hz_sin_dep_params[0], hz_sin_dep_params[1], hz_sin_dep_params[2]);
	printf("hz_sin_indep = %e\t%e\t%e\n", hz_sin_indep_params[0], hz_sin_indep_params[1], hz_sin_indep_params[2]);
	#endif
	
	return(0);
}	

/*TimeDepFun() calculates either the time-dependent field or current (depending on the value of the first parameter)
The function is given the above parameter, as well as time current time t and current j/Hz value, and a bit-wise or
of the functions to use, and it returns the result.
Uses external variables "sin_dep_params", "sin_indep_params", "user_indep_vals", "user_dep_vals" as the parameters
for the time-dependent functions
Params:
val_to_cal - Number telling us what function we are calcuating, current (CALC_CURRENT) or field (CALC_FIELD)
t - the current time
j_Hz - the j or Hz value that llg1() is using for its current run
what_funs - the bitwise-ORd value that tells us what functions to include in the result. See line 17 "bitwise values",
	what_funs is a bitwise-ORd number of these values.
Returns: time dependent field or current */
double TimeDepFun(unsigned int val_to_cal, double t, double j_Hz, unsigned int what_funs, double t_relax) {
	double output_value = 0.0;
	
	//Added for "relaxation time" - i.e. don't apply neither field nor current for the first t_relax seconds
	if(t < t_relax) return(0.0);
	
	//Superimpose time-independent value, if necessary
	if((what_funs & DC_FUN) != 0) {
		output_value += j_Hz;
	}
	//Superimpose sine, j/Hz dependent function, if necessary
	if((what_funs & SIN_DEP_FUN) != 0) {
		//Calculating field or current?
		if(val_to_cal == CALC_FIELD) {
			output_value += j_Hz * hz_sin_dep_params[0] * sin(2 * PI * t * hz_sin_dep_params[1] + hz_sin_dep_params[2]);
		}
		else if(val_to_cal == CALC_CURRENT) {
			output_value += j_Hz * j_sin_dep_params[0] * sin(2 * PI * t * j_sin_dep_params[1] + j_sin_dep_params[2]);
		}
		else {
			printf("ERROR: TimeDepFun got invaid parameters. Results likely invalid!\n");
			exit(-1);
		}
	}
	//Superimpose sine, j/Hz independent function, if necessary
	if((what_funs & SIN_INDEP_FUN) != 0) {
		//Calculating field or current?
		if(val_to_cal == CALC_FIELD) {
			output_value += hz_sin_indep_params[0] * sin(2 * PI * t * hz_sin_indep_params[1] + hz_sin_indep_params[2]);
		}
		else if(val_to_cal == CALC_CURRENT) {
			output_value += j_sin_indep_params[0] * sin(2 * PI * t * j_sin_indep_params[1] + j_sin_indep_params[2]);
		}
		else {
			printf("ERROR: TimeDepFun got invaid parameters. Results likely invalid!\n");
			exit(-1);
		}
	}
	//Superimpose user-defined j/Hz function, if necessary
	if(((what_funs & USER_DEP_FUN) != 0) || ((what_funs & USER_INDEP_FUN) != 0)) {
		//Calculating field or current? Dependent or independent?
		double * data_values = NULL;
		double value_to_add;
		if((what_funs & USER_DEP_FUN) != 0) {
			data_values = (val_to_cal == CALC_FIELD) ? hz_user_dep_params : ((val_to_cal == CALC_CURRENT) ? j_user_dep_params : NULL);
		}
		else if((what_funs & USER_INDEP_FUN) != 0) {
			data_values = (val_to_cal == CALC_FIELD) ? hz_user_indep_params : ((val_to_cal == CALC_CURRENT) ? j_user_indep_params : NULL);
		}
		
		if(data_values == NULL) {
			printf("ERROR: TimeDepFun got invaid parameters. Results likely invalid!\n");
			exit(-1);
		}
		
		//Find the closest time. If we can't find an exact match, get two nearest neighbors and do linear approximation
		double closest_times[2] = {0.0, 0.0};
		double closest_values[2] = {0.0, 0.0};
		int count;
		int structure_length = (int) data_values[0];
		
		//find the time either closes to or equal to the time t
		for(count = 0; count < structure_length; count++) {
			//offset reads by 2 because 0th element is structure lenght, 1st element is reserved
			//we find the first time index that is GREATER THAN t
			if(data_values[(2*(count+1))] > t) {
				//if we're past the first index, then use the previous index as the time index
				//less than or equal to t
				if(count > 0) {
					closest_times[0] = data_values[(2*count)];
					closest_values[0] = data_values[(2*count) + 1];
				}
				//otherwise, mark the less_than_or_equal_to time as -1, to denote that the very
				//first time index in the data series was already greater than current time t
				else {
					closest_times[0] = -1.0;
					closest_values[0] = -1.0;
				}
				//Use the current times as the greater than or equal to time values
				closest_times[1] = data_values[(2*(count+1))];
				closest_values[1] = data_values[(2*(count+1)) + 1];
				//Mark us as having "not" gone off the end of the structure
				count = structure_length + 1;
			}
		}
		
		#ifdef DEBUG
		printf("DEBUG: At time t=%e ",t);
		#endif
		
		//if closest_times[0] is set to -1, then the first time index was already greater than t. 
		//In this case, just use the first data point in the set
		if(closest_times[0] == -1.0) {
			value_to_add = closest_values[1];
			#ifdef DEBUG
			printf("decided to use value %e, since that was the first value (with time t=%e)\n", closest_values[1], closest_times[1]);
			#endif
		}
		//if count == structure length, then we went to the very end without finding a large enough time index
		//In this case, use the last data point in the set
		else if(count == structure_length) {
			value_to_add = data_values[(2*structure_length)+1];
			#ifdef DEBUG
			printf("decided to use value %e, since that was the last value (with time t=%e)\n", data_values[(2*structure_length)+1], data_values[(2*structure_length)]);
			#endif
		}
		//if neither of the above apply, and the time at closest_times[0] is exactly t, use that value
		else if(closest_times[0] == t) {
			value_to_add = closest_values[0];
			#ifdef DEBUG
			printf("decided to use value %e, since that was the exact time index (with time t=%e)\n", closest_values[0], closest_times[0]);
			#endif
		}
		//if the above don't apply, and t is between two time indices, use linear approximation between the two
		else if(closest_times[0] < t && closest_times[1] > t) {
			double time_difference = (t - closest_times[0]) / (closest_times[1] - closest_times[0]);
			double val_difference = closest_values[1] - closest_values[0];
			value_to_add = closest_values[0] + (val_difference * time_difference);
			#ifdef DEBUG
			printf("decided to do linear approximation %e between points (%e, %e) and (%e, %e)\n", 
				value_to_add, closest_times[0], closest_values[0], closest_times[1], closest_values[1]);
			#endif
		}
		else {
			printf("ERROR: TimeDepFun() did something horribly wrong. Exiting...\n");
			exit(-1);
		}
		
		//Finally, superimpose the result on the output
		if((what_funs & USER_INDEP_FUN) != 0) {
			output_value += value_to_add;
		}
		else if((what_funs & USER_DEP_FUN) != 0) {
			output_value += value_to_add * j_Hz;
		}
		else {
			printf("ERROR: TimeDepFun() did something terribly wrong. Exiting...\n");
			exit(-1);
		}
	}	//[/Superimpose user-defined j/Hz function, if necessary]
	
	return(output_value);
}
