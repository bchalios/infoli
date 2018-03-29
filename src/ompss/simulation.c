/*
 *
 * Copyright (c) 2012, Neurasmus B.V., The Netherlands,
 * web: www.neurasmus.com email: voltage@neurasmus.com
 *
 * Any use or reproduction in whole or in parts is prohibited
 * without the written consent of the copyright owner.
 *
 * All Rights Reserved.
 *
 *
 * Author: Sebastian Isaza
 * Created: 19-01-2012
 * Modified: 07-08-2012
 *
 * Description: Top source file of the Inferior Olive model, originally written
 * in Matlab by Jornt De Gruijl. It contains the implementation of all functions.
 * The main function allocates the necessary memory, initializes the system
 * state and runs the model calculations.
 *
 */

/*
 * we assume that dim1 refers to the number of ROWS
 * and dim2 refers to COLUMNS (cell network size).
 */


#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <time.h>
#include <omp.h>
#include "infoli.h"
#include "infoli_log.h"
#include "utils.h"
//#include <mic_power.h>

int main(int argc, char *argv[])
{
	infoli_conf_t infoli_config = { 0 };
	char outFileName[PATH_MAX] = { 0 };
	char paramsFileName[PATH_MAX] = { 0 };
	char exeLink[PATH_MAX] = { 0 };
	int total_simulation_steps;
	float simTime = 0;

	/* Process command line arguments
	 * If we have no input from file then a one-pulse input is stimulated.
	 * Otherwise we receive input from a specified file and simulation runs accordingly
	 * in the case of inputFromFile, we will also need a buffer to hold the stimulus
	 * for each cell on each step
	 */

	/* Default values for parameters */
	IO_NETWORK_SIZE = 1000;
	CONN_PROBABILITY = 0.5;
	simTime = 5000;


	sprintf(exeLink, "/proc/self/exe");
	if (readlink(exeLink, paramsFileName, PATH_MAX) == -1) {
		printf("Error: Couldn't read simlink to running executable.\n");
		exit(EXIT_FAILURE);
	}
	stopAtSubstring(paramsFileName, "infoli.x");
	strcat(paramsFileName, "default.conf");

	sprintf(outFileName,"InferiorOlive_Output.txt");

	/* Argument parsing
	*/ 

	int option = 0;
	while ((option = getopt(argc, argv,"hn:p:t:i:c:o:")) != -1) {
		switch (option) {
			case 'h' : print_usage();
				   return 0;
			case 'n' : IO_NETWORK_SIZE = atoi(optarg);
				   if (IO_NETWORK_SIZE<=0) {
					   printf("Incorrect network size argument. Using default value of 1000.\n");
					   IO_NETWORK_SIZE = 1000;
				   }
				   break;
			case 'p' : CONN_PROBABILITY = atof(optarg);
				   if (CONN_PROBABILITY<0) {
					   printf("Incorrect network density argument. Using default value of 0.5.\n");
					   CONN_PROBABILITY = 0.5;
				   }
				   if (CONN_PROBABILITY==0)
					   printf("Warning: Simulating non-connected network (network density = 0).\n");
				   if (CONN_PROBABILITY>1) {
					   printf("Warning: Setting network density to 1 (=100%).\n");
					   CONN_PROBABILITY = 1;
				   }
				   break;
			case 't' : simTime = atof(optarg);
				   if (simTime<=0) {
					   printf("Incorrect simulation time argument. Using default value of 5000ms.\n");
					   simTime = 5000;
				   }
				   break;
			case 'c' :
				   strncpy(paramsFileName, optarg, PATH_MAX);
				   break;
			case 'o' :
				   strncpy(outFileName, optarg, PATH_MAX);
				   break;
			default: print_usage();
				 exit(EXIT_FAILURE);
		}
	}

	/* Initialize infoli logging */
	infoli_init_log(outFileName);

	/* Initialize problem */
	infoli_init(paramsFileName, &infoli_config);


	/* Start the simulation */
	total_simulation_steps = (int)(simTime/DELTA);

	//mic_power_start(1000, 200);
	gettimeofday(&tic, NULL);

	simulate(total_simulation_steps, &infoli_config);

	gettimeofday(&toc, NULL);
	//double energy = mic_power_finish();
	//printf("Sire, we measure energy levels of %0.3lf kJ.\n", energy);

	/* SIM END: Free  memory and close files */

	/* Execution Time for the Sim */

	printf("Execution Time for Simulation: %0.2f ms.\n",
		((toc.tv_sec*1000.0 + ((float)(toc.tv_usec)/1000.0)) -
		(tic.tv_sec*1000.0 + ((float)(tic.tv_usec)/1000.0))));

#ifndef G_CAL_FROM_FILE
	infoli_log("Execution Time for Simulation in ms: %0.2f\n",
		((toc.tv_sec * 1000.0 + ((float)(toc.tv_usec) / 1000.0)) -
		 (tic.tv_sec * 1000.0 + ((float)(tic.tv_usec) / 1000.0))));
#endif

	infoli_close_log();

	return 0;
}
