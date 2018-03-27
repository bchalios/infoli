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
#include <time.h>
#include <omp.h>
#include "infoli.h"
#include "utils.h"
//#include <mic_power.h>

int main(int argc, char **argv){
	
	/* declaration of variables
	 * necessary for program flow
	 */

	int i=0, j, k=0, line_count, targetCore, print_granularity=1;
	char c;
	char *inFileName;
	char *outFileName;
	char *paramsFileName;
	char *exeLink;
	FILE *pInFile, *coreF;
	char conFile[200],core_file[100];
	FILE *pConFile,*pOutFile;
	mod_prec* iAppArray= NULL;
	int simStep = 0, inputFromFile, initSteps,total_simulation_steps;
	float simTime = 0;
	int seedvar;
	char tempbuf[100];		
	mod_prec iApp, voltage;
	mod_prec *gcal;
	int *init_steps;
	long double seconds;
	int simulation_array_ID,target_cell, requested_neighbour;
	mod_prec* initValues= NULL;

	cellCompParams cellParamsPtr;

	/* declaration of variables
	 * necessary for ion channel
	 * computations
	 */

	mod_prec q_inf, tau_q, dq_dt;
	mod_prec alpha_r, beta_r, r_inf, tau_r, dr_dt;
	mod_prec alpha_s, beta_s, s_inf, tau_s, ds_dt;
	mod_prec dCa_dt, f, V, I_c_storage;
	mod_prec I_sd, I_CaH_temp, I_K_Ca, I_ld, I_h, dVd_dt;

	mod_prec k_inf, l_inf, tau_k, tau_l, dk_dt, dl_dt;
	mod_prec m_inf, h_inf, tau_h, dh_dt;
	mod_prec n_inf, p_inf, tau_n, tau_p, dn_dt, dp_dt;
	mod_prec alpha_x_s, beta_x_s, x_inf_s, tau_x_s, dx_dt_s;
	mod_prec I_ds, I_CaL, I_Na_s, I_ls, I_Kdr_s, I_K_s, I_as, dVs_dt;

	mod_prec m_inf_a, h_inf_a, tau_h_a, dh_dt_a;
	mod_prec alpha_x_a, beta_x_a, x_inf_a, tau_x_a, dx_dt_a;
	mod_prec I_Na_a, I_la, I_sa, I_K_a, dVa_dt;

	/* end of declarations
	 */ 

	/* Process command line arguments
	 * If we have no input from file then a one-pulse input is stimulated.
	 * Otherwise we receive input from a specified file and simulation runs accordingly
	 * in the case of inputFromFile, we will also need a buffer to hold the stimulus
	 * for each cell on each step
  	 */

	/* We first set some default values
 	*/ 	

	IO_NETWORK_SIZE = 1000;
	CONN_PROBABILITY = 0.5;
	simTime = 5000;

	inFileName = (char*) malloc(PATH_MAX*sizeof(char));
	inputFromFile = 0;

	paramsFileName = (char*) malloc(PATH_MAX*sizeof(char));
	exeLink = (char*) malloc(PATH_MAX*sizeof(char));
	sprintf(exeLink, "/proc/self/exe");
	if (readlink(exeLink, paramsFileName, PATH_MAX) == -1) {
		printf("Error: Couldn't read simlink to running executable.\n");
		exit(EXIT_FAILURE);
	}
	stopAtSubstring(paramsFileName, "infoli.x");
	strcat(paramsFileName, "default.conf");

	outFileName = (char*) malloc(PATH_MAX*sizeof(char));
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
			case 'i' : inputFromFile = 1; inFileName = optarg;
				break;
			case 'c' : paramsFileName = optarg;
				break;
			case 'o' : outFileName = optarg;
				break;
			default: print_usage();
				exit(EXIT_FAILURE);
		}
	}
	
	if (inputFromFile) {
		pInFile = fopen(inFileName,"r");
		if(pInFile==NULL) {
			printf("Error: Couldn't open %s\n", inFileName);
			exit(EXIT_FAILURE);
		}
		iAppArray= (mod_prec*) _mm_malloc(cellCount*(sizeof(mod_prec)), 64);
	}	

        /* compute how many grid cells 
	 * are assigned to each core
	*/

        cellCount= IO_NETWORK_SIZE;

	/* PRINTING is true in case we want to write the
	 * output to a specified file (outFileName)
	 */

	if (PRINTING) {
		pOutFile = fopen(outFileName,"w+");
		if(pOutFile==NULL){
			printf("Error: Couldn't create %s\n", outFileName);
			exit(EXIT_FAILURE);
		}
		print_granularity=100;
		#ifdef G_CAL_FROM_FILE
			print_granularity=1;
		#endif
	}

	/* Channel declaration and mem allocation 
	 * 
	 */

	int* cellID = (int*) _mm_malloc(cellCount*sizeof(int), 64);
	mod_prec* V_dend = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Hcurrent_q = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Calcium_r = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Potassium_s = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* I_CaH = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Ca2Plus = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* iAppIn = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* I_c = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* V_soma = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* g_CaL = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Sodium_m = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Sodium_h = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Calcium_k = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Calcium_l = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Potassium_n = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Potassium_p = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Potassium_x_s = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* V_axon = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Sodium_m_a = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Sodium_h_a = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Potassium_x_a = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);

	if (!(Potassium_x_a)) {
		printf("Error: Couldn't allocate memory for channels\n");
		exit(EXIT_FAILURE);
	}

	/* Initialise simulation with appropriate values.
	*/
	initValues= (mod_prec*) malloc(42*sizeof(mod_prec));
	read_parameters(paramsFileName, initValues);

	/* Generic Simulation Values
	 */ 

	DELTA = initValues[0];
	CONDUCTANCE = initValues[1];
	C_M = initValues[2];
	G_NA_S = initValues[3];
	G_KDR_S = initValues[4];
	G_K_S = initValues[5];
	G_LS = initValues[6];
	G_K_CA = initValues[7];
	G_CAH = initValues[8];
	G_LD = initValues[9];
	G_H = initValues[10];
	G_NA_A = initValues[11];
	G_NA_R = initValues[12];
	G_K_A = initValues[13];
	G_LA = initValues[14];
	P1 = initValues[15];
	P2 = initValues[16];
	G_INT = initValues[17];
	V_NA = initValues[18];
	V_K = initValues[19];
	V_CA = initValues[20];
	V_H = initValues[21];
	V_L = initValues[22];

	/*  Network Specific Values
	 */

	for (i=0;i<cellCount;i++) {
		cellID[i] = i;
		//Initial dendritic parameters
		V_dend[i] = initValues[23];
		Calcium_r[i] = initValues[24];// High-threshold calcium
		Potassium_s[i] = initValues[25];// Calcium-dependent potassium
		Hcurrent_q[i] = initValues[26];// H current
		Ca2Plus[i] = initValues[27];// Calcium concentration
		I_CaH[i] = initValues[28];// High-threshold calcium current
		//Initial somatic parameters
		g_CaL[i] = initValues[29]; //default arbitrary value but it should be randomized per cell
		V_soma[i] = initValues[30];
		Sodium_m[i] = initValues[31];// Sodium (artificial)
		Sodium_h[i] = initValues[32];
		Potassium_n[i] = initValues[33];// Potassium (delayed rectifier)
		Potassium_p[i] = initValues[34];
		Potassium_x_s[i] = initValues[35];// Potassium (voltage-dependent)
		Calcium_k[i] = initValues[36];// Low-threshold calcium
		Calcium_l[i] = initValues[37];
		// Initial axonal parameters
		V_axon[i] = initValues[38];
		//sisaza: Sodium_m_a doesn't have a state, therefore this assignment doesn'thave any effect
		Sodium_m_a[i] = initValues[39];// Sodium (thalamocortical)
		Sodium_h_a[i] = initValues[40];
		Potassium_x_a[i] = initValues[41];// Potassium (transient)
	}

	#ifdef G_CAL_FROM_FILE
		read_g_CaL_from_file(g_CaL);
	#endif

	/*Value Initialization Complete
	 */ 

	/* cellCompParams contains
	 * connection information for the NW
	 */


	//1D-array storing how many connections each neuron has
	cellParamsPtr.total_amount_of_neighbours = (int*) _mm_malloc(cellCount*sizeof(int), 64);
	for (i=0;i<cellCount;i++)
		cellParamsPtr.total_amount_of_neighbours[i] = 0;	//needed bugfix initialization
	//2D-array storing the ids of each connection
	cellParamsPtr.neighId = (int**) _mm_malloc(cellCount*sizeof(int*), 64);
	//2D-array storing the conductances of each connection
	cellParamsPtr.neighConductances = (mod_prec**) _mm_malloc(cellCount*sizeof(mod_prec*), 64);

	if(cellParamsPtr.neighConductances==NULL){
		printf("Error: Couldn't malloc for cellParamsPtr\n");
		exit(EXIT_FAILURE);
	}

	mod_prec cond_value= CONDUCTANCE;	//default value for conductance in synapses (microSiemens)

	/* we shall generate a connectivity map
	 * based on the probability given by the user
	 */

	int sender_cell, receiver_cell;
	int* conn_gen_buffer=(int*) calloc(sizeof(int), IO_NETWORK_SIZE);	//temp buffer for storing bonds
	float rndm;

	/* we generate rng checks against the probability
	 * supplied by the user to calculate connection counts
 	 * We store the results of each rng check in a temp buffer
 	 * After we are done, we allocate memory for dendritic data structs
 	 * and fill it with neighbour (bonds) ids
 	 */

	srand ( time(NULL) );
	for (receiver_cell=0; receiver_cell<IO_NETWORK_SIZE; receiver_cell++) {

//		#pragma omp parallel for shared (cellParamsPtr, conn_gen_buffer, receiver_cell, IO_NETWORK_SIZE, CONN_PROBABILITY) private(sender_cell, rndm)
		for (sender_cell=0;sender_cell<IO_NETWORK_SIZE;sender_cell++) {

			if (sender_cell == receiver_cell)
				;		//no self-feeding connections allowed
			else {
				rndm = ((float) rand()) / ((float) RAND_MAX);	//generate rng and compare to probability
				if (rndm <= CONN_PROBABILITY) {
					cellParamsPtr.total_amount_of_neighbours[receiver_cell]++;  //increase neighbour count
					conn_gen_buffer[sender_cell]++;	//mark that we formed a bond with this cell
				}
			}
		}

		//allocate enough space now that we know how many neighbours this receiving cell has
		cellParamsPtr.neighConductances[receiver_cell] = (mod_prec*) _mm_malloc(cellParamsPtr.total_amount_of_neighbours[receiver_cell]*sizeof(mod_prec), 64);
		cellParamsPtr.neighId[receiver_cell] = (int*) _mm_malloc(cellParamsPtr.total_amount_of_neighbours[receiver_cell]*sizeof(int), 64);

		//run through the temporary buffer and fill the data structs with the bonds' info
		i=0;	//this temporary variable will help fill the data structs
//		#pragma omp parallel for shared (cellParamsPtr, conn_gen_buffer, receiver_cell, IO_NETWORK_SIZE, cond_value) private(sender_cell) firstprivate(i)
		for (sender_cell=0;sender_cell<IO_NETWORK_SIZE;sender_cell++) {
			if (i>cellParamsPtr.total_amount_of_neighbours[receiver_cell])
				;
			else if (conn_gen_buffer[sender_cell]==0)
				;	//skip this cell, it is not a bond
			else {
				cellParamsPtr.neighConductances[receiver_cell][i]=cond_value;
				cellParamsPtr.neighId[receiver_cell][i]=sender_cell;
				i++;
			}
		}

		//reset the buffer for the next receiver cell
		memset(conn_gen_buffer, 0, IO_NETWORK_SIZE*sizeof(int));
	}

	/* connections mapping
 	 * is now complete
 	 */


	//Initialize output file IF enabled
	#ifndef G_CAL_FROM_FILE
		if (PRINTING) {
			sprintf(tempbuf, "Execution Time for Simulation in ms:              \n");
				fputs(tempbuf, pOutFile);
			sprintf(tempbuf, "#simSteps Input(Iapp) Output(V_axon)");
				fputs(tempbuf, pOutFile);
			for (i=0;i<cellCount;i++) {
				sprintf(tempbuf, "%d ", cellID[i]);
					fputs(tempbuf, pOutFile);
			}
			sprintf(tempbuf, "\n");
			fputs(tempbuf, pOutFile);
		}
	#endif

	//Initialize g_CaL
	seedvar = time(NULL);		//seedvar will be time- and core-dependent now!
	srand(seedvar);

	if (RAND_INIT)
		for(i=0;i<cellCount;i++)
			g_CaL[i] = 0.6+(0.2*(rand()%100)/100);

	//random initialization process
	if (RAND_INIT) {
		for(i=0;i<cellCount;i++) {
			initSteps = rand()%(int)ceil(100/DELTA);
			initSteps = initSteps | 0x00000001;//make initialization steps odd

			for(j=0;j<initSteps;j++){
				//Arrange input
				iAppIn[i] = 0;//No stimulus
				
//				CompDend(&cellParamsPtr[i], 1);  PLACEHOLDER RANDOMIZATION, NEEDS RE-CODING
//				CompSoma(&cellParamsPtr[i]);
//				CompAxon(&cellParamsPtr[i]);
			}

		}
	}

	/* start of the simulation
	 * In case we want to read the stimulus from file inputFromFile = true
	 */
//	mic_power_start(1000, 200);
	gettimeofday(&tic, NULL);

	/* 	WARNING, recent changes have made inputFromFile currently unusable-
 	 *  	related code will be commented out until fixed and commited
 	 */

/*	if(inputFromFile){
		simStep = 0;


		while(ReadFileLine(pInFile, iAppArray)) {
			
			simulation_array_ID = simStep%2;
			if (PRINTING) {
				sprintf(tempbuf, "%d %.2f ", (simStep+1), simStep*0.05); // start @ 1 because skipping initial values
				fputs(tempbuf, pOutFile);
			}
			
		#pragma omp parallel for shared(cellParamsPtr, cellPtr) private(iAppArray, target_cell, tempbuf, i, requested_neighbour) firstprivate(simulation_array_ID)	
			for (target_cell=0;target_cell<cellCount;target_cell++) {
				
				for (i=0; i<cellParamsPtr[target_cell].total_amount_of_neighbours; i++){
                                        requested_neighbour = cellParamsPtr[target_cell].neighId[i];
                                        cellParamsPtr[target_cell].neighVdend[i] = cellPtr[simulation_array_ID][requested_neighbour].dend.V_dend;
				}

				cellParamsPtr[target_cell].iAppIn = iAppArray[target_cell];
				cellParamsPtr[target_cell].prevCellState = &cellPtr[simulation_array_ID][target_cell];
				cellParamsPtr[target_cell].newCellState = &cellPtr[simulation_array_ID^1][target_cell];

				CompDend(&cellParamsPtr[target_cell], 0);
				CompSoma(&cellParamsPtr[target_cell]);
				CompAxon(&cellParamsPtr[target_cell]);

				if (PRINTING) {
					sprintf(tempbuf, "%.16f ", cellPtr[(simulation_array_ID)^1][target_cell].axon.V_axon);
					fputs(tempbuf, pOutFile);
				}

			}
			if (PRINTING) {
				sprintf(tempbuf, "\n");
				fputs(tempbuf,pOutFile);
			}

			simStep++;
		}
		
	}
	else {	
*/		total_simulation_steps = (int)(simTime/DELTA);

		for(simStep=0;simStep<total_simulation_steps;simStep++) {
								
			if ((simStep>=20000)&&(simStep<20500-1))
				iApp = 6; 
			else
		   		iApp = 0;
			#ifndef G_CAL_FROM_FILE
				if (PRINTING&&((simStep%print_granularity)==0)) {
					sprintf(tempbuf, " %d %.2f ", simStep+1, iApp);
					fputs(tempbuf, pOutFile);
				}
			#endif

	/*	First Loop takes care of Gap Junction Functions
 	*/			

		#pragma omp parallel for shared (cellParamsPtr, iApp, iAppIn, V_dend, I_c, pOutFile) private(target_cell, i, \
		requested_neighbour, f, V, voltage, I_c_storage) firstprivate(simulation_array_ID)
			for (target_cell=0;target_cell<cellCount;target_cell++) {

				//Feeding of Input Current
				iAppIn[target_cell] = iApp;

	/*	Gathering of the data concerning Vdend of neighours, then
 	*	computing and storing the incoming Ic from neighbours
 	*/

				I_c_storage = 0;
				__assume_aligned(cellParamsPtr.neighId[target_cell], 64);
				__assume_aligned(cellParamsPtr.neighConductances[target_cell], 64);
				__assume_aligned(V_dend, 64);
			#pragma ivdep
				for (i=0; i<cellParamsPtr.total_amount_of_neighbours[target_cell]; i++){
					requested_neighbour = cellParamsPtr.neighId[target_cell][i];
					voltage = V_dend[requested_neighbour];
					V = V_dend[target_cell] - voltage;
					f = 0.8f * expf(-1*powf(V, 2)/100) + 0.2f;// SCHWEIGHOFER 2004 VERSION
					I_c_storage += cellParamsPtr.neighConductances[target_cell][i] * f * V;
                                }
				I_c[target_cell] = I_c_storage;
			}

	/* Second Loop computes the new state (voltages, channels etc.)
 	*/ 

			__assume_aligned(V_dend, 64);
			__assume_aligned(Hcurrent_q, 64);
			__assume_aligned(Calcium_r, 64);
			__assume_aligned(Potassium_s, 64);
			__assume_aligned(I_CaH, 64);
			__assume_aligned(Ca2Plus, 64);
			__assume_aligned(I_c, 64);
			__assume_aligned(iAppIn, 64);

			__assume_aligned(V_soma, 64);
			__assume_aligned(g_CaL, 64);
			__assume_aligned(Sodium_m, 64);
			__assume_aligned(Sodium_h, 64);
			__assume_aligned(Calcium_k, 64);
			__assume_aligned(Calcium_l, 64);
			__assume_aligned(Potassium_n, 64);
			__assume_aligned(Potassium_p, 64);
			__assume_aligned(Potassium_x_s, 64);

			__assume_aligned(V_axon, 64);
			__assume_aligned(Sodium_m_a, 64);
			__assume_aligned(Sodium_h_a, 64);
			__assume_aligned(Potassium_x_a, 64);
		#pragma omp parallel for simd shared (cellParamsPtr, V_dend, Hcurrent_q, Calcium_r, Potassium_s, I_CaH, Ca2Plus,\
		iApp, iAppIn, I_c, V_soma, g_CaL, Sodium_m, Sodium_h, Calcium_k, Calcium_l, Potassium_n, Potassium_p, Potassium_x_s,\
		V_axon, Sodium_m_a, Sodium_h_a, Potassium_x_a, pOutFile) private(target_cell, tempbuf, q_inf, tau_q, dq_dt,\
		alpha_r, beta_r, r_inf, tau_r, dr_dt, alpha_s, beta_s, s_inf, tau_s, ds_dt, dCa_dt, I_sd, I_CaH_temp, I_K_Ca, I_ld,\
		I_h, dVd_dt, k_inf, l_inf, tau_k, tau_l, dk_dt, dl_dt, m_inf, h_inf, tau_h, dh_dt, n_inf, p_inf, tau_n, tau_p, dn_dt,\
		dp_dt, alpha_x_s, beta_x_s, x_inf_s, tau_x_s, dx_dt_s, I_ds, I_CaL, I_Na_s, I_ls, I_Kdr_s, I_K_s, I_as, dVs_dt,\
		m_inf_a, h_inf_a, tau_h_a, dh_dt_a, alpha_x_a, beta_x_a, x_inf_a, tau_x_a, dx_dt_a, I_Na_a, I_la, I_sa, I_K_a,\
		dVa_dt)	firstprivate(simulation_array_ID)
			for (target_cell=0;target_cell<cellCount;target_cell++) {

//				~DENDRITIC COMPUTATIONS~

//				Dend H Current Calcs
				q_inf = 1 /(1 + expf((V_dend[target_cell] + 80) / 4));
			        tau_q = 1 /(expf(-0.086f * V_dend[target_cell] - 14.6f) + expf(0.070f * V_dend[target_cell] - 1.87f));
			        dq_dt = (q_inf - Hcurrent_q[target_cell]) / tau_q;
				Hcurrent_q[target_cell] = DELTA * dq_dt + Hcurrent_q[target_cell];

//				Dend Ca Current Calcs
				alpha_r = 1.7f / (1 + expf( -(V_dend[target_cell] - 5) / 13.9f));
		       		beta_r = 0.02f * (V_dend[target_cell] + 8.5f) / (expf((V_dend[target_cell] + 8.5f) / 5) - 1);
	  		     	r_inf = alpha_r / (alpha_r + beta_r);
	       	 		tau_r = 5 / (alpha_r + beta_r);
	       		 	dr_dt = (r_inf - Calcium_r[target_cell]) / tau_r;
	    			Calcium_r[target_cell] = DELTA * dr_dt + Calcium_r[target_cell];

//				Dend K Current Calcs
				alpha_s = min((0.00002f * Ca2Plus[target_cell]), 0.01f);
			        beta_s = 0.015f;
			        s_inf = alpha_s / (alpha_s + beta_s);
			        tau_s = 1 / (alpha_s + beta_s);
			        ds_dt = (s_inf - Potassium_s[target_cell]) / tau_s;
			        Potassium_s[target_cell] = DELTA * ds_dt + Potassium_s[target_cell];

//				Dend Cal Current Calcs
				dCa_dt = -3 * I_CaH[target_cell] - 0.075f * Ca2Plus[target_cell];
				Ca2Plus[target_cell] = DELTA * dCa_dt + Ca2Plus[target_cell];

//				Dendritic Voltage And Current Calcs
				I_sd = (G_INT / (1 - P1)) * (V_dend[target_cell] - V_soma[target_cell]);
				I_CaH_temp =G_CAH * powf(Calcium_r[target_cell], 2) * (V_dend[target_cell] - V_CA);
				I_K_Ca =G_K_CA * Potassium_s[target_cell] * (V_dend[target_cell] - V_K);
				I_ld =G_LD * (V_dend[target_cell] - V_L);
				I_h=G_H * Hcurrent_q[target_cell] * (V_dend[target_cell] - V_H);
				dVd_dt = (-(I_CaH_temp + I_sd+ I_ld + I_K_Ca + I_c[target_cell] + I_h) + iAppIn[target_cell]) / C_M;
				I_CaH[target_cell] = I_CaH_temp;


//				~SOMATIC COMPUTATIONS~

//				Soma Calcium Calcs
				k_inf = (1 / (1 + expf(-1 * (V_soma[target_cell] + 61) / 4.2f)));
			        l_inf = (1 / (1 + expf(( V_soma[target_cell] + 85.5f) / 8.5f)));
			        tau_k = 1;
			        tau_l = ((20 * expf((V_soma[target_cell] + 160) / 30) / (1 + expf((V_soma[target_cell] + 84) / 7.3f))) +35);
			        dk_dt = (k_inf - Calcium_k[target_cell]) / tau_k;
			        dl_dt = (l_inf - Calcium_l[target_cell]) / tau_l;
			        Calcium_k[target_cell] = DELTA * dk_dt + Calcium_k[target_cell];
			        Calcium_l[target_cell] = DELTA * dl_dt + Calcium_l[target_cell];

//				Soma Sodium Calcs
				m_inf = 1 / (1 + (expf((-30 - V_soma[target_cell])/ 5.5f)));
				h_inf = 1 / (1 + (expf((-70 - V_soma[target_cell])/-5.8f)));
				tau_h = 3 * expf((-40 - V_soma[target_cell])/33);
				dh_dt = (h_inf - Sodium_h[target_cell])/tau_h;
				Sodium_m[target_cell] = m_inf;
				Sodium_h[target_cell] = Sodium_h[target_cell] + DELTA * dh_dt;

//				Soma Potassium Calcs
				n_inf = 1 / (1 + expf( ( -3 - V_soma[target_cell]) /10));
			        p_inf = 1 / (1 + expf( (-51 - V_soma[target_cell]) / -12));
			        tau_n = 5 + (47 * expf( -(-50 - V_soma[target_cell]) /900));
			        tau_p = tau_n;
			        dn_dt = (n_inf - Potassium_n[target_cell]) / tau_n;
			        dp_dt = (p_inf - Potassium_p[target_cell]) / tau_p;
			        Potassium_n[target_cell] = DELTA * dn_dt + Potassium_n[target_cell];
				Potassium_p[target_cell] = DELTA * dp_dt + Potassium_p[target_cell];

//				Soma Potassium X Calcs
				alpha_x_s = 0.13f * (V_soma[target_cell] + 25) / (1 - expf(-(V_soma[target_cell] + 25) / 10));
			        beta_x_s= 1.69f * expf(-0.0125f * (V_soma[target_cell] + 35));
				x_inf_s = alpha_x_s / (alpha_x_s + beta_x_s);
				tau_x_s = 1 / (alpha_x_s + beta_x_s);
			        dx_dt_s = (x_inf_s - Potassium_x_s[target_cell]) / tau_x_s;
			        Potassium_x_s[target_cell] = 0.05f * dx_dt_s + Potassium_x_s[target_cell];

//				Somatic Voltage And Current Calcs
				I_ds= (G_INT / P1) * (V_soma[target_cell] - V_dend[target_cell]);
				I_CaL = g_CaL[target_cell] * powf(Calcium_k[target_cell], 3) * Calcium_l[target_cell] * (V_soma[target_cell] - V_CA);
				I_Na_s= G_NA_S * powf(Sodium_m[target_cell], 3) * Sodium_h[target_cell] * (V_soma[target_cell] - V_NA);
				I_ls= G_LS * (V_soma[target_cell] - V_L);
				I_Kdr_s = G_KDR_S * powf(Potassium_n[target_cell], 4) * (V_soma[target_cell] - V_K);
				I_K_s = G_K_S * powf(Potassium_x_s[target_cell], 4) * (V_soma[target_cell] - V_K);
				I_as= (G_INT / (1 - P2)) * (V_soma[target_cell] - V_axon[target_cell]);
				dVs_dt = (-(I_CaL + I_ds+ I_as + I_Na_s + I_ls + I_Kdr_s + I_K_s)) / C_M;


//				~AXONAL COMPUTATIONS~

//				Axon Sodium Calcs
				m_inf_a = 1 / (1 + (expf((-30 - V_axon[target_cell])/ 5.5f)));
				h_inf_a = 1 / (1 + (expf((-60 - V_axon[target_cell])/-5.8f)));
				tau_h_a = 1.5f * expf((-40 - V_axon[target_cell])/33);
				dh_dt_a = (h_inf_a - Sodium_h_a[target_cell])/tau_h_a;
				Sodium_m_a[target_cell] = m_inf_a;
				Sodium_h_a[target_cell] = Sodium_h_a[target_cell] + DELTA * dh_dt_a;

//				Axon Potassium Calcs
				alpha_x_a = 0.13f * (V_axon[target_cell] + 25) / (1 - expf(-(V_axon[target_cell] + 25) / 10));
			        beta_x_a= 1.69f * expf(-0.0125f * (V_axon[target_cell] + 35));
			        x_inf_a = alpha_x_a / (alpha_x_a + beta_x_a);
			        tau_x_a = 1 / (alpha_x_a + beta_x_a);
			        dx_dt_a = (x_inf_a - Potassium_x_a[target_cell]) / tau_x_a;
			        Potassium_x_a[target_cell] = 0.05f * dx_dt_a + Potassium_x_a[target_cell];

//				Axonal Voltage And Current Calcs
				I_Na_a= G_NA_A * powf(Sodium_m_a[target_cell], 3) * Sodium_h_a[target_cell] * (V_axon[target_cell] - V_NA);
				I_la= G_LA* (V_axon[target_cell] - V_L);
				I_sa= (G_INT / P2) * (V_axon[target_cell] - V_soma[target_cell]);
				I_K_a = G_K_A * powf(Potassium_x_a[target_cell], 4) * (V_axon[target_cell] - V_K);
				dVa_dt = (-(I_K_a + I_sa + I_la + I_Na_a)) / C_M;


//				~NEW VOLTAGES~

				V_dend[target_cell] = DELTA * dVd_dt + V_dend[target_cell];
				V_soma[target_cell] = DELTA * dVs_dt + V_soma[target_cell];
				V_axon[target_cell] = DELTA * dVa_dt + V_axon[target_cell];

//				~END OF COMPUTATIONS~
 
			}


			if (PRINTING&&((simStep%print_granularity)==0)) {

			#pragma omp parallel for simd shared (V_axon, pOutFile) private(target_cell, tempbuf)
				for (target_cell=0;target_cell<cellCount;target_cell++) {
					#ifndef G_CAL_FROM_FILE
						sprintf(tempbuf, "%d : %.8f ", target_cell+1, V_axon[target_cell]);
					#else
						sprintf(tempbuf, "%.8f ", V_axon[target_cell]);
					#endif
					fputs(tempbuf, pOutFile);
				}

				sprintf(tempbuf, "\n");
				fputs(tempbuf, pOutFile);
			}
		
		}

//	}

	gettimeofday(&toc, NULL);
//	double energy = mic_power_finish();
//	printf("Sire, we measure energy levels of %0.3lf kJ.\n", energy);
	
	/* SIM END
	 * Free  memory and close files
	 */

	if (PRINTING) {
		fclose (pOutFile);
		chmod(outFileName,0x01B6);
	}

	if(inputFromFile)
		fclose (pInFile);

	/* Execution Time for the Sim
	 */

	printf("Execution Time for Simulation: %0.2f ms.\n", ((toc.tv_sec*1000.0 + ((float)(toc.tv_usec)/1000.0)) - \
		(tic.tv_sec*1000.0 + ((float)(tic.tv_usec)/1000.0))));
	if (PRINTING) {
		pOutFile = fopen (outFileName, "rb+");
		sprintf(tempbuf, "Execution Time for Simulation in ms: %0.2f\n", \
			((toc.tv_sec*1000.0 + ((float)(toc.tv_usec)/1000.0)) - \
			(tic.tv_sec*1000.0 + ((float)(tic.tv_usec)/1000.0))));
		#ifndef G_CAL_FROM_FILE
			fputs(tempbuf, pOutFile);
		#endif
		fclose (pOutFile);
	}

	return 0;
}
