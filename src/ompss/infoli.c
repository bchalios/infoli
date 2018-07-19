#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "utils.h"
#include "infoli.h"
#include "infoli_log.h"
#include "infoli_perf.h"

static inline void *_malloc64(size_t size)
{
	void *ret = _mm_malloc(size, 64);
	if (!ret) {
		perror("_mm_malloc");
		exit(1);
	}

	return ret;
}

static inline void *_malloc(size_t size)
{
	void *ret = malloc(size);
	if (!ret) {
		perror("malloc");
		exit(1);
	}

	return ret;
}

static inline void *_calloc(size_t nmemb, size_t size)
{
	void *ret = calloc(nmemb, size);
	if (!ret) {
		perror("calloc");
		exit(1);
	}

	return ret;
}

void infoli_init(char *params_file_name, infoli_conf_t *config)
{
	int i;

	assert(config != NULL);
	assert(params_file_name != NULL);

	printf("Initializing problem\n");

	config->cellCount = IO_NETWORK_SIZE;

	config->cellID = _malloc64(config->cellCount*sizeof(int));
	config->iAppIn =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->V_dend =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->V_axon =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->V_soma =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->Hcurrent_q =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->Calcium_r =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->Calcium_k =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->Calcium_l =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->Ca2Plus =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->Potassium_s =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->Potassium_n =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->Potassium_p =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->Potassium_x_s =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->Potassium_x_a =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->Sodium_m =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->Sodium_h =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->Sodium_m_a =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->Sodium_h_a =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->g_CaL =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->I_CaH =  _malloc64(config->cellCount*sizeof(mod_prec));
	config->I_c =  _malloc64(config->cellCount*sizeof(mod_prec));


	mod_prec *initValues= _malloc(42 * sizeof(mod_prec));
	read_parameters(params_file_name, initValues);


	/* Generic Simulation Values */ 
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


	/*  Network Specific Values */
	int cell_count = IO_NETWORK_SIZE;
	for (i = 0; i < cell_count; i++) {
		config->cellID[i] = i;
		//Initial dendritic parameters
		config->V_dend[i] = initValues[23];
		config->Calcium_r[i] = initValues[24];// High-threshold calcium
		config->Potassium_s[i] = initValues[25];// Calcium-dependent potassium
		config->Hcurrent_q[i] = initValues[26];// H current
		config->Ca2Plus[i] = initValues[27];// Calcium concentration
		config->I_CaH[i] = initValues[28];// High-threshold calcium current
		//Initial somatic parameters
		config->g_CaL[i] = initValues[29]; //default arbitrary value but it should be randomized per cell
		config->V_soma[i] = initValues[30];
		config->Sodium_m[i] = initValues[31];// Sodium (artificial)
		config->Sodium_h[i] = initValues[32];
		config->Potassium_n[i] = initValues[33];// Potassium (delayed rectifier)
		config->Potassium_p[i] = initValues[34];
		config->Potassium_x_s[i] = initValues[35];// Potassium (voltage-dependent)
		config->Calcium_k[i] = initValues[36];// Low-threshold calcium
		config->Calcium_l[i] = initValues[37];

		/* Initial axonal parameters */
		config->V_axon[i] = initValues[38];

		/* sisaza: Sodium_m_a doesn't have a state,
		 * therefore this assignment doesn'thave any effect */
		config->Sodium_m_a[i] = initValues[39];// Sodium (thalamocortical)
		config->Sodium_h_a[i] = initValues[40];
		config->Potassium_x_a[i] = initValues[41];// Potassium (transient)
	}


#ifdef G_CAL_FROM_FILE
	read_g_CaL_from_file(g_CaL);
#endif

	/* 1D-array storing how many connections each neuron has */

	cellCompParams *params = &config->cell_params;
	params->total_amount_of_neighbours =
		_malloc64(config->cellCount * sizeof(int));
	for (i = 0; i < config->cellCount; i++) {
		//needed bugfix initialization
		params->total_amount_of_neighbours[i] = 0;
	}
	
	//2D-array storing the ids of each connection
	params->neighId = _malloc64(config->cellCount*sizeof(int*));
	//2D-array storing the conductances of each connection
	params->neighConductances = _malloc64(config->cellCount*sizeof(mod_prec*));


	int sender_cell, receiver_cell;
	float rndm;
	
	/* temp buffer for storing bonds */
	int* conn_gen_buffer = _calloc(sizeof(int), IO_NETWORK_SIZE);

	srand(42);
	for (receiver_cell=0; receiver_cell < IO_NETWORK_SIZE; receiver_cell++) {

		/* #pragma omp parallel for \
		 * shared (cellParamsPtr, conn_gen_buffer, receiver_cell, \
	 	 *	IO_NETWORK_SIZE, CONN_PROBABILITY) \
		 * private(sender_cell, rndm) */
		for (sender_cell = 0;
			sender_cell < IO_NETWORK_SIZE;
			sender_cell++) {

			/* no self-feeding connections allowed */
			if (sender_cell == receiver_cell)
				continue;

			/*generate rng and compare to probability */
			rndm = ((float) rand()) / ((float) RAND_MAX);
			if (rndm <= CONN_PROBABILITY) {
				params->total_amount_of_neighbours[receiver_cell]++;  
				/* mark that we formed a bond with this cell */
				conn_gen_buffer[sender_cell]++;
			}
		}

		/* allocate enough space now that we know how many
		 * neighbours this receiving cell has */
		params->neighConductances[receiver_cell] = 
			_malloc64(params->total_amount_of_neighbours[receiver_cell]
					* sizeof(mod_prec));
		params->neighId[receiver_cell] =
			_malloc(params->total_amount_of_neighbours[receiver_cell]
					* sizeof(int));


		/* run through the temporary buffer and fill the data 
		 * structs with the bonds' info */
		i=0;

		/* #pragma omp parallel for \
		 * shared(cellParamsPtr, conn_gen_buffer, receiver_cell, \
		 * 	IO_NETWORK_SIZE, CONDUCTANCE) \
		 * private(sender_cell) firstprivate(i) */
		for (sender_cell = 0; sender_cell < IO_NETWORK_SIZE; sender_cell++) {
			if (i > params->total_amount_of_neighbours[receiver_cell])
				continue;
			
			/* skip this cell, it is not a bond */
			if (conn_gen_buffer[sender_cell] == 0)
				continue;
			
			params->neighConductances[receiver_cell][i] = CONDUCTANCE;
			params->neighId[receiver_cell][i] = sender_cell;
			i++;
		}

		//reset the buffer for the next receiver cell
		memset(conn_gen_buffer, 0, IO_NETWORK_SIZE*sizeof(int));
	}


#ifndef G_CAL_FROM_FILE
	/* Initialize output file */
	infoli_log("#simSteps Input(Iapp) Output(V_axon)\n");
	for (i = 0; i < config->cellCount; ++i)
		infoli_log("%d%c", config->cellID[i],
			(i == config->cellCount - 1) ? '\n' : ' ');
#endif

	/* Initialize g_CaL, with time and core-dependent value */
	int seedvar = time(NULL);
	srand(seedvar);

#ifdef RAND_INIT
	for(i = 0; i < config->cellCount; i++)
		config->g_CaL[i] = 0.6 + (0.2 * (rand() % 100) / 100);

	for(i = 0; i < config->cellCount; i++) {
		int initSteps = rand() % (int) ceil(100 / DELTA);
		/* make initialization steps odd */
		initSteps = initSteps | 0x00000001;

		int j;
		for(j = 0; j < initSteps; j++){
			/* No stimulus */
			config->iAppIn[i] = 0;
		}
	}
#endif /* RAND_INIT */
}


void gap_junction_functions(mod_prec iApp, int cellCount, mod_prec *iAppIn,
		mod_prec *V_dend, int **neighId, mod_prec **neighConductances,
		int *total_amount_of_neighbours, mod_prec *I_c)
{
	int target_cell, i;
	#pragma omp parallel for \
		shared(iAppIn, iApp, V_dend, I_c, neighId, \
			total_amount_of_neighbours, \
			neighConductances) \
		private(target_cell) \
		firstprivate(cellCount)
	for (target_cell = 0; target_cell < cellCount; ++target_cell) {
		/* Feeding of input current */
		iAppIn[target_cell] = iApp;

		/* Gather data concerning Vdend of neighbours, then
		 * comput and store the incoming Ic for neighbours */
		mod_prec I_c_storage = 0;
		__assume_aligned(neighId[target_cell], 64);
		__assume_aligned(neighConductances[target_cell], 64);
		__assume_aligned(V_dend, 64);

		int neighbours = total_amount_of_neighbours[target_cell];

		for (int i = 0; i < neighbours; ++i) {
			int requested_neighbour = neighId[target_cell][i];
			mod_prec voltage = V_dend[requested_neighbour];
			mod_prec V = V_dend[target_cell] - voltage;
			/* SCHWEIGHOFER 2004 VERSION */
			mod_prec f = 0.8f * expf(-1 * powf(V, 2) / 100) + 0.2f;
			I_c_storage += neighConductances[target_cell][i] * f * V;
		}

		I_c[target_cell] = I_c_storage;
	}
}

void update_state(int cellCount, mod_prec *V_dend, mod_prec *Hcurrent_q,
		mod_prec *Calcium_r, mod_prec *Ca2Plus, mod_prec *Potassium_s,
		mod_prec *I_CaH, mod_prec *V_soma, mod_prec *I_c,
		mod_prec *iAppIn, mod_prec *Calcium_k, mod_prec *Calcium_l,
		mod_prec *Sodium_h, mod_prec *Sodium_m, mod_prec *Potassium_n,
		mod_prec *Potassium_p, mod_prec *Potassium_x_s, mod_prec *g_CaL,
		mod_prec *Sodium_m_a, mod_prec *Sodium_h_a, mod_prec *V_axon,
		mod_prec *Potassium_x_a)
{
	int target_cell;


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

	for (target_cell = 0; target_cell < cellCount; ++target_cell) {
		/* ~DENDRITIC COMPUTATIONS~ */

		/* Dend H current calcs */
		mod_prec q_inf = 1 / (1 + expf((V_dend[target_cell] + 80) / 4));
		mod_prec tau_q = 1 / (expf(-0.086f * V_dend[target_cell] - 14.6f)
					+ expf(0.07f * V_dend[target_cell]
					- 1.87f));
		mod_prec dq_dt = (q_inf - Hcurrent_q[target_cell]) / tau_q;
		Hcurrent_q[target_cell] = DELTA * dq_dt
					+ Hcurrent_q[target_cell];
	
		/* Dend Ca current calcs */
		mod_prec alpha_r = 1.7f / (1 + expf(-(V_dend[target_cell] - 5)
					/ 13.9f));
		mod_prec beta_r = 0.02f * (V_dend[target_cell] + 8.5f)
				/ (expf((V_dend[target_cell] + 8.5f) / 5) - 1);
		mod_prec r_inf = alpha_r / (alpha_r + beta_r);
		mod_prec tau_r = 5 / (alpha_r + beta_r);
		mod_prec dr_dt = (r_inf - Calcium_r[target_cell]) / tau_r;
		Calcium_r[target_cell] = DELTA * dr_dt + Calcium_r[target_cell];

		/* Dend K current calcs */
		mod_prec alpha_s = min((0.00002f * Ca2Plus[target_cell]), 0.01f);
		mod_prec beta_s = 0.015f;
		mod_prec s_inf = alpha_s / (alpha_s + beta_s);
		mod_prec tau_s = 1 / (alpha_s + beta_s);
		mod_prec ds_dt = (s_inf - Potassium_s[target_cell]) / tau_s;
		Potassium_s[target_cell] = DELTA * ds_dt
					+ Potassium_s[target_cell];

		/* Dend Cal current calcs */
		mod_prec dCa_dt = -3 * I_CaH[target_cell] - 0.075f
				* Ca2Plus[target_cell];
		Ca2Plus[target_cell] = DELTA * dCa_dt + Ca2Plus[target_cell];

		/* Dendritic voltage and current calcs */
		mod_prec I_sd = (G_INT / (1 - P1)) * 
				(V_dend[target_cell] - V_soma[target_cell]);
		mod_prec I_CaH_temp = G_CAH * powf(Calcium_r[target_cell], 2)
					* (V_dend[target_cell] - V_CA);
		mod_prec I_K_Ca = G_K_CA * Potassium_s[target_cell] * (V_dend[target_cell] - V_K);
		mod_prec I_ld = G_LD * (V_dend[target_cell] - V_L);
		mod_prec I_h = G_H * Hcurrent_q[target_cell]
				* (V_dend[target_cell] - V_H);
		mod_prec dVd_dt = (-(I_CaH_temp + I_sd
					+ I_ld + I_K_Ca
					+ I_c[target_cell]
					+ I_h)
				+ iAppIn[target_cell]) / C_M;
		I_CaH[target_cell] = I_CaH_temp;

		/* ~SOMATIC COMPUTATIONS~ */

		/* Soma calcium calcs */
		mod_prec k_inf = (1 / (1 + expf(-1 * (V_soma[target_cell] + 61)
						/ 4.2f)));
		mod_prec l_inf = (1 / (1 + expf((V_soma[target_cell] + 85.5f)
						/ 8.5f)));
		mod_prec tau_k = 1;
		mod_prec tau_l = ((20 * expf((V_soma[target_cell] + 160) / 30)
					/ (1 + expf((V_soma[target_cell] + 84) 
							/ 7.3f))) + 35);
		mod_prec dk_dt = (k_inf - Calcium_k[target_cell]) / tau_k;
		mod_prec dl_dt = (l_inf - Calcium_l[target_cell]) / tau_l;
		Calcium_k[target_cell] = DELTA * dk_dt + Calcium_k[target_cell];
		Calcium_l[target_cell] = DELTA * dl_dt + Calcium_l[target_cell];

		/* Soma sodium calcs */
		mod_prec m_inf = 1 / (1 + (expf((-30 - V_soma[target_cell])
						/ 5.5f)));
		mod_prec h_inf = 1 / (1 + (expf((-70 - V_soma[target_cell])
						/-5.8f)));
		mod_prec tau_h = 3 * expf((-40 - V_soma[target_cell]) / 33);
		mod_prec dh_dt = (h_inf - Sodium_h[target_cell]) / tau_h;
		Sodium_m[target_cell] = m_inf;
		Sodium_h[target_cell] = Sodium_h[target_cell] + DELTA * dh_dt;

		/* Soma potassium calcs */
		mod_prec n_inf = 1 /
			(1 + expf(( -3 - V_soma[target_cell]) / 10));
		mod_prec p_inf = 1 /
			(1 + expf((-51 - V_soma[target_cell]) / -12));
		mod_prec tau_n = 5 +
			(47 * expf( -(-50 - V_soma[target_cell]) / 900));
		mod_prec tau_p = tau_n;
		mod_prec dn_dt = (n_inf - Potassium_n[target_cell]) / tau_n;
		mod_prec dp_dt = (p_inf - Potassium_p[target_cell]) / tau_p;
		Potassium_n[target_cell] = DELTA * dn_dt
			+ Potassium_n[target_cell];
		Potassium_p[target_cell] = DELTA * dp_dt
			+ Potassium_p[target_cell];

		/* Soma potassium X calcs */
		mod_prec alpha_x_s = 0.13f * (V_soma[target_cell] + 25)
				/ (1 - expf(-(V_soma[target_cell] + 25) / 10));
		mod_prec beta_x_s = 1.69f * expf(-0.0125f
				* (V_soma[target_cell] + 35));
		mod_prec x_inf_s = alpha_x_s / (alpha_x_s + beta_x_s);
		mod_prec tau_x_s = 1 / (alpha_x_s + beta_x_s);
		mod_prec dx_dt_s = (x_inf_s - Potassium_x_s[target_cell])
				/ tau_x_s;
		Potassium_x_s[target_cell] = 0.05f * dx_dt_s
				+ Potassium_x_s[target_cell];

		/* somatic voltage and current calcs */
		mod_prec I_ds = (G_INT / P1)
				* (V_soma[target_cell] - V_dend[target_cell]);
		mod_prec I_CaL = g_CaL[target_cell]
				* powf(Calcium_k[target_cell], 3)
				* Calcium_l[target_cell]
				* (V_soma[target_cell] - V_CA);
		mod_prec I_Na_s = G_NA_S * powf(Sodium_m[target_cell], 3)
				* Sodium_h[target_cell] * (V_soma[target_cell] - V_NA);
		mod_prec I_ls = G_LS * (V_soma[target_cell] - V_L);
		mod_prec I_Kdr_s = G_KDR_S * powf(Potassium_n[target_cell], 4)
				* (V_soma[target_cell] - V_K);
		mod_prec I_K_s = G_K_S * powf(Potassium_x_s[target_cell], 4)
				* (V_soma[target_cell] - V_K);
		mod_prec I_as = (G_INT / (1 - P2))
				* (V_soma[target_cell] - V_axon[target_cell]);
		mod_prec dVs_dt = (-(I_CaL + I_ds + I_as + I_Na_s + I_ls
					+ I_Kdr_s + I_K_s)) / C_M;

		/* ~AXONAL COMPUTATIONS~ */

		/* Axon sodium calcs */
		mod_prec m_inf_a = 1 /
			(1 + (expf((-30 - V_axon[target_cell]) / 5.5f)));
		mod_prec h_inf_a = 1 /
			(1 + (expf((-60 - V_axon[target_cell]) / -5.8f)));
		mod_prec tau_h_a = 1.5f * expf((-40 - V_axon[target_cell]) / 33);
		mod_prec dh_dt_a = (h_inf_a - Sodium_h_a[target_cell]) / tau_h_a;
		Sodium_m_a[target_cell] = m_inf_a;
		Sodium_h_a[target_cell] = Sodium_h_a[target_cell] + DELTA
			* dh_dt_a;

		/* Axon potassium calcs */
		mod_prec alpha_x_a = 0.13f * (V_axon[target_cell] + 25)
			/ (1 - expf(-(V_axon[target_cell] + 25) / 10));
		mod_prec beta_x_a = 1.69f * expf(-0.0125f * (V_axon[target_cell]
					+ 35));
		mod_prec x_inf_a = alpha_x_a / (alpha_x_a + beta_x_a);
		mod_prec tau_x_a = 1 / (alpha_x_a + beta_x_a);
		mod_prec dx_dt_a = (x_inf_a - Potassium_x_a[target_cell])
			/ tau_x_a;
		Potassium_x_a[target_cell] = 0.05f * dx_dt_a
			+ Potassium_x_a[target_cell];

		/* Axonal voltage and current calcs */
		mod_prec I_Na_a = G_NA_A * powf(Sodium_m_a[target_cell], 3)
			* Sodium_h_a[target_cell] * (V_axon[target_cell] - V_NA);
		mod_prec I_la = G_LA * (V_axon[target_cell] - V_L);
		mod_prec I_sa = (G_INT / P2) * (V_axon[target_cell]
				- V_soma[target_cell]);
		mod_prec I_K_a = G_K_A * powf(Potassium_x_a[target_cell], 4) *
			(V_axon[target_cell] - V_K);
		mod_prec dVa_dt = (-(I_K_a + I_sa + I_la + I_Na_a)) / C_M;

		/* ~NEW VOLTAGES~ */
		V_dend[target_cell] = DELTA * dVd_dt + V_dend[target_cell];
		V_soma[target_cell] = DELTA * dVs_dt + V_soma[target_cell];
		V_axon[target_cell] = DELTA * dVa_dt + V_axon[target_cell];

	}
}

void simulate(int simulation_steps, infoli_conf_t *config)
{
	mod_prec iApp;
	cellCompParams cell_params = config->cell_params;
	int sim_step;

	struct infoli_perf_region *gjf_stats = infoli_perf_create_region("gap juction functions");
	struct infoli_perf_region *upd_stats = infoli_perf_create_region("update state");

	printf("simulation_steps: %d\n", simulation_steps);
	for (sim_step = 0; sim_step < simulation_steps; ++sim_step) {
		if (sim_step >= 20000 && sim_step < 20500 - 1)
			iApp = 6;
		else
			iApp = 0;

#ifndef G_CAL_FROM_FILE
		infoli_log_periodic(sim_step, " %d %.2f ", sim_step + 1, iApp);
#endif

		infoli_perf_start(gjf_stats);
		
		/* compute gap junction functions */
		gap_junction_functions(config->iApp, config->cellCount,
			config->iAppIn, config->V_dend, cell_params.neighId,
			cell_params.neighConductances,
			cell_params.total_amount_of_neighbours,
			config->I_c);

		infoli_perf_stop(gjf_stats);
		infoli_perf_start(upd_stats);
		
		/* compute the new state (voltages, channels, etc.) */
		update_state(config->cellCount, config->V_dend,
			config->Hcurrent_q, config->Calcium_r, config->Ca2Plus,
			config->Potassium_s, config->I_CaH, config->V_soma,
			config->I_c, config->iAppIn, config->Calcium_k,
			config->Calcium_l, config->Sodium_h, config->Sodium_m,
			config->Potassium_n, config->Potassium_p,
			config->Potassium_x_s, config->g_CaL,
			config->Sodium_m_a, config->Sodium_h_a, config->V_axon,
			config->Potassium_x_a
		);

		infoli_perf_stop(upd_stats);
		
		infoli_print_results_periodic(sim_step, config->cellCount, config->V_axon);
	}

	infoli_perf_print(gjf_stats);
	infoli_perf_print(upd_stats);
}
