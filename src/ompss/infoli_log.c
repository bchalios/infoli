#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "infoli_log.h"

static FILE *infoli_out = NULL;
static int print_granularity = 1;
static const char *_logName = NULL;

void infoli_init_log(const char *log_filename)
{
	printf("Initializing logger\n");
	_logName = log_filename;
	infoli_out = fopen(log_filename, "w+");
	if(infoli_out == NULL){
		printf("Error: Couldn't create %s\n", log_filename);
		exit(EXIT_FAILURE);
	}

#ifdef G_CAL_FROM_FILE
	print_granularity=1;
#else
	print_granularity=1;
#endif
}

void infoli_close_log(void)
{
	printf("Closing logger\n");
	fclose(infoli_out);
	/* TODO: Check out the following */
	chmod(_logName, 0x01B6);
}

void infoli_log(const char *format, ...)
{
	va_list ap;

	assert(infoli_out);
	va_start(ap, format);
	vfprintf(infoli_out, format, ap);
}

void infoli_log_periodic(int step, const char *format, ...)
{
	if (step % print_granularity)
		return;

	va_list ap;

	assert(infoli_out);
	va_start(ap, format);
	vfprintf(infoli_out, format, ap);
}

void infoli_print_results(int cellCount, mod_prec *V_axon)
{
	int target_cell;

	for (target_cell = 0; target_cell < cellCount; target_cell++) {
#ifndef G_CAL_FROM_FILE
		infoli_log("%d : %.8f ", target_cell + 1, V_axon[target_cell]);
#else
		infoli_log("%.8f ", V_axon[target_cell]); 
#endif
	}
	infoli_log("\n");
}

void infoli_print_results_periodic(int step, int cellCount, mod_prec *V_axon)
{
	if (step % print_granularity)
		return;

	infoli_print_results(cellCount, V_axon);
}
