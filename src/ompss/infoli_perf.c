#include "infoli_perf.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define _POSIX_C_SOURCE 199309L
#include <time.h>

#define INFOLI_MAX_REGION_NAME 128

struct infoli_perf_region {
	char name[INFOLI_MAX_REGION_NAME];
	struct timespec tp;
	long running_ns;
};


struct infoli_perf_region *
infoli_perf_create_region(const char *name)
{
	struct infoli_perf_region *ret;
	
	ret = malloc(sizeof(*ret));
	if (!ret) {
		perror("Could not allocate stats region");	
		exit(1);
	}

	strncpy(ret->name, name, INFOLI_MAX_REGION_NAME);
	ret->running_ns = 0;
	
	return ret;
}

void infoli_perf_start(struct infoli_perf_region *region)
{
	clock_gettime(CLOCK_MONOTONIC_RAW, &region->tp);
}

void infoli_perf_stop(struct infoli_perf_region *region)
{
	struct timespec tp;
	
	clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
	
	region->running_ns += (tp.tv_sec - region->tp.tv_sec) * 1e9
				+ (tp.tv_nsec - region->tp.tv_nsec);
}

void infoli_perf_print(struct infoli_perf_region *region)
{
	printf("===== INFOLI_STATS_REGION: %s =====\n", region->name);
	printf("Time: %lf msec\n", region->running_ns / 1e6);
	printf("===================================\n");
}
