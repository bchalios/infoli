#ifndef __INFOLI_LOG_H__
#define __INFOLI_LOG_H__

#include "infoli.h"

/* flag enabling or disabling output , axon's voltage is the output at every step */
#define PRINTING 1

#ifdef PRINTING

/* initialize logging */
void infoli_init_log(const char *log_filename);

/* shutdown logging */
void infoli_close_log(void);

/* Print out a message to the output file */
void infoli_log(const char *format, ...);

/* Print out a message to the output file periodically
 * with a period of 'print_granularity' */
void infoli_log_periodic(int _step, const char *fmt, ...);

/* print simulation results */
void infoli_print_results(int cellCount, mod_prec *V_axon);

/* print simulation results periodically */
void infoli_print_results_periodic(int step, int cellCount, mod_prec *V_axon);

#else /* PRINTING */

#define infoli_init_log(log_filename) {}
#define infoli_close_log() {}
#define infoli_log(fmt...) {}
#define infoli_log_periodic(_step, fmt...) {}
#define infoli_print_results(_cellCount, _vaxon) {}
#define infoli_print_results_periodic(_step, _cellcount, _vaxon) {}

#endif /* PRINTING */



#endif /* __INFOLI_LOG_H__ */
