#ifndef __INFOLI_PERF_H__
#define __INFOLI_PERF_H__

struct infoli_perf_region;

/** \brief Create a new named stats region
 * 
 * \param name is the name of the new region
 * 
 * \return a handle to be used for statistics-related operations
 */
struct infoli_perf_region *
infoli_perf_create_region(const char *name);

/** \brief Start measurements for the region
 * 
 * \param region is the stats region to start
 */
void infoli_perf_start(struct infoli_perf_region *region);

/** \brief Stops measurements for the region
 * 
 * \param region is the stats region to stop
 */
void infoli_perf_stop(struct infoli_perf_region *region);

/** \brief Print measurements for the region
 * 
 * \param region is the stats region to print
 */
void infoli_perf_print(struct infoli_perf_region *region);


#endif /* __INFOLI_PERF_H__ */
