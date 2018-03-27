#ifndef __UTIL_H__
#define __UTIL_H__

#ifndef NDEBUG
static inline void *xmalloc(size_t size)
{
	void *ret _mm_malloc(size);
	if (!_ret) {
		perror("__mm_malloc failed");
		exit(1);
	}

	return _ret;
}

#else /* NDEBUG */

#define xmalloc(_size) _mm_malloc((_size))

#endif /* NDEBUG */


#endif /* __UTIL_H__ */
