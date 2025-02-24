/* ----------------------------------------------------------------------
** Audio Weaver Target Heap Used Header File
** Created on 23-Feb-2025 21:18:13
** ------------------------------------------------------------------- */

#if defined(SINGLE_HEAP)
#define MASTER_HEAP_SIZE (3070)
#define FASTB_HEAP_SIZE (0)
#define SLOW_HEAP_SIZE (0)
#define SHARED_HEAP_SIZE (45)
#else
#define MASTER_HEAP_SIZE (1885)
#define FASTB_HEAP_SIZE (988)
#define SLOW_HEAP_SIZE (197)
#define SHARED_HEAP_SIZE (45)
#endif

