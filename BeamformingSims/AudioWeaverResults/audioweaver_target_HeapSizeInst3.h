/* ----------------------------------------------------------------------
** Audio Weaver Target Heap Used Header File
** Created on 23-Feb-2025 21:18:13
** ------------------------------------------------------------------- */

#if defined(SINGLE_HEAP)
#define MASTER_HEAP_SIZE (51)
#define FASTB_HEAP_SIZE (0)
#define SLOW_HEAP_SIZE (0)
#define SHARED_HEAP_SIZE (45)
#else
#define MASTER_HEAP_SIZE (17)
#define FASTB_HEAP_SIZE (17)
#define SLOW_HEAP_SIZE (17)
#define SHARED_HEAP_SIZE (45)
#endif

