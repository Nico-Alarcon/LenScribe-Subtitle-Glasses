/* ----------------------------------------------------------------------
** Audio Weaver Target Heap Used Header File
** Created on 04-Mar-2025 00:43:44
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

