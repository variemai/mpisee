#undef assert
#ifdef NDEBUG
#define assert(e) ((void)0)
#else
#include "except.h"
#include <stdint.h>
extern void assert(int32_t e);
#define assert(e) ((void)((e)||(RAISE(Assert_Failed),0)))
#endif
extern const Except_T Assert_Failed;
