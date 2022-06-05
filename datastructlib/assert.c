#include "assert.h"
#include <stdint.h>
const Except_T Assert_Failed = {"Assertion Failed"};

void (assert)(int32_t e){
    assert(e);
}
