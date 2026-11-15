/* Private to the pinned C target; never included by host/library clients.
 * -fexceptions permits a fatal library error to unwind to its C++ owner.
 */
#pragma once
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static inline uint64_t ramag_wfa_load_u64(const void* p) {
  uint64_t value; memcpy(&value,p,sizeof(value)); return value;
}
static inline void ramag_wfa_store_u64(void* p,uint64_t value) {
  memcpy(p,&value,sizeof(value));
}
void* ramag_wfa_malloc(size_t);
void* ramag_wfa_calloc(size_t,size_t);
void* ramag_wfa_realloc(void*,size_t);
void ramag_wfa_free(void*);
_Noreturn void ramag_wfa_exit(int);
#define malloc ramag_wfa_malloc
#define calloc ramag_wfa_calloc
#define realloc ramag_wfa_realloc
#define free ramag_wfa_free
#define exit ramag_wfa_exit
