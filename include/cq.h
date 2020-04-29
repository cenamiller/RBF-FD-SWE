#pragma once
///////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2018-2019 CacheQ Systems Inc. All rights reserved.
//
// cq.h -- general cacheq includes
///////////////////////////////////////////////////////////////////////////////

// include function declarations before we #define to replace them
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifdef __cplusplus
#include <map>
#include "mem.h"

extern "C"
{
#else
#include <stdbool.h>
#endif
///////////////////////////////////////////////////////////////////////////////
// cq_*() functions are defined in mem.cpp
#ifdef __clang__
#define CQ_POOL(p) __attribute__((address_space(p)))
#else
#define CQ_POOL(p)
#endif

#define CQ_POOL_NONE    (0)
#define CQ_POOL_CONST   (1)
#define CQ_POOL_GLOBAL  (2)
#define CQ_POOL_USER    (3)

#ifdef _WIN32
    #define CLANG_NOINLINE    __declspec(noinline)
    #define FORCE_INLINE      __forceinline
    #ifdef __always_inline
        #undef __always_inline
    #endif
#elif defined(__APPLE__)
    #define CLANG_NOINLINE    __attribute__((noinline))
    // no equivalent for __always_inline on mac
    #define FORCE_INLINE   static inline
#else
    #define CLANG_NOINLINE    __attribute__((noinline))
    #define FORCE_INLINE      __attribute__((always_inline))
    // __always_inline should already exist on Linux
#endif

#define __never_inline CLANG_NOINLINE
#ifndef __always_inline
    #define __always_inline FORCE_INLINE
#endif

void*  cq_malloc(int pool, size_t size);
void*  cq_calloc(int pool, size_t n, size_t size);
void*  cq_realloc(void* p, size_t newsize);
void   cq_free(void* p);
void   cq_mstats(int verbosity);
size_t cq_msize(void* p);
void   cq_malloc_verify(bool verbose);
void   cq_malloc_debug(int level);
char*  cq_strdup(int pool, const char* from);
const  char* cq_getenv(const char* name);

uint64_t cq_canonicalize_pointer(double p);
void*  cq_localize_pointer(uint64_t p);

void   cq_exit(int status);
void __cq__exit(int status); // just calls cq_exit.  Here for hysterical reasons

void   cq_halt();                       // Trapped in infinite loop

// number of nanoseconds expired since the first call to this function
uint64_t cq_nstime();

// Use our own ver4ion of printf for FPGA simulation runs...
int cq_printf(const char * format, ...);
#ifdef CQ_SIMULATION
#define printf cq_printf
#endif

// Overridden std lib memory functions to allocate from pool 0
#define malloc(sz) cq_malloc(CQ_POOL_GLOBAL, (sz))
#define calloc(n,sz) cq_calloc(CQ_POOL_GLOBAL, (n), (sz))
#define realloc cq_realloc
#define free cq_free
#define getenv cq_getenv

#if defined(strdup)
#undef strdup
#endif
#define strdup(s) cq_strdup(CQ_POOL_GLOBAL, s)

#define exit cq_exit


///////////////////////////////////////////////////////////////////////////////
// psuedo function calls to tell the compiler to do something.  Someday to
// be replaced by #pragmas

// Allow specification of unroll in the c code.
#ifdef CACHEQ
extern void _cq_unroll(int unrollcount, ...);   // second argument of 1 will
                                                // prevent from checking for
                                                // leftover iterations.

#else
#define _cq_unroll(a, ...)
#endif

// Tell the compiler to ignore loop-carried dependencies
extern void _cq_ignorelcd_pool(int poolid);     // Ignore only LCDs on a 
                                                // specific memory pool
#define _cq_ignorelcd() _cq_ignorelcd_pool(0);  // Ignore any memory LCDs

///////////////////////////////////////////////////////////////////////////////
// Random number generators
uint64_t cq_trng();     // TRNG on the FPGA.  calls cq_rand64() on CPU.

__always_inline uint32_t ___rotateleft(uint32_t v, int n)
{
   n = n % 32;
   v = v << n | v >> (32-n);
   return v;
}

__always_inline uint32_t ___rotateright(uint32_t v, int n)
{
   n = n % 32;
   v = v >> n | v << (32-n);
   return v;
}

__always_inline uint32_t dfsr_helper(uint32_t z, uint32_t seed)
{
   uint32_t r = ___rotateleft(z, 18*seed);
   r = ___rotateright((r << 18) ^ r, 13*seed);
   uint32_t b = ((r << 6) ^ r) >> 13;
   r = ((r & 4294967294U) << 18) ^ b;
   return r;
}

__always_inline uint64_t cq_zzrand64(uint64_t seed)
{
    const uint64_t PRIME64_1 = 11400714785074694791ULL;
    const uint64_t PRIME64_2 = 14029467366897019727ULL;
    const uint64_t PRIME64_3 =  1609587929392839161ULL;
    const uint64_t PRIME64_4 =  9650029242287828579ULL;
    const uint64_t PRIME64_5 =  2870177450012600261ULL;

    //first arg is 'sizeof(seed)', which will always be 8
    uint64_t start_1 = 8 + PRIME64_5;
    uint64_t start_2 = seed * PRIME64_2;
    uint64_t result_1 = ((start_2 << 31) | (start_2 >> 33)) * PRIME64_1;
    uint64_t start_3 = start_1 ^ result_1;
    uint64_t result_2 = (start_3 << 27) | (start_3 >> 37);

    uint64_t result = result_2 * PRIME64_1 + PRIME64_4;
    result ^= result >> 33;
    result *= PRIME64_2;
    result ^= result >> 29;
    result *= PRIME64_3;
    result ^= result >> 32;
    return result;
}

// These are ideal ones to use if you can be sure you give a different
// seed every time
__always_inline uint32_t cq_rand32(uint32_t seed)
{
    unsigned int z1 = dfsr_helper(seed, seed+1);

    const uint32_t PRIME32_2 = 2246822519U;
    const uint32_t PRIME32_3 = 3266489917U;
    const uint32_t PRIME32_4 =  668265263U;
    const uint32_t PRIME32_5 =  374761393U;

    uint32_t start_state = 4 + PRIME32_5 + seed * PRIME32_3;
    uint32_t z4 = (start_state << 17) | (start_state >> 15);
#ifdef CQ_STRONGER_RAND
    z4 *= PRIME32_4;
#endif
    z4 ^= z4 >> 15;
    z4 *= PRIME32_2;
    z4 ^= z4 >> 13;
#ifdef CQ_STRONGER_RAND
    z4 *= PRIME32_3;
    z4 ^= z4 >> 16;
#endif
    return (z1 ^ z4);
}

__always_inline uint64_t cq_rand64(uint64_t seed)
{
   uint32_t z1 = dfsr_helper(seed, seed+1);
   uint64_t z4 = cq_zzrand64(seed);
   return (z4 ^ z1);
}

__always_inline double cq_drand(uint32_t seed)
{
    uint32_t y = cq_rand32(seed);
    return y*(1.0/4294967295.0);
}

// These three should be okay within an inner loop but will create an LCD in the
// outer loop
__always_inline uint32_t cq_rand_r(unsigned int* seedp)
{
    uint32_t r = cq_rand32(*seedp);
    *seedp += 1;
    return r;
}

__always_inline uint64_t cq_rand64_r(uint64_t* seedp)
{
    uint64_t r = cq_rand64(*seedp);
    *seedp += 1;
    return r;
}

__always_inline double cq_drand_r(unsigned int* seed)
{
    uint32_t y = cq_rand_r(seed);
    return y*(1.0/4294967295.0);
}

// This is not really recommended, but probably will do okay in many cases
static inline int32_t cq_rand(int zero_or_seed)
{
   // I do it this way so I can have this static declared within the dataflow
   static unsigned int seed = 19;
   if (zero_or_seed != 0)
      return (seed = zero_or_seed);
   return cq_rand_r(&seed);
}

static inline void cq_srand(int seed)
{
   cq_rand(seed ? seed : 1);     // zero seed gets changed to 1
}

// Use this one for TRNG in FPGA.  In CPU it is a predictable PRNG.
//extern uint32_t cq_trand();

#define rand()  ((int)(cq_rand(0) & 0x7fffffff))
#define rand_r(p) ((int)cq_rand_r(p) & 0x7fffffff))
#define srand(i) cq_srand(i)


#ifdef CACHEQ
#endif

#ifdef __cplusplus
}

namespace cacheq {
///////////////////////////////////////////////////////////////////////////////
// This needs to be called all the time now.  It should only be called from the
// "main" CE.  It returns a pointer to the newly allocated copy of argv.
const MemoryPool::CanonPtr* init_shared_memory(int argc, const char* argv[]);
const MemoryPool::CanonPtr* init_shared_memory(
                        const void* const_begin = nullptr, size_t const_len = 0,
                        const void* globs_begin = nullptr, size_t globs_len = 0,
                        int argc = 0, const char* argv[] = nullptr);

///////////////////////////////////////////////////////////////////////////////
// Memory pool definition and creation

class MemInterface;          // Not currently implemented

class MemPoolDef
{
public:
    // Type of implentation in the hardware
    typedef enum {IT_CONST, IT_QRAM, IT_QCACHE, IT_SUPERCACHE, IT_DDR} ImplType;

public:
    ImplType             type;          // Implementation
    size_t               max_size;      // 0 means unknown and may only be used
                                        // for DRAM
    MemInterface*        interface;     // Hardware interface for cache or other

public:
    MemPoolDef(ImplType t, int64_t ms, MemInterface* mi=nullptr)
     : type(t)
     , max_size(ms)
     , interface(mi)
    {}
};

typedef std::map<int, MemPoolDef> MemPoolDefs;

void initialize_memory_pools(const MemPoolDefs& defs);

}

#endif
