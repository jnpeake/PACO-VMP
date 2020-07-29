#pragma once

#include <immintrin.h>

#ifndef _PLATFORM_INC_
#define _PLATFORM_INC_
//#define __AVX512F__
#define __AVX2__


#ifdef __AVX2__
#define AVX2
#define BITS 32
#define _VECSIZE 8
#elif defined __AVX512F__
#define _VECSIZE 16
#define BITS 64
#define AVX512
#else
#define SISD
#define BITS 32
#define _VECSIZE 8
#endif
#ifdef __AVX2__
#define ALLOC( _x ) _mm_malloc( (_x), BITS )
#define FREE( _x ) _mm_free( _x )
#define USE_OMP
#else
#define ALLOC( _x ) malloc( (_x) )
#define FREE( _x ) free( _x )
#endif

#ifdef __GNUC__
#define ALIGN(a) a __attribute__ ((aligned (BITS)))
#else
#define ALIGN(a) __declspec(align(BITS)) a
#endif

#ifdef USE_OMP
#include <omp.h>
#endif

#endif
