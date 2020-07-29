#pragma once

#ifndef _VECTOR_INC_
#define _VECTOR_INC_
#include <immintrin.h>
#include <cstdio> 
#include "platform.h"


class Vector
{

/* 
	The specific types of fields of the Vector class change depending on the 
	AVX instructions being used, but Vectors have three possible uses: storing
	floats, storing integers, or storing a mask. Each Vector, regardless of
	it's intended purpose, contains all three of the following:

	A float array/vector (values/AVXVec)
	An integer array/vector (iValues/AVXIntVec)
	A mask value. The SISD Vector stores masks as an iValues array,
	AVX512 uses the exclusive mmask16 type, and AVX2 uses an integer.
*/
public:
#ifdef SISD //fields used for the non-AVX vector implementation
	float AVXVec[8];
	int AVXIntVec[8];
	int vectorSize = 8;
#elif defined AVX512 //fields used for AVX512 compatible hardware
	__m512 AVXVec;
	__m512i AVXIntVec;
	__mmask16 maskVec;
#elif defined AVX2 //fields used for AVX2 compatible hardware
	__m256 AVXVec;
	__m256i AVXIntVec;
	int mask;
#endif
/*
	The + function for Vectors, adding the values of two Vector
	objects together lane-by-lane

	@param v | The Vector that will be added to the Vector that called the method
	@return | Returns a vector containing the results of adding the vectors
*/
	Vector operator+(const Vector& v) const
	{
		Vector vector;
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			vector.AVXVec[i] = this->AVXVec[i] + v.AVXVec[i];
		}

		return vector;
#elif defined AVX512
		vector.AVXVec = _mm512_add_ps(this->AVXVec, v.AVXVec);
		return vector;
#elif (defined AVX || defined AVX2)
		vector.AVXVec = _mm256_add_ps(this->AVXVec, v.AVXVec);
		return vector;

#endif
	}

/*
	The - function for Vectors, subtracting the values of the parameter Vector
	from the values of the calling Vector lane-by-lane

	@param v | The Vector that will be subtracted from the Vector that called the method
	@return | Returns a vector containing the results of subtracting the vectors
*/
	Vector operator-(const Vector& v) const
	{
		Vector vector;
#ifdef SISD
		for (int i = 0; i < vectorSize; i++)
		{
			vector.AVXVec[i] = (this->AVXVec[i] - v.AVXVec[i]);
		}

		return vector;
#elif defined AVX512
		vector.AVXVec = _mm512_sub_ps(this->AVXVec, v.AVXVec);
#elif (defined AVX || defined AVX2)
		vector.AVXVec = _mm256_sub_ps(this->AVXVec, v.AVXVec);
		return vector;

#endif
	}

/*
	The * function for Vectors, multiplying the values of two Vector
	objects together lane-by-lane

	@param v | The Vector that will be multiplied by the Vector that called the method
	@return | Returns a vector containing the results of multiplying the vectors
*/
	Vector operator*(const Vector& v) const
	{
		Vector vector;
#ifdef SISD	
		for (int i = 0; i < vectorSize; i++)
		{
			vector.AVXVec[i] = (this->AVXVec[i] * v.AVXVec[i]);
		}

		return vector;
#elif defined AVX512
		vector.AVXVec = _mm512_mul_ps(this->AVXVec, v.AVXVec);
		return vector;
#elif (defined AVX || defined AVX2)
		vector.AVXVec = _mm256_mul_ps(this->AVXVec, v.AVXVec);
		return vector;

#endif
	}

	Vector operator*(const int v) const
	{
		Vector vector;
		Vector intVec;

		intVec.set1(v);
#ifdef SISD	
		for (int i = 0; i < vectorSize; i++)
		{
			vector.AVXVec[i] = (this->AVXVec[i] * intVec.AVXVec[i]);
		}

		return vector;
#elif defined AVX512
		vector.AVXVec = _mm512_mul_ps(this->AVXVec, intVec.AVXVec);
		return vector;
#elif (defined AVX || defined AVX2)
		
		vector.AVXVec = _mm256_mul_ps(this->AVXVec, intVec.AVXVec);
		return vector;

#endif
	}


	Vector operator/(const Vector& v) const
	{
		Vector vector;
#ifdef SISD	
		for (int i = 0; i < vectorSize; i++)
		{
			vector.AVXVec[i] = (this->AVXVec[i] / v.AVXVec[i]);
		}

		return vector;
#elif defined AVX512
		vector.AVXVec = _mm512_div_ps(this->AVXVec, v.AVXVec);
		return vector;
#elif (defined AVX || defined AVX2)
		vector.AVXVec = _mm256_div_ps(this->AVXVec, v.AVXVec);
		return vector;

#endif
	}


/*
	Enforces a maximum value upon a vector, reducing all values that are above
	the maximum to the maximum.

	@param maxVal | The maximum allowed value
*/
	void vecMax(float maxVal)
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			if (this->AVXVec[i] > maxVal)
			{
				this->AVXVec[i] = maxVal;
			}
		}
#elif defined AVX512
		__declspec(align(64)) float maxValAr[16] = { maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal };
		__m512 maxValVec = _mm512_load_ps(maxValAr);
		this->AVXVec = _mm512_min_ps(this->AVXVec, maxValVec);
#elif (defined AVX || defined AVX2)
		float maxValAr[8] __attribute__ ((aligned (32))) = { maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal};
		__m256 maxValVec = _mm256_load_ps(maxValAr);
		this->AVXVec = _mm256_min_ps(this->AVXVec, maxValVec);
#endif

	}
/*
	Enforces a minimum value upon a vector, increasing all values that are below
	the minimum to the minimum.

	@param minVal | The minimum allowed value
*/
	void vecMin(float minVal)
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			if (this->AVXVec[i] < minVal)
			{
				this->AVXVec[i] = minVal;
			}
		}
#elif defined AVX512
		__declspec(align(64)) float minValAr[16] = { minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal };
		__m512 minValVec = _mm512_load_ps(minValAr);
		this->AVXVec = _mm512_max_ps(this->AVXVec, minValVec);
#elif (defined AVX || defined AVX2)
		float minValAr[8] __attribute__ ((aligned (32))) = { minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal};
		__m256 minValVec = _mm256_load_ps(minValAr);
		this->AVXVec = _mm256_max_ps(this->AVXVec, minValVec);
#endif
	}

/*
	Sets each Vector lane to one value

	@param setValue | The value that the Vector lanes will be set to
*/
	void set1(float setValue)
	{
#ifdef SISD
		for (int i = 0; i < vectorSize; i++)
		{
			this->AVXVec[i] = setValue;
		}
#elif defined AVX512
		this->AVXVec = _mm512_set1_ps(setValue);
#elif (defined AVX || defined AVX2)
		this->AVXVec = _mm256_set1_ps(setValue);

#endif
	}

/*
	Loads float data from a location in memory (normally an array) into a Vector

	@param source | The memory location of the data to be loaded
*/
	void load(float* source)
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{

			this->AVXVec[i] = source[i];


		}
#elif defined AVX512
		this->AVXVec = _mm512_load_ps(source);
#elif (defined AVX || defined AVX2)
		this->AVXVec = _mm256_load_ps(source);

#endif
	}

	void load_u(float* source)
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{

			this->AVXVec[i] = source[i];


		}
#elif defined AVX512
		this->AVXVec = _mm512_load_ps(source);
#elif (defined AVX || defined AVX2)
		this->AVXVec = _mm256_loadu_ps(source);

#endif
	}

/*
	Loads integer data from a location in memory (normally an array) into a Vector

	@param source | The memory location of the data to be loaded
*/
	
	void load(int* source)
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			this->AVXVec[i] = source[i];
		}
#elif defined AVX512
		this->AVXIntVec = _mm512_load_epi32(source);
#elif (defined AVX || defined AVX2)
		this->AVXIntVec = _mm256_load_si256(((const __m256i *)source));
#endif
	}

/*
	Loads unsigned integer data from a location in memory (normally an array) into a Vector

	@param source | The memory location of the data to be loaded
*/

	void load(unsigned int * source)
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			this->AVXVec[i] = source[i];
		}
#elif defined AVX512
		this->AVXIntVec = _mm512_load_epi32(source);
#elif defined AVX2
#endif
	}

/*
	Loads float data from a location in memory (normally an array) into a Vector. Applies
	a mask to the data, loading data from the calling Vector when the corresponding
	mask value is 1, and data from the src Vector when the corresponding value is 0.

	@param src | A vector containing data that will be added to the new vector is the
	corresponding mask value is 0
	@param mask | A vector of mask values (0 or 1) that determines which values will come
	from the calling Vector, and which values will come from src
	@param mem_addr | The memory location of the data to be loaded
*/
	
	void maskLoad(Vector &src, Vector &mask, float* mem_addr) 
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			if (mask.AVXVec[i] == 1)
				this->AVXVec[i] = mem_addr[i];

			else
				this->AVXVec[i] = src.AVXVec[i];
		}
#elif defined AVX512
		this->AVXVec = _mm512_mask_load_ps(src.AVXVec, mask.maskVec, mem_addr);
#elif (defined AVX || defined AVX2)

		Vector resultVector;
		resultVector.AVXVec = _mm256_load_ps(mem_addr);
		this->AVXVec = _mm256_blendv_ps(resultVector.AVXVec, src.AVXVec, mask.AVXVec);

#endif
	}

	


private:


};

Vector int2mask(int maskInt);
Vector mask_mov(const Vector& v1, const Vector& bitMask, const Vector& v2);
Vector gtMask(const Vector& v1, const Vector& v2);
Vector ltMask(const Vector& v1, const Vector& v2);
Vector vecRandom(Vector& rC0, Vector& rC1, Vector& factors, Vector& rSeed);
Vector pow(const Vector& v1, const Vector& v2);
Vector max (const Vector& v1, const Vector& v2);
Vector abs(const Vector& v1);
void seedVecRandom(Vector& rC0, Vector& rC1, Vector& factors, int *seeds, Vector& rSeed);
void maxLocStep(Vector &oldWeights, Vector &oldIndices, Vector &newWeights, Vector &newIndices);
int reduceMax(Vector &curWeights, Vector &curIndices);
void store(float* loc, const Vector& v1);
void store(int* loc, const Vector& v1);
void printVec(Vector v1);

#endif