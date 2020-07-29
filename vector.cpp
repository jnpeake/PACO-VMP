/**
    vector.cpp
    Purpose: Contains all vector instructions for AVX512 and AVX2, as well as
	SISD variants for incompatible hardware.

    @author Joshua Peake & Huw Lloyd
    @version 1.0
*/

#include "vector.h"
#include <iostream>
#include "platform.h"

/*
	Converts an integer to either an 8-wide or 16-wide mask vector, which
	is equivalent to the binary representation of that integer, for example:
	14291 becomes:
	0011011111010011 (AVX512) or 11010011 (AVX2)

	@param maskInt | The integer to be converted to a mask
	@return | The vector object containing the mask

*/
Vector int2mask(int maskInt) // NO AVX --------------------------
{
#ifdef SISD

	Vector mask;
	for (int i = 0; i < _VECSIZE; ++i) {
		mask.AVXVec[i] = (maskInt >> i) & 1;
	}


	return mask;

#elif defined AVX512
	Vector resultVector;
	resultVector.maskVec = _mm512_int2mask(maskInt);
	return resultVector;

#elif (defined AVX || defined AVX2)
	Vector resultVector;
	float mask[8];
	for (int i = 0; i < 8; i++) {
		if ((maskInt >> i) & 1)
		{
			mask[i] = 1.0f;
		}

		else
		{
			mask[i] = -1.0f;
		}
	}

	resultVector.AVXVec = _mm256_setr_ps(mask[0], mask[1], mask[2], mask[3], mask[4], mask[5], mask[6], mask[7]);

	return resultVector;

#endif
}

/*
	Moves a vector from one vector to another, using a mask. Values are moved from
	v1 to resultVector if corresponding lane in bitmask is 0, or moved from v2 to
	resultVector if corresponding lane in bitmask is 1.

	@param v1 | A vector of values, usually weight data
	@param bitMask | A vector containing a bitMask, used to filter values from v1
	@param v2 | A vector usually containing one value multiple times, used for
	when values from v1 are masked 
	@return | The filtered vector containing a combination of values from v1 and v2

*/
Vector mask_mov(const Vector& v1, const Vector& bitMask, const Vector& v2) // NO AVX -----------------
{
	Vector resultVector;
#ifdef SISD
	for (int i = 0; i < _VECSIZE; i++)
	{
		if (bitMask.AVXVec[i] == 0)
		{
			resultVector.AVXVec[i] = v1.AVXVec[i];
		}

		else
		{
			resultVector.AVXVec[i] = v2.AVXVec[i];
		}
	}

	return resultVector;

#elif defined AVX512
	resultVector.AVXVec = _mm512_mask_mov_ps(v1.AVXVec, bitMask.maskVec, v2.AVXVec);
	return resultVector;
#elif (defined AVX || defined AVX2)
	resultVector.AVXVec = _mm256_blendv_ps(v2.AVXVec, v1.AVXVec, bitMask.AVXVec);
	return resultVector;
#endif 
}

/*
	Combines two vectors into a new vector, comparing vectors lane-by-lane and retaining the
	largest value.

	@param v1 | A vector of float values
	@param v2 | A vector of float values
	@return | The filtered vector containing a combination of values from v1 and v2

*/
Vector gtMask(const Vector& v1, const Vector& v2) // NO AVX -------------------
{
	Vector resultMask;
#ifdef SISD

	for (int i = 0; i < _VECSIZE; i++)
	{
		if (v1.AVXVec[i] > v2.AVXVec[i])
		{
			resultMask.AVXVec[i] = 1;
		}

		else
		{
			resultMask.AVXVec[i] = 0;
		}
	}

	return resultMask;

#elif defined AVX512
	resultMask.maskVec = _mm512_cmp_ps_mask(v1.AVXVec, v2.AVXVec, _MM_CMPINT_GT);
	return resultMask;
#elif (defined AVX || defined AVX2)
	resultMask.AVXVec = _mm256_cmp_ps(v1.AVXVec, v2.AVXVec, _CMP_GT_OS);
	return resultMask;
#endif
}

/*
	Combines two vectors into a new vector, comparing vectors lane-by-lane and retaining the
	lowest value.

	@param v1 | A vector of float values
	@param v2 | A vector of float values
	@return | The filtered vector containing a combination of values from v1 and v2

*/
Vector ltMask(const Vector& v1, const Vector& v2) // NO AVX --------------------------
{
	Vector resultMask;
#ifdef SISD
	for (int i = 0; i < _VECSIZE; i++)
	{
		if (v1.AVXVec[i] < v2.AVXVec[i])
		{
			resultMask.AVXVec[i] = 1;
		}

		else
		{
			resultMask.AVXVec[i] = 0;
		}
	}
	return resultMask;
#elif defined AVX512
	resultMask.maskVec = _mm512_cmp_ps_mask(v1.AVXVec, v2.AVXVec, _MM_CMPINT_LT);
	return resultMask;
#elif (defined AVX || defined AVX2)
	resultMask.AVXVec = _mm256_cmp_ps(v1.AVXVec, v2.AVXVec, _CMP_LE_OS);
	return resultMask;
#endif

}

/*
	Initialises variables required for random number generation. The constant values 1664525 and 1013904223
	cause integer overflow, allowing for true random number generation. The value 2.3283064e-10f will be used
	in a multiplication to reduce very large random numbers to a value between 0 and 1.


	@param rC0 | An empty vector that will be filled with the first constant int, 1664525
	@param rC1 | An empty vector that will be filled with the second constant int, 1013904223
	@param factor | An empty vector that will be filled with the constant float, 2.3283064e-10f
	@param seeds | An array of seeds generated from the input seed by ranluxgen
	@param rSeed | An empty vector that will be filled with seeds from the seeds parameter
*/
void seedVecRandom(Vector& rC0, Vector& rC1, Vector& factor, int *seeds, Vector& rSeed)
{
#ifdef SISD
	for (int i = 0; i < 8; i++)
	{
		rSeed.AVXIntVec[i] = seeds[i];
	}


#elif defined AVX512
	__declspec(align(64)) unsigned int c0[16] = { 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L,1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L };
	__declspec(align(64)) unsigned int c1[16] = { 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 	1013904223L, 1013904223L, 1013904223L,1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L };
	__declspec(align(64)) float factors[16] = { 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f,2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f };

	rSeed.AVXIntVec = _mm512_load_epi32(seeds);
	rC0.AVXIntVec = _mm512_load_epi32(&c0); //m512i
	rC1.AVXIntVec = _mm512_load_epi32(&c1); //m512i
	factor.AVXVec = _mm512_load_ps(factors);


#elif (defined AVX || defined AVX2)
	int c0 = 1664525L;
	int c1 = 1013904223L;
	float factors[8] __attribute__ ((aligned (32))) = { 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f};
	
	rSeed.AVXIntVec = _mm256_load_si256((__m256i*)seeds);
	rC0.AVXIntVec = _mm256_set1_epi32(c0); //m512i
	rC1.AVXIntVec = _mm256_set1_epi32(c1); //m512i
	factor.AVXVec = _mm256_load_ps(factors);

	


#endif // AVX512

}

/*
	Generates a vector of random numbers between 0 and 1, created using the seeds created by the ranluxgen and
	3 constants (rC0, rC1 and factor). Used by ant.cpp when performing edge selection in csRoulette().


	@param rC0 | A vector containing the first constant int, 1664525
	@param rC1 | A vector containing the second constant int, 1013904223
	@param factor | A vector containing the constant float, 2.3283064e-10f
	@param rSeed | A vector containing seeds generated by ranluxgen

	@return | A vector contained with 8 or 16 random numbers between 0 and 1
*/

Vector vecRandom(Vector& rC0, Vector& rC1, Vector& factor, Vector& rSeed)
{
	Vector r;
#ifdef SISD
	for (int i = 0; i < 8; i++)
	{

		rSeed.AVXIntVec[i] = rSeed.AVXIntVec[i] * 1664525L + 1013904223L;
		r.AVXVec[i] = (float)rSeed.AVXIntVec[i] * 2.328306437087974e-10;
		r.AVXVec[i] += 0.5f;
	}
	return r;

#elif defined AVX512
	rSeed.AVXIntVec = _mm512_mullo_epi32(rC0.AVXIntVec, rSeed.AVXIntVec);
	rSeed.AVXIntVec = _mm512_add_epi32(rC1.AVXIntVec, rSeed.AVXIntVec);
	// convert to float in range 0 to 1 and return
	__m512 returnValue = _mm512_cvt_roundepu32_ps(rSeed.AVXIntVec, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
	r.AVXVec = _mm512_mul_ps(returnValue, factor.AVXVec);
	return r;

#elif defined AVX2
	__m256 addedVal;
	addedVal =_mm256_set1_ps(0.5f);
	rSeed.AVXIntVec = _mm256_mullo_epi32(rC0.AVXIntVec, rSeed.AVXIntVec);
	rSeed.AVXIntVec = _mm256_add_epi32(rC1.AVXIntVec, rSeed.AVXIntVec);
	// convert to float in range 0 to 1 and return
	__m256 returnValue = _mm256_cvtepi32_ps(rSeed.AVXIntVec);
	
	returnValue = _mm256_mul_ps(returnValue, factor.AVXVec);
	r.AVXVec = _mm256_add_ps(returnValue,addedVal);
	return r;
	

#endif
}

/*
	Compares values in the oldWeights vector (the vector of the previous largest values) with
	values in the newWeights vector (the vector of values from the current ACO iteration) and
	replaces values in oldWeights with larger values in the same vector lane from newWeights.
	The corresponding indexes for each weight are stored in oldIndices and newIndices
	respectively.

	@param oldWeights | A vector containing weight values for the previous largest weights
	@param oldIndices | A vector containing the corresponding indices for weights in oldWeights
	@param newWeights | A vector containing weight values for the largest weights in the current
	ACO iteration
	@param newIndices | A vector containing the corresponding indices for weights in newWeights

*/
void maxLocStep(Vector &oldWeights, Vector &oldIndices, Vector &newWeights, Vector &newIndices)
{
#ifdef SISD
	Vector maxMask = gtMask(newWeights, oldWeights);
	oldWeights = mask_mov(oldWeights, maxMask, newWeights);
	oldIndices = mask_mov(oldIndices, maxMask, newIndices);
#elif defined AVX512
	Vector maxMask;
	maxMask.maskVec = _mm512_cmp_ps_mask(newWeights.AVXVec, oldWeights.AVXVec, _MM_CMPINT_GT);
	oldWeights.AVXVec = _mm512_mask_mov_ps(oldWeights.AVXVec, maxMask.maskVec, newWeights.AVXVec);
	oldIndices.AVXVec = _mm512_mask_mov_ps(oldIndices.AVXVec, maxMask.maskVec, newIndices.AVXVec);
#elif (defined AVX || defined AVX2)
	Vector maxMask;
	maxMask = gtMask(newWeights, oldWeights);
	oldWeights = mask_mov(newWeights, maxMask, oldWeights);
	oldIndices = mask_mov(newIndices, maxMask, oldIndices);
#endif


}

/*
	Reduces the vector containing the largest weight values from every iteration (curWeights), 
	determining the largest value in the vector and changing all vector lanes to that value, with 
	identical operations being performed on the vector of the indeces associated with those weights

	@param curWeights | A vector containing weight values for largest weights from all iterations
	@param curIndices | A vector containing the corresponding indices for weights in curWeights
	@return | A vector with each lane containing the largest weight value in curWeights
*/
int reduceMax(Vector &curWeights, Vector &curIndices)
{
#ifdef SISD

	int highestIndex = -1;
	float highestValue = -1.0f;
	for (int i = 0; i < _VECSIZE; i++)
	{
		if (curWeights.AVXVec[i] > highestValue)
		{
			highestValue = curWeights.AVXVec[i];
			highestIndex = curIndices.AVXVec[i];
		}
	}

	return highestIndex;
#elif defined AVX512
	// return a vector with all elements equal to ivec[imax] where
	// valvec[imax] is largest element of valvec
	__m512 permVal;
	__m512 permInd;
	__mmask16 maxMask;
	// swap with neighbour 
	permVal = _mm512_swizzle_ps(curWeights.AVXVec, _MM_SWIZ_REG_CDAB);
	permInd = _mm512_swizzle_ps(curIndices.AVXVec, _MM_SWIZ_REG_CDAB);
	maxMask = _mm512_cmp_ps_mask(curWeights.AVXVec, permVal, _MM_CMPINT_GT);
	curWeights.AVXVec = _mm512_mask_mov_ps(permVal, maxMask, curWeights.AVXVec);
	curIndices.AVXVec = _mm512_mask_mov_ps(permInd, maxMask, curIndices.AVXVec);
	// swap pairs
	permVal = _mm512_swizzle_ps(curWeights.AVXVec, _MM_SWIZ_REG_BADC);
	permInd = _mm512_swizzle_ps(curIndices.AVXVec, _MM_SWIZ_REG_BADC);
	maxMask = _mm512_cmp_ps_mask(curWeights.AVXVec, permVal, _MM_CMPINT_GT);
	curWeights.AVXVec = _mm512_mask_mov_ps(permVal, maxMask, curWeights.AVXVec);
	curIndices.AVXVec = _mm512_mask_mov_ps(permInd, maxMask, curIndices.AVXVec);
	// swap lanes
	permVal = _mm512_permute4f128_ps(curWeights.AVXVec, 0xB1); // 2, 3, 0, 1
	permInd = _mm512_permute4f128_ps(curIndices.AVXVec, 0xB1);
	maxMask = _mm512_cmp_ps_mask(curWeights.AVXVec, permVal, _MM_CMPINT_GT);
	curWeights.AVXVec = _mm512_mask_mov_ps(permVal, maxMask, curWeights.AVXVec);
	curIndices.AVXVec = _mm512_mask_mov_ps(permInd, maxMask, curIndices.AVXVec);
	// swap pairs of lanes
	permVal = _mm512_permute4f128_ps(curWeights.AVXVec, 0x4E); // 1, 0, 3, 2
	permInd = _mm512_permute4f128_ps(curIndices.AVXVec, 0x4E);
	maxMask = _mm512_cmp_ps_mask(curWeights.AVXVec, permVal, _MM_CMPINT_GT);
	curIndices.AVXVec = _mm512_mask_mov_ps(permInd, maxMask, curIndices.AVXVec);
	// all elements of ivec now contain index of maximum
	return curIndices.AVXVec[0];
#elif (defined AVX || defined AVX2)
	
	__m256 permVal;
	__m256 permInd;
	__m256 maxMask;
	float result[8] __attribute__ ((aligned (32)));
	permVal = _mm256_permute_ps(curWeights.AVXVec, _MM_SHUFFLE(2,3,0,1)); //01001110
	permInd = _mm256_permute_ps(curIndices.AVXVec, _MM_SHUFFLE(2, 3, 0, 1)); //01001110
	maxMask = _mm256_cmp_ps(curWeights.AVXVec, permVal, _CMP_GT_OS);
	curWeights.AVXVec = _mm256_blendv_ps(permVal, curWeights.AVXVec, maxMask);
	curIndices.AVXVec = _mm256_blendv_ps(permInd, curIndices.AVXVec, maxMask);

	permVal = _mm256_permute_ps(curWeights.AVXVec, _MM_SHUFFLE(1,0,3,2)); //01001110
	permInd = _mm256_permute_ps(curIndices.AVXVec, _MM_SHUFFLE(1, 0, 3, 2)); //01001110
	maxMask = _mm256_cmp_ps(curWeights.AVXVec, permVal, _CMP_GT_OS);
	curWeights.AVXVec = _mm256_blendv_ps(permVal, curWeights.AVXVec, maxMask);
	curIndices.AVXVec = _mm256_blendv_ps(permInd, curIndices.AVXVec, maxMask);

	permVal = _mm256_permute2f128_ps(curWeights.AVXVec, curWeights.AVXVec, 0b0001);
	permInd = _mm256_permute2f128_ps(curIndices.AVXVec, curIndices.AVXVec, 0b0001);
	maxMask = _mm256_cmp_ps(curWeights.AVXVec, permVal, _CMP_GT_OS);
	curWeights.AVXVec = _mm256_blendv_ps(permVal, curWeights.AVXVec, maxMask);
	curIndices.AVXVec = _mm256_blendv_ps(permInd, curIndices.AVXVec, maxMask);
	store(result, curIndices);
	return result[0];
#endif

}

/*
	Takes float values that are currently stored in a vector (most likely an array) and stores them
	in a memory (most likely an array).

	@param loc | The memory location where the data is to be stored
	@param v1 | A vector containing data that will be stored
*/
void store(float* loc, const Vector& v1)
{
#ifdef SISD

	for (int i = 0; i < _VECSIZE; i++)
	{
		loc[i] = v1.AVXVec[i];
	}

#elif defined AVX512
	_mm512_store_ps(loc, v1.AVXVec);
#elif (defined AVX || defined AVX2)
	_mm256_store_ps(loc, v1.AVXVec);
#endif
}

/*
	Takes integer values that are currently stored in a vector (most likely an array) and stores them
	in a memory (most likely an array).

	@param loc | The memory location where the data is to be stored
	@param v1 | A vector containing data that will be stored
*/
void store(int* loc, const Vector& v1)
{
#ifdef SISD

	for (int i = 0; i < _VECSIZE; i++)
	{
		loc[i] = v1.AVXVec[i];
	}

#elif defined AVX512
	_mm512_store_ps(loc, v1.AVXVec);
#elif (defined AVX || defined AVX2)
	_mm256_store_si256((__m256i *)loc, v1.AVXIntVec);
#endif
}

Vector max (const Vector& v1, const Vector& v2)
{
	#ifdef SISD
	Vector result;
	for(int i = 0; i < _VECSIZE; i++)
	{
		if(v1.AVXVec[i] > v2.AVXVec[i])
		{
			result.AVXVec[i] = v1.AVXVec[i];
		}
		else
		{
			result.AVXVec[i] = v2.AVXVec[i];
		}
	}

	return result;
	#elif (defined AVX || defined AVX2)
	Vector result;
	result.AVXVec = _mm256_max_ps(v1.AVXVec, v2.AVXVec);
	return result;
	#endif
}


Vector abs(const Vector& v1)
{
	#ifdef SISD
	Vector result;
	for(int i = 0; i < _VECSIZE; i++)
	{
		result.AVXVec[i] = std::abs(v1.AVXVec[i]);
	}
	return result;

	#elif (defined AVX || defined AVX2)
	Vector signMask, resultVector;
    signMask.set1(-0.0f);
	resultVector.AVXVec =  _mm256_andnot_ps(signMask.AVXVec, v1.AVXVec);
	return resultVector;
	#endif
}


