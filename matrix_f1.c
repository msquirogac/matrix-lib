#include "matrix.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#pragma GCC push_options
#pragma GCC optimize ("-O3")

uint32_t Matrix_AssignMemory(uint16_t d1, uint16_t d2, MFlt_t (**matrix1)[d2])
{
	*matrix1 = malloc(d1 * sizeof(**matrix1));
	return  (*matrix1 == NULL)? 1: 0;
}

uint32_t Matrix_AssignClearMemory(uint16_t d1, uint16_t d2, MFlt_t (**matrix1)[d2])
{
	*matrix1 = calloc(d1, sizeof(**matrix1));
  return  (*matrix1 == NULL)? 1: 0;
}

void Matrix_FreeMemory(uint16_t d1, uint16_t d2, MFlt_t (**matrix1)[d2])
{
	free(*matrix1);
	*matrix1 = NULL;
}

void Matrix_Copy(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2])
{
	memcpy(matrix2, matrix1, sizeof(*matrix1) * d1);
}

uint32_t Matrix_Compare(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], MFlt_t eps)
{
	MFlt_t *pIn1 = (MFlt_t*)matrix1;
	MFlt_t *pIn2 = (MFlt_t*)matrix2;
	uint32_t blkCnt;

	uint32_t result = 0;

	blkCnt = d1*d2;
	while(blkCnt > 0u)
	{
		result += C(fabs)(*pIn1 - *pIn2) > eps;
		pIn1++;
		pIn2++;
		blkCnt--;
	}

	return result;
}

void Matrix_GetRow(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], uint16_t row)
{
	memcpy(matrix2, &matrix1[row][0], sizeof(*matrix1));
}

void Matrix_SetRow(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], uint16_t row)
{
	memcpy(matrix1, &matrix2[row][0], sizeof(*matrix1));
}

void Matrix_GetColumn(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], uint16_t col)
{
	MFlt_t *pIn = (MFlt_t*)&matrix1[0][col];
	MFlt_t *pOut = (MFlt_t*)matrix2;
	uint32_t blkCnt;

	blkCnt = d1;
	while(blkCnt > 0u)
	{
		*pOut++ = *pIn;
		pIn += d2;
		blkCnt--;
	}
}

void Matrix_SetColumn(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], uint16_t col)
{
	MFlt_t *pIn = (MFlt_t*)matrix1;
	MFlt_t *pOut = (MFlt_t*)&matrix2[0][col];
	uint32_t blkCnt;

	blkCnt = d1;
	while(blkCnt > 0u)
	{
		*pOut = *pIn;
		pIn++;
		pOut += d2;
		blkCnt--;
	}
}


void Matrix_Eye(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2])
{
	MFlt_t *pOut = (MFlt_t*)matrix1;
	MFlt_t *pOutA = (MFlt_t*)matrix1;
	uint32_t numSamples;
	uint32_t blkCnt;

	numSamples = (uint32_t) d1 * d2;

	/* Loop unrolling */
	blkCnt = numSamples >> 2u;
	while(blkCnt > 0u)
	{
		pOut[0] = C(0.0);
		pOut[1] = C(0.0);
		pOut[2] = C(0.0);
		pOut[3] = C(0.0);

		pOut += 4u;
		blkCnt--;
	}

	/* If the numSamples is not a multiple of 4, compute any remaining output samples here.
	 **	No	loop	unrolling	is	used.	*/
	blkCnt = numSamples % 0x4u;
	while(blkCnt > 0u)
	{
		*pOut++ = C(0.0);
		blkCnt--;
	}

	pOut = pOutA;
	blkCnt = (uint32_t) d1;
	while(blkCnt > 0u)
	{
		*pOut = C(1.0);
		pOut += d2+1;
		blkCnt--;
	}
}

void Matrix_Zeros(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2])
{
	MFlt_t *pOut = (MFlt_t*)matrix1;
	uint32_t numSamples;
	uint32_t blkCnt;

	numSamples = (uint32_t) d1 * d2;

	/* Loop unrolling */
	blkCnt = numSamples >> 2u;
	while(blkCnt > 0u)
	{
		pOut[0] = C(0.0);
		pOut[1] = C(0.0);
		pOut[2] = C(0.0);
		pOut[3] = C(0.0);

		pOut += 4u;
		blkCnt--;
	}

	/* If the numSamples is not a multiple of 4, compute any remaining output samples here.
	 **	No	loop	unrolling	is	used.	*/
	blkCnt = numSamples % 0x4u;
	while(blkCnt > 0u)
	{
		*pOut++ = C(0.0);
		blkCnt--;
	}
}


void Matrix_Addition(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], MFlt_t (*matrix3)[d2])
{
	MFlt_t *pIn1 = (MFlt_t*)matrix1;
	MFlt_t *pIn2 = (MFlt_t*)matrix2;
	MFlt_t *pOut = (MFlt_t*)matrix3;
	MFlt_t inA1, inA2, inB1, inB2, out1, out2;
	uint32_t numSamples;
	uint32_t blkCnt;

	numSamples = (uint32_t) d1 * d2;

	/* Loop unrolling */
	blkCnt = numSamples >> 2u;
	while(blkCnt > 0u)
	{
		inA1 = pIn1[0];
    inA2 = pIn1[1];
    inB1 = pIn2[0];
    inB2 = pIn2[1];
    out1 = inA1 + inB1;
    out2 = inA2 + inB2;
		pOut[0] = out1;
		pOut[1] = out2;

		inA1 = pIn1[2];
    inA2 = pIn1[3];
    inB1 = pIn2[2];
    inB2 = pIn2[3];
    out1 = inA1 + inB1;
    out2 = inA2 + inB2;
		pOut[2] = out1;
		pOut[3] = out2;

		pIn1 += 4u;
		pIn2 += 4u;
		pOut += 4u;
		blkCnt--;
	}

	/* If the numSamples is not a multiple of 4, compute any remaining output samples here.
	 **	No	loop	unrolling	is	used.	*/
	blkCnt = numSamples % 0x4u;
	while(blkCnt > 0u)
	{
		*pOut++ = (*pIn1++) + (*pIn2++);
		blkCnt--;
	}
}

void Matrix_Subtraction(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], MFlt_t (*matrix3)[d2])
{
	MFlt_t *pIn1 = (MFlt_t*)matrix1;
	MFlt_t *pIn2 = (MFlt_t*)matrix2;
	MFlt_t *pOut = (MFlt_t*)matrix3;
	MFlt_t inA1, inA2, inB1, inB2, out1, out2;
	uint32_t numSamples;
	uint32_t blkCnt;

	numSamples = (uint32_t) d1 * d2;

	/* Loop unrolling */
	blkCnt = numSamples >> 2u;
	while(blkCnt > 0u)
	{
		inA1 = pIn1[0];
    inA2 = pIn1[1];
    inB1 = pIn2[0];
    inB2 = pIn2[1];
    out1 = inA1 - inB1;
    out2 = inA2 - inB2;
		pOut[0] = out1;
		pOut[1] = out2;

		inA1 = pIn1[2];
    inA2 = pIn1[3];
    inB1 = pIn2[2];
    inB2 = pIn2[3];
    out1 = inA1 - inB1;
    out2 = inA2 - inB2;
		pOut[2] = out1;
		pOut[3] = out2;

		pIn1 += 4u;
		pIn2 += 4u;
		pOut += 4u;
		blkCnt--;
	}

	/* If the numSamples is not a multiple of 4, compute any remaining output samples here.
	 **	No	loop	unrolling	is	used.	*/
	blkCnt = numSamples % 0x4u;
	while(blkCnt > 0u)
	{
		*pOut++ = (*pIn1++) - (*pIn2++);
		blkCnt--;
	}
}

void Matrix_MultiplyScalar(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t scale, MFlt_t (*matrix2)[d2])
{
	MFlt_t *pIn = (MFlt_t*)matrix1;
	MFlt_t *pOut = (MFlt_t*)matrix2;
	MFlt_t in1, in2, in3, in4;
	MFlt_t out1, out2, out3, out4;
	uint32_t numSamples;
	uint32_t blkCnt;

	numSamples = (uint32_t) d1 * d2;

	/* Loop Unrolling */
	blkCnt = numSamples >> 2;
	while(blkCnt > 0u)
	{
		in1 = pIn[0];
		in2 = pIn[1];
		in3 = pIn[2];
		in4 = pIn[3];
		out1 = in1 * scale;
		out2 = in2 * scale;
		out3 = in3 * scale;
		out4 = in4 * scale;
		pOut[0] = out1;
		pOut[1] = out2;
		pOut[2] = out3;
		pOut[3] = out4;

		pIn += 4u;
		pOut += 4u;
		blkCnt--;
	}

	/* If the numSamples is not a multiple of 4, compute any remaining output samples here.
	 **	No	loop	unrolling	is	used.	*/
	blkCnt = numSamples % 0x4u;
	while(blkCnt > 0u)
	{
		*pOut++ = (*pIn++) * scale;
		blkCnt--;
	}
}

void Matrix_MultiplyElementWise(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], MFlt_t (*matrix3)[d2])
{
	MFlt_t *pIn1 = (MFlt_t*)matrix1;
	MFlt_t *pIn2 = (MFlt_t*)matrix2;
	MFlt_t *pOut = (MFlt_t*)matrix3;
	MFlt_t inA1, inA2, inB1, inB2, out1, out2;
	uint32_t numSamples;
	uint32_t blkCnt;

	numSamples = (uint32_t) d1 * d2;

	/* Loop unrolling */
	blkCnt = numSamples >> 2u;
	while(blkCnt > 0u)
	{
		inA1 = pIn1[0];
		inB1 = pIn2[0];
		inA2 = pIn1[1];
		out1 = inA1 * inB1;
		inB2 = pIn2[1];
		inA1 = pIn1[2];
		out2 = inA2 * inB2;
		inB1 = pIn2[2];
		pOut[0] = out1;
		pOut[1] = out2;

		inA2 = pIn1[3];
		inB2 = pIn2[3];
		out1 = inA1 * inB1;
		out2 = inA2 * inB2;
		pOut[2] = out1;
		pOut[3] = out2;

		pIn1 += 4u;
		pIn2 += 4u;
		pOut += 4u;
		blkCnt--;
	}

	/* If the numSamples is not a multiple of 4, compute any remaining output samples here.
	 **	No	loop	unrolling	is	used.	*/
	blkCnt = numSamples % 0x4u;
	while(blkCnt > 0u)
	{
		*pOut++ = (*pIn1++) * (*pIn2++);
		blkCnt--;
	}
}

void Matrix_Multiply(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3])
{
	int32_t	i, k, p;
	uint32_t dimRows = d1;
	uint32_t dimCols = d3;
  uint32_t dimCommon = d2;

	MFlt_t *pIn1 = (MFlt_t*)matrix1;
	MFlt_t *pIn2 = (MFlt_t*)matrix2;
	MFlt_t *pInA = (MFlt_t*)matrix1;
	MFlt_t *pInB = (MFlt_t*)matrix2;
	MFlt_t *pOut = (MFlt_t*)matrix3;
	MFlt_t in1, in2, in3, in4;
	MFlt_t sum;

	uint32_t repCnt1 = dimCommon >> 2u;
  uint32_t repCnt2 = dimCommon % 0x4u;
	for(k = dimRows; k > 0u; k--)
	{
		for(i = 0u; i < dimCols; i++)
		{
			pIn1 = pInA; pIn2 = &pInB[i];
			sum = C(0.0);
			for(p = repCnt1; p > 0u; p--)
			{
				in1 = *pIn1++;
        in2 = *pIn1++;
				in3 = *pIn2;
				pIn2 += dimCols;
				in4 = *pIn2;
				pIn2 += dimCols;
        sum += in1 * in3;
        sum += in2 * in4;

				in1 = *pIn1++;
        in2 = *pIn1++;
				in3 = *pIn2;
				pIn2 += dimCols;
				in4 = *pIn2;
				pIn2 += dimCols;
        sum += in1 * in3;
        sum += in2 * in4;
			}
			for(p = repCnt2; p > 0u; p--)
			{
				sum += *pIn1++ * *pIn2;
        pIn2 += dimCols;
			}
      pOut[i] = sum;
		}
    pOut += dimCols; pInA += dimCommon;
	}
}

void Matrix_Transpose(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d1])
{
	int32_t	k, p;
	uint32_t dimRows = d2;
	uint32_t dimCols = d1;

	MFlt_t *pIn = (MFlt_t*)matrix1;
  MFlt_t *pInA = (MFlt_t*)matrix1;
	MFlt_t *pOut = (MFlt_t*)matrix2;

	uint32_t repCnt1 = dimCols >> 2u;
  uint32_t repCnt2 = dimCols % 0x4u;
	for(k = 0u; k < dimRows; k++)
	{
		pIn = &pInA[k];
		for(p = repCnt1; p > 0u; p--)
		{
			*pOut++ = *pIn;
			pIn += dimRows;
			*pOut++ = *pIn;
			pIn += dimRows;
			*pOut++ = *pIn;
			pIn += dimRows;
			*pOut++ = *pIn;
			pIn += dimRows;
		}
		for(p = repCnt2; p > 0u; p--)
		{
			*pOut++ = *pIn;
			pIn += dimRows;
		}
	}
}
