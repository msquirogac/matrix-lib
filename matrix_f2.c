#include "matrix.h"

#pragma GCC push_options
#pragma GCC optimize ("-O3")

void Matrix_MultiplyAndAddition(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3], MFlt_t (*matrix4)[d3])
{
	int32_t	i, k, p;
	uint32_t dimRows = d1;
	uint32_t dimCols = d3;
	uint32_t dimCommon = d2;

	MFlt_t *pIn1 = (MFlt_t*)matrix1;
	MFlt_t *pIn2 = (MFlt_t*)matrix2;
	MFlt_t *pIn3 = (MFlt_t*)matrix3;
	MFlt_t *pInA = (MFlt_t*)matrix1;
	MFlt_t *pInB = (MFlt_t*)matrix2;
	MFlt_t *pOut = (MFlt_t*)matrix4;
	MFlt_t in1, in2, in3, in4;
	MFlt_t sum;

	uint32_t repCnt1 = dimCommon >> 2u;
	uint32_t repCnt2 = dimCommon % 0x4u;
	for(k = 0u; k < dimRows; k++)
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
			pOut[i] = sum + pIn3[i];
		}
		pOut += dimCols; pInA += dimCommon, pIn3 += dimCols;
	}
}

void Matrix_MultiplyAndSubtraction1(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3], MFlt_t (*matrix4)[d3])
{
	int32_t	i, k, p;
	uint32_t dimRows = d1;
	uint32_t dimCols = d3;
	uint32_t dimCommon = d2;

	MFlt_t *pIn1 = (MFlt_t*)matrix1;
	MFlt_t *pIn2 = (MFlt_t*)matrix2;
	MFlt_t *pIn3 = (MFlt_t*)matrix3;
	MFlt_t *pInA = (MFlt_t*)matrix1;
	MFlt_t *pInB = (MFlt_t*)matrix2;
	MFlt_t *pOut = (MFlt_t*)matrix4;
	MFlt_t in1, in2, in3, in4;
	MFlt_t sum;

	uint32_t repCnt1 = dimCommon >> 2u;
	uint32_t repCnt2 = dimCommon % 0x4u;
	for(k = 0u; k < dimRows; k++)
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
			pOut[i] = sum - pIn3[i];
		}
		pOut += dimCols; pInA += dimCommon, pIn3 += dimCols;
	}
}

void Matrix_MultiplyAndSubtraction2(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3], MFlt_t (*matrix4)[d3])
{
	int32_t	i, k, p;
	uint32_t dimRows = d1;
	uint32_t dimCols = d3;
	uint32_t dimCommon = d2;

	MFlt_t *pIn1 = (MFlt_t*)matrix1;
	MFlt_t *pIn2 = (MFlt_t*)matrix2;
	MFlt_t *pIn3 = (MFlt_t*)matrix3;
	MFlt_t *pInA = (MFlt_t*)matrix1;
	MFlt_t *pInB = (MFlt_t*)matrix2;
	MFlt_t *pOut = (MFlt_t*)matrix4;
	MFlt_t in1, in2, in3, in4;
	MFlt_t sum;

	uint32_t repCnt1 = dimCommon >> 2u;
	uint32_t repCnt2 = dimCommon % 0x4u;
	for(k = 0u; k < dimRows; k++)
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
			pOut[i] = pIn3[i] - sum;
		}
		pOut += dimCols; pInA += dimCommon, pIn3 += dimCols;
	}
}

void Matrix_MultiplyTranspose1(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], MFlt_t (*matrix3)[d3])
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
	MFlt_t sum;

	uint32_t repCnt1 = dimCommon >> 2u;
	uint32_t repCnt2 = dimCommon % 0x4u;
	for(k = 0u; k < dimRows; k++)
	{
		pIn2 = pInB;
		for(i = 0u; i < dimCols; i++)
		{
			pIn1 = pInA;
			sum = C(0.0);
			for(p = repCnt1; p > 0u; p--)
			{			
				sum += (*pIn1++) * (*pIn2++);
				sum += (*pIn1++) * (*pIn2++);
				sum += (*pIn1++) * (*pIn2++);
				sum += (*pIn1++) * (*pIn2++);			
			}
			for(p = repCnt2; p > 0u; p--)
			{
				sum += (*pIn1++) * (*pIn2++);
			}
			pOut[i] = sum;
		}
		pOut += dimCols; pInA += dimCommon;
	}
}

void Matrix_MultiplyTranspose2(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d1], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3])
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
	for(k = 0u; k < dimRows; k++)
	{
		for(i = 0u; i < dimCols; i++)
		{
			pIn1 = &pInA[k]; pIn2 = &pInB[i];
			sum = C(0.0);
			for(p = repCnt1; p > 0u; p--)
			{
				in1 = *pIn1;
				pIn1 += dimRows;
				in2 = *pIn1;
				pIn1 += dimRows;
				in3 = *pIn2;
				pIn2 += dimCols;
				in4 = *pIn2;
				pIn2 += dimCols;
				sum += in1 * in3;
				sum += in2 * in4;

				in1 = *pIn1;
				pIn1 += dimRows;
				in2 = *pIn1;
				pIn1 += dimRows;
				in3 = *pIn2;
				pIn2 += dimCols;
				in4 = *pIn2;
				pIn2 += dimCols;
				sum += in1 * in3;
				sum += in2 * in4;
			}
			for(p = repCnt2; p > 0u; p--)
			{
				sum += (*pIn1) * (*pIn2);
				pIn1 += dimRows;
				pIn2 += dimCols;
			}
			pOut[i] = sum;
		}
		pOut += dimCols;
	}
}

void Matrix_MultiplyTranspose3(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d1], MFlt_t (*matrix2)[d2], MFlt_t (*matrix3)[d3])
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
	for(k = 0u; k < dimRows; k++)
	{
		pIn2 = pInB;
		for(i = 0u; i < dimCols; i++)
		{
			pIn1 = pInA;
			sum = C(0.0);
			for(p = repCnt1; p > 0u; p--)
			{
				in1 = *pIn1;
				pIn1 += dimRows;
				in2 = *pIn1;
				pIn1 += dimRows;
				in3 = *pIn2++;
				in4 = *pIn2++;
				sum += in1 * in3;
				sum += in2 * in4;

				in1 = *pIn1;
				pIn1 += dimRows;
				in2 = *pIn1;
				pIn1 += dimRows;
				in3 = *pIn2++;
				in4 = *pIn2++;
				sum += in1 * in3;
				sum += in2 * in4;
			}
			for(p = repCnt2; p > 0u; p--)
			{
				sum += (*pIn1) * (*pIn2);
				pIn1 += dimRows;
				pIn2 ++;
			}
			pOut[i] = sum;
		}
		pOut += dimCols; pInA++;
	}
}

void Matrix_AdditionIdentity(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2])
{
	uint32_t dimRows = d1;
	uint32_t dimCols = d2;

	MFlt_t *pIn = (MFlt_t*)matrix1;
	MFlt_t *pOut = (MFlt_t*)matrix2;
	MFlt_t *pOutA = (MFlt_t*)matrix2;

	uint32_t dimColsNext = dimCols + 1u;
	uint32_t numElements = dimRows * dimCols;
	uint32_t repCnt1 = numElements >> 2u;
	uint32_t repCnt2 = numElements % 0x4u;
	while(repCnt1 > 0u)
	{
		*pOut++ = *pIn++;
		*pOut++ = *pIn++;
		*pOut++ = *pIn++;
		*pOut++ = *pIn++;
		repCnt1--;
	}
	while(repCnt2 > 0u)
	{
		*pOut++ = *pIn++;
		repCnt2--;
	}

	pOut = pOutA;
	repCnt1 = dimRows >> 2u;
	repCnt2 = dimRows % 0x4u;
	while(repCnt1 > 0u)
	{
		*pOut = *pOut + C(1.0);
		pOut += dimColsNext;
		*pOut = *pOut + C(1.0);
		pOut += dimColsNext;
		*pOut = *pOut + C(1.0);
		pOut += dimColsNext;
		*pOut = *pOut + C(1.0);
		pOut += dimColsNext;
		repCnt1--;
	}
	while(repCnt2 > 0u)
	{
		*pOut = *pOut + C(1.0);
		pOut += dimColsNext;
		repCnt2--;
	}
}

void Matrix_SubtractionIdentity1(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2])
{
	uint32_t dimRows = d1;
	uint32_t dimCols = d2;

	MFlt_t *pIn = (MFlt_t*)matrix1;
	MFlt_t *pOut = (MFlt_t*)matrix2;
	MFlt_t *pOutA = (MFlt_t*)matrix2;

	uint32_t dimColsNext = dimCols + 1u;
	uint32_t numElements = dimRows * dimCols;
	uint32_t repCnt1 = numElements >> 2u;
	uint32_t repCnt2 = numElements % 0x4u;
	while(repCnt1 > 0u)
	{
		*pOut++ = *pIn++;
		*pOut++ = *pIn++;
		*pOut++ = *pIn++;
		*pOut++ = *pIn++;
		repCnt1--;
	}
	while(repCnt2 > 0u)
	{
		*pOut++ = *pIn++;
		repCnt2--;
	}

	pOut = pOutA;
	repCnt1 = dimRows >> 2u;
	repCnt2 = dimRows % 0x4u;
	while(repCnt1 > 0u)
	{
		*pOut = *pOut - C(1.0);
		pOut += dimColsNext;
		*pOut = *pOut - C(1.0);
		pOut += dimColsNext;
		*pOut = *pOut - C(1.0);
		pOut += dimColsNext;
		*pOut = *pOut - C(1.0);
		pOut += dimColsNext;
		repCnt1--;
	}
	while(repCnt2 > 0u)
	{
		*pOut = *pOut - C(1.0);
		pOut += dimColsNext;
		repCnt2--;
	}
}

void Matrix_SubtractionIdentity2(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2])
{
	uint32_t dimRows = d1;
	uint32_t dimCols = d2;

	MFlt_t *pIn = (MFlt_t*)matrix1;
	MFlt_t *pOut = (MFlt_t*)matrix2;
	MFlt_t *pOutA = (MFlt_t*)matrix2;

	uint32_t dimColsNext = dimCols + 1u;
	uint32_t numElements = dimRows * dimCols;
	uint32_t repCnt1 = numElements >> 2u;
	uint32_t repCnt2 = numElements % 0x4u;
	while(repCnt1 > 0u)
	{
		*pOut++ = -(*pIn++);
		*pOut++ = -(*pIn++);
		*pOut++ = -(*pIn++);
		*pOut++ = -(*pIn++);
		repCnt1--;
	}
	while(repCnt2 > 0u)
	{
		*pOut++ = -(*pIn++);
		repCnt2--;
	}

	pOut = pOutA;
	repCnt1 = dimRows >> 2u;
	repCnt2 = dimRows % 0x4u;
	while(repCnt1 > 0u)
	{
		*pOut = *pOut + C(1.0);
		pOut += dimColsNext;
		*pOut = *pOut + C(1.0);
		pOut += dimColsNext;
		*pOut = *pOut + C(1.0);
		pOut += dimColsNext;
		*pOut = *pOut + C(1.0);
		pOut += dimColsNext;
		repCnt1--;
	}
	while(repCnt2 > 0u)
	{
		*pOut = *pOut + C(1.0);
		pOut += dimColsNext;
		repCnt2--;
	}
}
