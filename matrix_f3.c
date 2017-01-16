#include "matrix.h"
#include <math.h>

#pragma GCC push_options
#pragma GCC optimize ("-O3")

uint32_t Matrix_SolveLowerTriangular(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3])
{
	int32_t	i, k, p;
	uint32_t dimRows = d2;
	uint32_t dimCols = d3;
	uint32_t dimCommon = d1;
	uint32_t dimMin = (dimCommon < dimRows)? dimCommon: dimRows;

	MFlt_t *pIn1 = (MFlt_t*)matrix1;
	MFlt_t *pIn2 = (MFlt_t*)matrix2;
	MFlt_t *pOut = (MFlt_t*)matrix3;
	MFlt_t *pOutA = (MFlt_t*)matrix3;
	MFlt_t *px;

	MFlt_t sum, idiv;

	for(k = 0u; k < dimMin; k++)
	{
		if(pIn1[k] == C(0.0)) return 1;
		idiv = C(1.0) / pIn1[k];
		for(i = 0u; i < dimCols; i++)
		{
			px = &pOutA[i];
			for(p = 0u, sum = C(0.0); p < k; p++)
			{
				sum += pIn1[p] * (*px);
				px += dimCols;
			}
			pOut[i] = (pIn2[i] - sum) * idiv;
		}
		pIn1 += dimRows; pIn2 += dimCols; pOut += dimCols;
	}
	return 0;
}

uint32_t Matrix_SolveUpperTriangular(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3])
{
	int32_t	i, k, p;
	uint32_t dimRows = d2;
	uint32_t dimCols = d3;
	uint32_t dimCommon = d1;
	uint32_t dimMin = (dimCommon < dimRows)? dimCommon: dimRows;

	MFlt_t *pIn1 = (MFlt_t*)matrix1[dimMin-1];
	MFlt_t *pIn2 = (MFlt_t*)matrix2[dimMin-1];
	MFlt_t *pOut = (MFlt_t*)matrix3[dimMin-1];
	MFlt_t *pOutA = (MFlt_t*)matrix3[dimMin-1];
	MFlt_t *px;

	MFlt_t sum, idiv;

	for(k = dimMin-1; k >= 0; k--)
	{
		if(pIn1[k] == C(0.0)) return 1;
		idiv = C(1.0) / pIn1[k];
		for(i = dimCols-1; i >= 0; i--)
		{
			px = &pOutA[i];
			for(p = dimMin-1, sum = C(0.0); p > k; p--)
			{
				sum += pIn1[p] * (*px);
				px -= dimCols;
			}
			pOut[i] = (pIn2[i] - sum) * idiv;
		}
		pIn1 -= dimRows; pIn2 -= dimCols; pOut -= dimCols;
	}
	return 0;
}

uint32_t Matrix_SolveUnitLowerTriangular(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3])
{
	int32_t	i, k, p;
	uint32_t dimRows = d2;
	uint32_t dimCols = d3;
	uint32_t dimCommon = d1;
	uint32_t dimMin = (dimCommon < dimRows)? dimCommon: dimRows;

	MFlt_t *pIn1 = (MFlt_t*)matrix1;
	MFlt_t *pIn2 = (MFlt_t*)matrix2;
	MFlt_t *pOut = (MFlt_t*)matrix3;
	MFlt_t *pOutA = (MFlt_t*)matrix3;
	MFlt_t *px;

	MFlt_t sum;

	for(k = 0u; k < dimMin; k++)
	{
		if(pIn1[k] == C(0.0)) return 1;
		for(i = 0u; i < dimCols; i++)
		{
			px = &pOutA[i];
			for(p = 0u, sum = C(0.0); p < k; p++)
			{
				sum += pIn1[p] * (*px);
				px += dimCols;
			}
			pOut[i] = (pIn2[i] - sum);
		}
		pIn1 += dimRows; pIn2 += dimCols; pOut += dimCols;
	}
	return 0;
}

uint32_t Matrix_SolveUnitUpperTriangular(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3])
{
	int32_t	i, k, p;
	uint32_t dimRows = d2;
	uint32_t dimCols = d3;
	uint32_t dimCommon = d1;
	uint32_t dimMin = (dimCommon < dimRows)? dimCommon: dimRows;

	MFlt_t *pIn1 = (MFlt_t*)matrix1[dimMin-1];
	MFlt_t *pIn2 = (MFlt_t*)matrix2[dimMin-1];
	MFlt_t *pOut = (MFlt_t*)matrix3[dimMin-1];
	MFlt_t *pOutA = (MFlt_t*)matrix3[dimMin-1];
	MFlt_t *px;

	MFlt_t sum;

	for(k = dimMin-1; k >= 0; k--)
	{
		if(pIn1[k] == C(0.0)) return 1;
		for(i = dimCols-1; i >= 0; i--)
		{
			px = &pOutA[i];
			for(p = dimMin-1, sum = C(0.0); p > k; p--)
			{
				sum += pIn1[p] * (*px);
				px -= dimCols;
			}
			pOut[i] = (pIn2[i] - sum);
		}
		pIn1 -= dimRows; pIn2 -= dimCols; pOut -= dimCols;
	}
	return 0;
}

uint32_t Matrix_InvertLowerTriangular(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d1])
{
	int32_t	i, k, p;
	uint32_t dimRows = d2;
	uint32_t dimCols = d1;
	uint32_t dimMin = (dimCols < dimRows)? dimCols: dimRows;

	MFlt_t *pIn = (MFlt_t*)matrix1;
	MFlt_t *pOut = (MFlt_t*)matrix2;
	MFlt_t *pOutA = (MFlt_t*)matrix2;
	MFlt_t *pOutB = (MFlt_t*)matrix2;
	MFlt_t *px;

	MFlt_t sum, idiv;

	for(k = 0u; k < dimMin; k++)
	{
		if(pIn[k] == C(0.0)) return 1;
		pOut[k] = C(1.0) / pIn[k];
		idiv = pOut[k];
		pOutB = pOutA;
		for(i = 0u; i < k; i++)
		{
			px = &pOutB[i];
			for(p = i, sum = C(0.0); p < k; p++)
			{
				sum += pIn[p] * (*px);
				px += dimCols;
			}
			pOut[i] = -sum * idiv;
			pOutB += dimCols;
		}
		pIn += dimRows; pOut += dimCols;
	}
	return 0;
}

uint32_t Matrix_InvertUpperTriangular(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d1])
{
	int32_t	i, k, p;
	uint32_t dimRows = d2;
	uint32_t dimCols = d1;
	uint32_t dimMin = (dimCols < dimRows)? dimCols: dimRows;

	MFlt_t *pIn = (MFlt_t*)matrix1[dimMin-1];
	MFlt_t *pOut = (MFlt_t*)matrix2[dimMin-1];
	MFlt_t *pOutA = (MFlt_t*)matrix2[dimMin-1];
	MFlt_t *pOutB = (MFlt_t*)matrix2[dimMin-1];
	MFlt_t *px;
	
	MFlt_t sum, idiv;

	for(k = dimMin-1; k >= 0; k--)
	{
		if(pIn[k] == C(0.0)) return 1;
		pOut[k] = C(1.0) / pIn[k];
		idiv = pOut[k];
		pOutB = pOutA;
		for(i = dimMin-1; i > k; i--)
		{
			px = &pOutB[i];
			for(p = i, sum = C(0.0); p > k; p--)
			{
				sum += pIn[p] * (*px);
				px -= dimCols;
			}
			pOut[i] = -sum * idiv;
			pOutB -= dimCols;
		}
		pIn -= dimRows; pOut -= dimCols;
	}
	return 0;
}

uint32_t Matrix_CroutLU_Decomposition(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d1], MFlt_t (*matrix3)[d2])
{
	int32_t	i, j, k, p;
	uint32_t dimRows = d1;
  uint32_t dimCols = d2;
	uint32_t dimMin = (dimCols < dimRows)? dimCols: dimRows;

	MFlt_t *pIn = (MFlt_t*)matrix1;
	MFlt_t *pInA = (MFlt_t*)matrix1;
	MFlt_t *pOut1 = (MFlt_t*)matrix2;
	MFlt_t *pOut2 = (MFlt_t*)matrix3;
	MFlt_t *pOutA1 = (MFlt_t*)matrix2;
	MFlt_t *pOutA2 = (MFlt_t*)matrix3;
	MFlt_t *pU = (MFlt_t*)matrix3;
	MFlt_t *px;

	MFlt_t sum, idiv;

	for(k = 0u; k < dimMin; k++)
	{
		pIn = pInA; pOut1 = pOutA1;
		for(i = k; i < dimRows; i++)
		{
			px = &pU[k];
			for(p = 0u, sum = pIn[k]; p < k; p++)
			{
				sum -= pOut1[p] * px[0];
				px += dimCols;
			}
			pOut1[k] = (sum);
			pIn += dimCols; pOut1 += dimRows;
		}

		pIn = pInA; pOut1 = pOutA1; pOut2 = pOutA2;
		if(pOut1[k] == C(0.0)) return 1;
		idiv = C(1.0) / pOut1[k];
		for(j = k; j < dimCols; j++)
		{
			px = &pU[j];
			for(p = 0u, sum = pIn[j]; p < k; p++)
			{
				sum -= pOut1[p] * px[0];
				px += dimCols;
			}
			pOut2[j] = (sum) * idiv;
		}
		pInA += dimCols; pOutA1 += dimRows; pOutA2 += dimCols;
	}
	return	0;
}

uint32_t Matrix_CholeskiLLT_Decomposition(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2])
{
	int32_t	i, k, p;
	uint32_t dimRows = d1;
  uint32_t dimCols = d2;

	MFlt_t *pIn = (MFlt_t*)matrix1;
	MFlt_t *pInA = (MFlt_t*)matrix1;
	MFlt_t *pOut = (MFlt_t*)matrix2;
	MFlt_t *pOutA = (MFlt_t*)matrix2;

	MFlt_t sum, idiv;

	for(k = 0u; k < dimRows; k++)
	{
		pIn = pInA; pOut = pOutA;
		for(p = 0u, sum = pIn[k]; p < k; p++)
		{
			sum -= pOut[p] * pOut[p];
		}
		if(sum <= C(0.0)) return 1;
		pOut[k] = C(sqrt)(sum);

		idiv = C(1.0) / pOut[k];
		for(i = k+1; i < dimRows; i++)
		{
			pIn += dimCols; pOut += dimCols;
			for(p = 0u, sum = pIn[k]; p < k; p++)
			{
				sum -= pOut[p] * pOutA[p];
			}
			pOut[k] = (sum) * idiv;
		}
		pInA += dimCols; pOutA += dimCols;
	}
	return	0;
}

uint32_t Matrix_HouseholderQTR_Decomposition(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d1], MFlt_t (*matrix3)[d2])
{
	int32_t	i, k, p;
	uint32_t dimRows = d1;
  uint32_t dimCols = d2;
	uint32_t dimMin = (dimRows-1u < dimCols)? dimRows-1u: dimCols;

	MFlt_t *pIn = (MFlt_t*)matrix1;
	MFlt_t *pOut1 = (MFlt_t*)matrix2;
	MFlt_t *pOut2 = (MFlt_t*)matrix3;
	MFlt_t *pOutA1 = (MFlt_t*)matrix2;
	MFlt_t *pOutA2 = (MFlt_t*)matrix3;
  MFlt_t *px;
	MFlt_t sum1, sum2, sum3, isum2, fact;
	MFlt_t	Vec[dimCols];
	uint32_t sign;

	uint32_t repCnt1, repCnt2;
	Matrix_Copy(d1, d2, (MFlt_t(*)[dimCols])pIn, (MFlt_t(*)[dimCols])pOut2);
	Matrix_Eye(d1, d1, (MFlt_t(*)[dimRows])pOut1);
	for(k = 0u; k < dimMin; k++)
	{
		pOut2 = &pOutA2[k], px = &Vec[k+1];
		sign =  signbit(*pOut2);
		repCnt1 = (dimRows-(k+1)) >> 2u;
    repCnt2 = (dimRows-(k+1)) % 0x4u;
    sum1 = C(0.0); sum3 = *pOut2;
		for(p = 0u; p < repCnt1; p++)
		{
      pOut2 += dimCols;
			*px++ = *pOut2;
      sum1 += (*pOut2) * (*pOut2);
      pOut2 += dimCols;
			*px++ = *pOut2;
      sum1 += (*pOut2) * (*pOut2);
      pOut2 += dimCols;
			*px++ = *pOut2;
      sum1 += (*pOut2) * (*pOut2);
      pOut2 += dimCols;
			*px++ = *pOut2;
      sum1 += (*pOut2) * (*pOut2);
		}
		for(p = 0u; p < repCnt2; p++)
		{
      pOut2 += dimCols;
			*px++ = *pOut2;
      sum1 += (*pOut2) * (*pOut2);
		}
		px = &Vec[k];
		sum2 = sum1;
    sum1 += sum3 * sum3;
    *px = sum3 + ((sign)? -C(sqrt)(sum1): C(sqrt)(sum1));
    sum2 += (*px) * (*px);
    isum2 = C(2.0) / sum2;

		repCnt1 = (dimRows-k) >> 2u;
		repCnt2 = (dimRows-k) % 0x4u;
		for(i = k; i < dimCols; i++)
		{
			pOut2 = &pOutA2[i]; px = &Vec[k];
			sum3 = C(0.0);
			for(p = 0u; p < repCnt1; p++)
			{
				sum3 += (*pOut2) * (*px++);
				pOut2 += dimCols;
				sum3 += (*pOut2) * (*px++);
				pOut2 += dimCols;
				sum3 += (*pOut2) * (*px++);
				pOut2 += dimCols;
				sum3 += (*pOut2) * (*px++);
				pOut2 += dimCols;
			}
			for(p = 0u; p < repCnt2; p++)
			{
				sum3 += (*pOut2) * (*px++);
				pOut2 += dimCols;
			}

			fact = sum3 * isum2;
      pOut2 = &pOutA2[i]; px = &Vec[k];
			for(p = 0u; p < repCnt1; p++)
			{
				*pOut2 = *pOut2 - fact * (*px++);
				pOut2 += dimCols;
				*pOut2 = *pOut2 - fact * (*px++);
				pOut2 += dimCols;
				*pOut2 = *pOut2 - fact * (*px++);
				pOut2 += dimCols;
				*pOut2 = *pOut2 - fact * (*px++);
				pOut2 += dimCols;
			}
			for(p = 0u; p < repCnt2; p++)
			{
				*pOut2 = *pOut2 - fact * (*px++);
				pOut2 += dimCols;
			}
		}
		for(i = 0u; i < dimRows; i++)
		{
			pOut1 = &pOutA1[i]; px = &Vec[k];
      sum3 = C(0.0);
			for(p = 0u; p < repCnt1; p++)
			{
				sum3 += (*px++) * (*pOut1);
				pOut1 += dimRows;
				sum3 += (*px++) * (*pOut1);
				pOut1 += dimRows;
				sum3 += (*px++) * (*pOut1);
				pOut1 += dimRows;
				sum3 += (*px++) * (*pOut1);
				pOut1 += dimRows;
			}
      for(p = 0u; p < repCnt2; p++)
			{
				sum3 += (*px++) * (*pOut1);
				pOut1 += dimRows;
			}

			fact = sum3 * isum2;
      pOut1 = &pOutA1[i]; px = &Vec[k];
      for(p = 0u; p < repCnt1; p++)
			{
				*pOut1 = *pOut1 - fact * (*px++);
				pOut1 += dimRows;
				*pOut1 = *pOut1 - fact * (*px++);
				pOut1 += dimRows;
				*pOut1 = *pOut1 - fact * (*px++);
				pOut1 += dimRows;
				*pOut1 = *pOut1 - fact * (*px++);
				pOut1 += dimRows;
			}
      for(p = 0u; p < repCnt2; p++)
			{
				*pOut1 = *pOut1 - fact * (*px++);
				pOut1 += dimRows;
			}
		}
    pOutA1 += dimRows; pOutA2 += dimCols;
	}
	return	0;
}
