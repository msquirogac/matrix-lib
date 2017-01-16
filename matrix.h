#include <stdint.h>

#ifdef MATRIX_USE_DOUBLE
  #define C(x) x
	typedef double MFlt_t;
#else
	#define MATRIX_USE_SINGLE
  #define C(x) x##f
	typedef float MFlt_t;
#endif

#define Matrix_Add					Matrix_Addition
#define Matrix_Sub					Matrix_Subtraction
#define Matrix_Scale				Matrix_MultiplyScalar
#define Matrix_Mul					Matrix_Multiply

#define Matrix_Trans				Matrix_Transpose
#define Matrix_Inv					Matrix_Invert1
#define Matrix_Inv1					Matrix_Invert1
#define Matrix_Inv2					Matrix_Invert2
#define Matrix_Inv3					Matrix_Invert3

#define Matrix_SolLT				Matrix_SolveLowerTriangular
#define Matrix_SolUT				Matrix_SolveUpperTriangular
#define Matrix_SolULT				Matrix_SolveUnitLowerTriangular
#define Matrix_SolUUT				Matrix_SolveUnitUpperTriangular
#define Matrix_InvLT				Matrix_InvertLowerTriangular
#define Matrix_InvUT				Matrix_InvertUpperTriangular

#define Matrix_AddI					Matrix_AdditionIdentity
#define Matrix_SubI1				Matrix_SubtractionIdentity1
#define Matrix_SubI2				Matrix_SubtractionIdentity2

#define Matrix_MulTrans1		Matrix_MultiplyTranspose1
#define Matrix_MulTrans2		Matrix_MultiplyTranspose2
#define Matrix_MulTrans3		Matrix_MultiplyTranspose3

#define Matrix_MulAdd				Matrix_MultiplyAndAddition
#define Matrix_MulSub1			Matrix_MultiplyAndSubtraction1
#define Matrix_MulSub2			Matrix_MultiplyAndSubtraction2

#define Matrix_MulElemWise	Matrix_MultiplyElementWise


uint32_t Matrix_AssignMemory(uint16_t d1, uint16_t d2, MFlt_t (**matrix1)[d2]);
uint32_t Matrix_AssignClearMemory(uint16_t d1, uint16_t d2, MFlt_t (**matrix1)[d2]);
void Matrix_FreeMemory(uint16_t d1, uint16_t d2, MFlt_t (**matrix1)[d2]);
void Matrix_Copy(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2]);
uint32_t Matrix_Compare(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], MFlt_t eps);
void Matrix_GetRow(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], uint16_t row);
void Matrix_SetRow(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], uint16_t row);
void Matrix_GetColumn(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], uint16_t col);
void Matrix_SetColumn(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], uint16_t col);

void Matrix_Eye(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2]);
void Matrix_Zeros(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2]);

void Matrix_Addition(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], MFlt_t (*matrix3)[d2]);
void Matrix_Subtraction(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], MFlt_t (*matrix3)[d2]);
void Matrix_MultiplyScalar(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t scale, MFlt_t (*matrix2)[d2]);
void Matrix_MultiplyElementWise(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], MFlt_t (*matrix3)[d2]);
void Matrix_Multiply(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3]);

void Matrix_Transpose(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d1]);
void Matrix_Invert1(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d1]);
void Matrix_Invert2(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d1]);
void Matrix_Invert3(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d1]);

uint32_t Matrix_SolveLowerTriangular(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3]);
uint32_t Matrix_SolveUpperTriangular(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3]);
uint32_t Matrix_SolveUnitLowerTriangular(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3]);
uint32_t Matrix_SolveUnitUpperTriangular(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3]);
uint32_t Matrix_InvertLowerTriangular(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d1]);
uint32_t Matrix_InvertUpperTriangular(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d1]);
uint32_t Matrix_CroutLU_Decomposition(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], MFlt_t (*matrix3)[d2]);
uint32_t Matrix_CholeskiLLT_Decomposition(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2]);
uint32_t Matrix_HouseholderQTR_Decomposition(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d1], MFlt_t (*matrix3)[d2]);

void Matrix_MultiplyAndAddition(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3], MFlt_t (*matrix4)[d3]);
void Matrix_MultiplyAndSubtraction1(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3], MFlt_t (*matrix4)[d3]);
void Matrix_MultiplyAndSubtraction2(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3], MFlt_t (*matrix4)[d3]);
void Matrix_MultiplyTranspose1(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2], MFlt_t (*matrix3)[d3]);
void Matrix_MultiplyTranspose2(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d1], MFlt_t (*matrix2)[d3], MFlt_t (*matrix3)[d3]);
void Matrix_MultiplyTranspose3(uint16_t d1, uint16_t d2, uint16_t d3, MFlt_t (*matrix1)[d1], MFlt_t (*matrix2)[d2], MFlt_t (*matrix3)[d3]);
void Matrix_AdditionIdentity(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2]);
void Matrix_SubtractionIdentity1(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2]);
void Matrix_SubtractionIdentity2(uint16_t d1, uint16_t d2, MFlt_t (*matrix1)[d2], MFlt_t (*matrix2)[d2]);
