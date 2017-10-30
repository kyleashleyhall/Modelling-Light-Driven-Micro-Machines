/* File: linalg.h
 * $Date:: 2013-03-18 17:23:24 +0700 #$
 * Descr: definitions for linear algebra operations on large vectors; see source (linalg.c) for details
 *
 * Copyright (C) 2006,2008,2010-2013 ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#ifndef __linalg_h
#define __linalg_h

// project headers
#include "timing.h"   // for TIME_TYPE
#include "types.h"    // for doublecomplex

void nInit(doublecomplex * restrict a);
void nCopy(doublecomplex * restrict a,doublecomplex * restrict b);
double nNorm2(doublecomplex * restrict a,TIME_TYPE *comm_timing);
void nDotProd(doublecomplex * restrict a,doublecomplex * restrict b,doublecomplex c,TIME_TYPE *comm_timing);
void nDotProd_conj(doublecomplex * restrict a,doublecomplex * restrict b,doublecomplex c,TIME_TYPE *comm_timing);
void nDotProdSelf_conj(doublecomplex * restrict a,doublecomplex c,TIME_TYPE *comm_timing);
void nDotProdSelf_conj_Norm2(doublecomplex * restrict a,doublecomplex c,double * restrict norm,TIME_TYPE *comm_timing);
void nIncrem110_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,doublecomplex * restrict c,
	const doublecomplex c1,const doublecomplex c2);
void nIncrem011_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,doublecomplex * restrict c,
	const doublecomplex c1,const doublecomplex c2);
void nIncrem111_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,doublecomplex * restrict c,
	const doublecomplex c1,const doublecomplex c2,const doublecomplex c3);
void nIncrem(doublecomplex * restrict a,doublecomplex * restrict b,double * restrict inprod,TIME_TYPE *comm_timing);
void nDecrem(doublecomplex * restrict a,doublecomplex * restrict b,double * restrict inprod,TIME_TYPE *comm_timing);
void nIncrem01(doublecomplex * restrict a,doublecomplex * restrict b,const double c,double * restrict inprod,
	TIME_TYPE *comm_timing);
void nIncrem10(doublecomplex * restrict a,doublecomplex * restrict b,const double c,double * restrict inprod,
	TIME_TYPE *comm_timing);
void nIncrem11_d_c(doublecomplex * restrict a,doublecomplex * restrict b,const double c1,const doublecomplex c2,
	double * restrict inprod,TIME_TYPE *comm_timing);
void nIncrem01_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,const doublecomplex c,
	double * restrict inprod,TIME_TYPE *comm_timing);
void nIncrem110_d_c_conj(doublecomplex * restrict a,doublecomplex * restrict b,doublecomplex * restrict c,
	const double c1,const doublecomplex c2,double * restrict inprod,TIME_TYPE *comm_timing);
void nIncrem10_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,const doublecomplex c,
	double * restrict inprod,TIME_TYPE *comm_timing);
void nLinComb_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,doublecomplex * restrict c,
	const doublecomplex c1,const doublecomplex c2,double * restrict inprod,TIME_TYPE *comm_timing);
void nLinComb1_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,doublecomplex * restrict c,
	const doublecomplex c1,double * restrict inprod,TIME_TYPE *comm_timing);
void nLinComb1_cmplx_conj(doublecomplex * restrict a,doublecomplex * restrict b,doublecomplex * restrict c,
	const doublecomplex c1,double * restrict inprod,TIME_TYPE *comm_timing);
void nSubtr(doublecomplex * restrict a,doublecomplex * restrict b,doublecomplex * restrict c,
	double * restrict inprod,TIME_TYPE *comm_timing);
void nMult(doublecomplex * restrict a,doublecomplex * restrict b,const double c);
void nMult_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,const doublecomplex c);
void nMultSelf(doublecomplex * restrict a,const double c);
void nMultSelf_conj(doublecomplex * restrict a,const double c);
void nMultSelf_cmplx(doublecomplex * restrict a,const doublecomplex c);
void nMult_mat(doublecomplex * restrict a,doublecomplex * restrict b,doublecomplex (* restrict c)[3]);
void nMultSelf_mat(doublecomplex * restrict a,doublecomplex (* restrict c)[3]);
void nConj(doublecomplex * restrict a);

#endif // __linalg_h
