/* File: oclkernels.cl
 * $Date:: 2013-03-18 17:23:24 +0700 #$
 * Descr: Kernel File for OpenCL kernels. Includes all subfunctions of the Matvec routine as OpenCL kernels
 *
 * Copyright (C) 2010-2013 ADDA contributors
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
// Somehow current AMD drivers do not define CL_VERSION_1_0 when compiling OpenCL kernels
//#ifndef CL_VERSION_1_0
//#	error "OpenCL version at least 1.0 is required"
//#endif

#ifdef USE_DOUBLE
#	ifdef DOUBLE_AMD
#		pragma OPENCL EXTENSION cl_amd_fp64 : enable
#	else
#		pragma OPENCL EXTENSION cl_khr_fp64 : enable
#	endif
#endif

// defines type to work with variables defined as size_t in the main ADDA code
#ifdef SIZET_UINT
typedef uint in_sizet;
#elif defined(SIZET_ULONG)
typedef ulong in_sizet;
#else
#	error "size_t alternative is not defined"
#endif

//======================================================================================================================
// functions used in kernels

void cMult(__constant double2 *a,__global const double2 *b,__global double2 *c)
// complex multiplication; c=ab; !!! c should be different from a and b !!!
{
	(*c).s0=(*a).s0*(*b).s0 - (*a).s1*(*b).s1;
	(*c).s1=(*a).s1*(*b).s0 + (*a).s0*(*b).s1;
}

//======================================================================================================================

void cMult2(__constant double2 *a,__global const double2 *b,double2 *c)
// complex multiplication; c=ab; !!! c should be different from a and b !!! c is private
{
	(*c).s0=(*a).s0*(*b).s0 - (*a).s1*(*b).s1;
	(*c).s1=(*a).s1*(*b).s0 + (*a).s0*(*b).s1;
}

//======================================================================================================================

double cvNorm2(__global const double2 *a)
// square of the norm of a complex vector[3]
{
	return ( a[0].s0*a[0].s0 + a[0].s1*a[0].s1 + a[1].s0*a[1].s0 + a[1].s1*a[1].s1
	       + a[2].s0*a[2].s0 + a[2].s1*a[2].s1 );
}

//======================================================================================================================
// Arith1 kernels

__kernel void nConj(__global double2 *a)
{
	const size_t id=get_global_id(0);

	a[id].s1= -a[id].s1;
}

//======================================================================================================================

__kernel void clzero(__global double2 *input)
{
	const size_t id = get_global_id(0);

	input[id] = 0.0;
}

//======================================================================================================================

__kernel void arith1(__global const uchar *material,__global const ushort *position,__constant double2 *cc_sqrt,
	__global const double2 *argvec, __global double2 *Xmatrix,const in_sizet local_Nsmall,const in_sizet smallY,
	const in_sizet gridX)
{
	const size_t id=get_global_id(0);
	const size_t j=3*id;
	const uchar mat=material[id];
	size_t index;
	int xcomp;

	index = ((position[j+2]*smallY+position[j+1])*gridX+position[j]);
	for (xcomp=0;xcomp<3;xcomp++) cMult(&cc_sqrt[mat*3+xcomp],&argvec[j+xcomp],&Xmatrix[index+xcomp*local_Nsmall]);
}

//======================================================================================================================
// Arith2 kernel

__kernel void arith2(__global const double2 *Xmatrix,__global double2 *slices,const in_sizet gridZ,
	const in_sizet smallY,const in_sizet gridX,const in_sizet gridYZ,const in_sizet local_Nsmall,const in_sizet x)
{
	const size_t y=get_global_id(0);
	const size_t z=get_global_id(1);
	size_t i;
	size_t j;
	int xcomp;

	i = y*gridZ+z;
	j = (z*smallY+y)*gridX+x;
	for (xcomp=0;xcomp<3;xcomp++) slices[i+xcomp*gridYZ]=Xmatrix[j+xcomp*local_Nsmall];
}

//======================================================================================================================
// Arith3 kernels and functions

void cSymMatrVec(const double2 *matr,const double2 *vec,double2 *res)
// multiplication of complex symmetric matrix[6] by complex vec[3]; res=matr.vec
{
	res[0].s0 = matr[0].s0*vec[0].s0 - matr[0].s1*vec[0].s1
	          + matr[1].s0*vec[1].s0 - matr[1].s1*vec[1].s1
	          + matr[2].s0*vec[2].s0 - matr[2].s1*vec[2].s1;
	res[0].s1 = matr[0].s0*vec[0].s1 + matr[0].s1*vec[0].s0
	          + matr[1].s0*vec[1].s1 + matr[1].s1*vec[1].s0
	          + matr[2].s0*vec[2].s1 + matr[2].s1*vec[2].s0;

	res[1].s0 = matr[1].s0*vec[0].s0 - matr[1].s1*vec[0].s1
	          + matr[3].s0*vec[1].s0 - matr[3].s1*vec[1].s1
	          + matr[4].s0*vec[2].s0 - matr[4].s1*vec[2].s1;
	res[1].s1 = matr[1].s0*vec[0].s1 + matr[1].s1*vec[0].s0
	          + matr[3].s0*vec[1].s1 + matr[3].s1*vec[1].s0
	          + matr[4].s0*vec[2].s1 + matr[4].s1*vec[2].s0;

	res[2].s0 = matr[2].s0*vec[0].s0 - matr[2].s1*vec[0].s1
	          + matr[4].s0*vec[1].s0 - matr[4].s1*vec[1].s1
	          + matr[5].s0*vec[2].s0 - matr[5].s1*vec[2].s1;
	res[2].s1 = matr[2].s0*vec[0].s1 + matr[2].s1*vec[0].s0
	          + matr[4].s0*vec[1].s1 + matr[4].s1*vec[1].s0
	          + matr[5].s0*vec[2].s1 + matr[5].s1*vec[2].s0;
}

//======================================================================================================================

__kernel void arith3(__global double2 *slices_tr,__global const double2 *Dmatrix,const in_sizet local_x0,
	const in_sizet smallY,const in_sizet smallZ,const in_sizet gridX,const in_sizet DsizeY,const in_sizet DsizeZ,
	const char NDCOMP,const char reduced_FFT,const char transposed,const in_sizet x)
{
	size_t const z = get_global_id(0);
	size_t const y = get_global_id(1);
	size_t const gridZ = get_global_size(0);
	size_t const gridY = get_global_size(1);
	double2 xv[3];
	double2 yv[3];
	double2 fmat[6];
	int Xcomp;
	const size_t i=z *gridY + y; // indexSliceZY
	size_t j;
	size_t xa=x;
	size_t ya=y;
	size_t za=z;

	// works, because of the double2 vector type
	for (Xcomp=0;Xcomp<3;Xcomp++) xv[Xcomp]=slices_tr[i+Xcomp*gridY*gridZ];

	if (transposed==1) { // used only for G_SO
		if (x>0) xa=gridX-x;
		if (y>0) ya=gridY-y;
		if (z>0) za=gridZ-z;
	}
	else {
		if (y>=DsizeY) ya=gridY-y;
		if (z>=DsizeZ) za=gridZ-z;
	}
	j=NDCOMP*(xa*DsizeY*DsizeZ+za*DsizeY+ya);

	for (int f=0;f<6;f++) fmat[f]=Dmatrix[j+f]; // maybe could be done more efficient TODO

	if (reduced_FFT==1) {
		if (y>smallY) {
			fmat[1]*=-1;
			if (z>smallZ) fmat[2]*=-1;
			else fmat[4]*=-1;
		}
		else if (z>smallZ) {
			fmat[2]*=-1;
			fmat[4]*=-1;
		}
	}
	cSymMatrVec(fmat,xv,yv); // yv=fmat*xv
	for (Xcomp=0;Xcomp<3;Xcomp++) slices_tr[i+Xcomp*gridY*gridZ]=yv[Xcomp];
}

//======================================================================================================================
// Arith4 kernel

__kernel void arith4(__global double2 *Xmatrix,__global const double2 *slices,const in_sizet gridZ,
	const in_sizet smallY,const in_sizet gridX,const in_sizet gridYZ,const in_sizet local_Nsmall,const in_sizet x)
{
	const size_t y =get_global_id(0);
	const size_t z =get_global_id(1);
	size_t i;
	size_t j;
	int xcomp;

	i = y*gridZ+z;
	j = (z*smallY+y)*gridX+x;
	for (xcomp=0;xcomp<3;xcomp++) Xmatrix[j+xcomp*local_Nsmall]=slices[i+xcomp*gridYZ];
}

//======================================================================================================================
// Arith5 kernel

__kernel void arith5(__global const uchar *material,__global const ushort *position,__constant double2 *cc_sqrt,
	__global const double2 *argvec,__global const double2 *Xmatrix,const in_sizet local_Nsmall,const in_sizet smallY,
	const in_sizet gridX,__global double2 *resultvec)
{
	const size_t id = get_global_id(0);
	const size_t j=3*id;
	const uchar mat = material[id];
	size_t index;
	double2 temp;
	int xcomp;

	index = ((position[j+2]*smallY+position[j+1])*gridX+position[j]);
	for (xcomp=0;xcomp<3;xcomp++) {
		cMult2(&cc_sqrt[mat*3+xcomp],&Xmatrix[index+xcomp*local_Nsmall],&temp);
		resultvec[j+xcomp]=argvec[j+xcomp]+temp;
	}
}

//======================================================================================================================

__kernel void inpr(__global double *inprod, __global const double2 *resultvec)
{
	const size_t id = get_global_id(0);

	inprod[id]=cvNorm2(resultvec+(id*3));
}

//======================================================================================================================
// transpose kernel

__kernel void transpose(__global const double2 *input,__global double2 *output,const in_sizet width,
	const in_sizet height)
{
	const size_t idz = get_global_id(0);
	const size_t idy = get_global_id(1);
	const size_t wth = width*height;

	for (int k=0;k<3;k++) output[idz*height+idy+k*wth]=input[idy*width+idz+k*wth];
}

//======================================================================================================================

// This corresponds to value of tblock in TransposeYZ() in fft.c
#define BLOCK_DIM 16
__kernel void transposeo(__global const double2 *idata,__global double2 *odata,const in_sizet width,
	const in_sizet height,__local double2 *block)
//optimised transpose kernel with cache and removed bank conflicts obtained from Nvidia SDK samples
{
	// read tiles into local memory
	size_t xIndex = get_global_id(0);
	size_t yIndex = get_global_id(1);
	const size_t zIndex = get_global_id(2);
	const size_t htw = height*width;
	if ((xIndex < width) && (yIndex < height)) {
		size_t index_in = yIndex * width + xIndex + htw * zIndex;
		block[get_local_id(1)*(BLOCK_DIM+1)+get_local_id(0)] = idata[index_in];
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	// write transposed tile back to global memory
	xIndex = get_group_id(1) * BLOCK_DIM + get_local_id(0);
	yIndex = get_group_id(0) * BLOCK_DIM + get_local_id(1);
	if ((xIndex < height) && (yIndex < width)) {
		size_t index_out = yIndex * height + xIndex + htw * zIndex;
		odata[index_out] = block[get_local_id(0)*(BLOCK_DIM+1)+get_local_id(1)];
	}
}
