/* File: timing.c
 * $Date:: 2013-05-04 12:16:05 +0700 #$
 * Descr: basic timing and statistics routines
 *
 * Copyright (C) 2006,2008-2013 ADDA contributors
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
#include "const.h" // keep this first
#include "timing.h" // corresponding header
// project headers
#include "comm.h"
#include "io.h"
#include "vars.h"
// system headers
#include <time.h>
#include <stdio.h>

#ifdef ADDA_MPI
#	define TO_SEC(p) (p)
#else
#	define TO_SEC(p) ((p) / (double) CLOCKS_PER_SEC)
#endif

// SEMI-GLOBAL VARIABLES

// used in CalculateE.c
TIME_TYPE Timing_EPlane,Timing_EPlaneComm,    // for Eplane calculation: total and comm
          Timing_IntField,Timing_IntFieldOne, // for internal fields: total & one calculation
          Timing_IncBeam,                     // for generation (or reading) and saving (if needed) of incident beam
          Timing_ScatQuan;                    // for integral scattering quantities
size_t TotalEFieldPlane; // total number of planes for scattered field calculations
// used in calculator.c
TIME_TYPE Timing_Init; // for total initialization of the program (before CalculateE)
size_t TotalEval;      // total number of orientation evaluations
#ifdef OPENCL
TIME_TYPE Timing_OCL_Init; // for initialization of OpenCL (including building program)
#endif
// used in comm.c
TIME_TYPE Timing_InitDmComm; // communication time for initialization of D-matrix
// used in crosssec.c
	// total time for all_dir and scat_grid calculations
TIME_TYPE Timing_EFieldAD,Timing_EFieldADComm,  // time for all_dir: total & comm
          Timing_EFieldSG,Timing_EFieldSGComm,  // time for scat_dir: total & comm
          Timing_ScatQuanComm;                  // time for comm of scat.quantities
// used in fft.c
TIME_TYPE Timing_FFT_Init, // for initialization of FFT routines
          Timing_Dm_Init;  // for building Dmatrix
// used in iterative.c
TIME_TYPE Timing_OneIter,Timing_OneIterComm,    // for one iteration: total & comm
          Timing_InitIter,Timing_InitIterComm,  // for initialization of iterations: total & comm
          Timing_IntFieldOneComm;               // comm for one calculation of the internal fields
size_t TotalIter;                               // total number of iterations performed
// used in make_particle.c
TIME_TYPE Timing_Particle,                 // for particle construction
          Timing_Granul,Timing_GranulComm; // for granule generation: total & comm
// used in matvec.c
size_t TotalMatVec; // total number of matrix-vector products

#define FFORMT "%.4f" // format for timing results

//======================================================================================================================

void StartTime(void)
// start global time
{
	time(&wt_start);
	last_chp_wt=wt_start;
#ifndef ADDA_MPI  // otherwise this initialization is performed immediately after MPI_Init
	tstart_main = GET_TIME();
#endif
}

//======================================================================================================================

void InitTiming(void)
// init timing variables and counters
{
	TotalIter=TotalMatVec=TotalEval=TotalEFieldPlane=0;
	Timing_EField=Timing_FileIO=Timing_IntField=Timing_ScatQuan=Timing_Integration=0;
	Timing_ScatQuanComm=Timing_InitDmComm=0;
#ifdef SPARSE
	Timing_Dm_Init=Timing_Granul=Timing_FFT_Init=Timing_GranulComm=0;
#endif	
}

//======================================================================================================================

void FinalStatistics(void)
// print final output and statistics
{
	time_t wt_end;
	TIME_TYPE Timing_TotalTime;

	// wait for all processes to show correct execution time
	Synchronize();
	if (IFROOT) {
		// last time measurements
		Timing_TotalTime = GET_TIME() - tstart_main;
		time(&wt_end);
		// log statistics
		fprintf(logfile,
			"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
			"                Timing Results             \n"
			"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
		if (!prognosis) {
			if (orient_avg) fprintf(logfile,
				"Total number of single particle evaluations: %zu\n",TotalEval);
			fprintf(logfile,
				"Total number of iterations: %zu\n"
				"Total number of matrix-vector products: %zu\n"
				"Total planes of E field calculation (each %d points): %zu\n\n",
				TotalIter,TotalMatVec,nTheta,TotalEFieldPlane);
		}
		fprintf(logfile,
			"Total wall time:     %.0f\n",difftime(wt_end,wt_start));
		fprintf(logfile,
#ifdef ADDA_MPI
			"--Everything below is also wall times--\n"
			"Time since MPI_Init: "FFORMT"\n",
#else
			"--Everything below is processor times--\n"
			"Total time:          "FFORMT"\n",
#endif
			TO_SEC(Timing_TotalTime));
		fprintf(logfile,
			"  Initialization time: "FFORMT"\n",TO_SEC(Timing_Init));
		if (!prognosis) {
#ifdef OPENCL
			fprintf(logfile,
				"    init OpenCL          "FFORMT"\n",TO_SEC(Timing_OCL_Init));

#endif
#ifndef SPARSE
			fprintf(logfile,
				"    init Dmatrix         "FFORMT"\n",TO_SEC(Timing_Dm_Init));
#	ifdef PARALLEL
			fprintf(logfile,
				"      communication:       "FFORMT"\n",TO_SEC(Timing_InitDmComm));
#	endif
			fprintf(logfile,
				"    FFT setup:           "FFORMT"\n",TO_SEC(Timing_FFT_Init));
#endif // !SPARSE
		}
		fprintf(logfile,
			"    make particle:       "FFORMT"\n",TO_SEC(Timing_Particle));
		if (sh_granul) {
			fprintf(logfile,
				"      granule generator:   "FFORMT"\n",TO_SEC(Timing_Granul));
#ifdef PARALLEL
			fprintf(logfile,
				"        communication:       "FFORMT"\n",TO_SEC(Timing_GranulComm));
#endif
		}
		if (!prognosis) {
			fprintf(logfile,
				"  Internal fields:     "FFORMT"\n"
				"    one solution:        "FFORMT"\n",
				TO_SEC(Timing_IntField),TO_SEC(Timing_IntFieldOne));
#ifdef PARALLEL
			fprintf(logfile,
				"      communication:       "FFORMT"\n",TO_SEC(Timing_IntFieldOneComm));
#endif
			fprintf(logfile,
				"      incident beam:       "FFORMT"\n",TO_SEC(Timing_IncBeam));
			fprintf(logfile,
				"      init solver:         "FFORMT"\n",TO_SEC(Timing_InitIter));
#ifdef PARALLEL
			fprintf(logfile,
				"        communication:       "FFORMT"\n",TO_SEC(Timing_InitIterComm));
#endif
			fprintf(logfile,
				"      one iteration:       "FFORMT"\n",TO_SEC(Timing_OneIter));
#ifdef PARALLEL
			fprintf(logfile,
				"        communication:       "FFORMT"\n",TO_SEC(Timing_OneIterComm));
#endif
			fprintf(logfile,
				"  Scattered fields:    "FFORMT"\n",TO_SEC(Timing_EField));
			if (yzplane) {
				fprintf(logfile,
					"    one plane:           "FFORMT"\n",TO_SEC(Timing_EPlane));
#ifdef PARALLEL
				fprintf(logfile,
					"      communication:       "FFORMT"\n",TO_SEC(Timing_EPlaneComm));
#endif
			}
			if (all_dir) {
				fprintf(logfile,
					"    one alldir:          "FFORMT"\n",TO_SEC(Timing_EFieldAD));
#ifdef PARALLEL
				fprintf(logfile,
					"      communication:       "FFORMT"\n",TO_SEC(Timing_EFieldADComm));
#endif
			}
			if (scat_grid) {
				fprintf(logfile,
					"    one scat_grid:          "FFORMT"\n",TO_SEC(Timing_EFieldSG));
#ifdef PARALLEL
				fprintf(logfile,
					"      communication:       "FFORMT"\n",TO_SEC(Timing_EFieldSGComm));
#endif
			}
			fprintf (logfile,
				"  Other sc.quantities: "FFORMT"\n",TO_SEC(Timing_ScatQuan));
#ifdef PARALLEL
			fprintf(logfile,
				"    communication:       "FFORMT"\n",TO_SEC(Timing_ScatQuanComm));
#endif
		}
		fprintf (logfile,
				"File I/O:            "FFORMT"\n",TO_SEC(Timing_FileIO));
		if (!prognosis) fprintf (logfile,
				"Integration:         "FFORMT"\n",TO_SEC(Timing_Integration));
		// close logfile
		FCloseErr(logfile,F_LOG,ONE_POS);
	}
}
