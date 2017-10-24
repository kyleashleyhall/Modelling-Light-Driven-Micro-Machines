Notes for SPHRN_PBC example

                        Geometry

ddscat.par is set up to calculate absorption and scattering by a "target unit
cell" (TUC) composed of a cluater of spheres, with the periodic repetition
of the TUC in the y and z directions.

The arrangement of the spheres in the TUC is specified by the file
'BAM2.16.1.targ' giving the locations and sizes of 16 spheres
In BAM2.16.1.targ the spheres all have radii a = 0.5*D
where D is some as-yet unknown length.

The spheres do not overlap, hence the effective radius of the material in
one TUC is aeff = 16^{1/3}*(D/2)

Inspecting BAM1.16.1.targ one finds that the target extends in the
x-direction from 
(-0.433343-0.5)*D = -0.933343*D 
to 
(+1.356062+0.5)*D = 1.856062*D

ddscat.par specifies SHPAR1,SHPAR2,SHPAR3 = 32 50 5
SHPAR1 = maximum extent/d of the TUC in the TF x-direction
SHPAR2 = PYAEFF = periodicity in y direction / aeff
SHPAR3 = PZAEFF = periodicity in z direction / aeff

Thus we expect that 32*d = (1.856062+0.933343)D = 2.78941*D
Thus D = 32*d/2.78941 = 11.47198 d
     a = D/2 = 5.73599 d

Thus individual spheres are each represented by approximately

N_D = (4*pi/3)*a^3 / d^3 dipoles
N_D = (4*pi/3)*(a/d)^3 = (4*pi/3)*(5.73599)^3 = 791 dipoles

and the cluster of 16 spheres should contain approximately 
N = 16*791 = 12648 dipoles

The actual number of dipoles is (see ddscat.log) N=11188 .

ddscat.par specifies aeff = 0.1 unit.
aeff^3 = 16*a^3 
hence individual sphere radius a = 0.1/16^{1/3} unit = 0.03969 unit
(4*pi/3)*aeff^3 = N d^3
hence the interdipole spacing
d = (4*pi/3N)^{1/3}*aeff = (4*pi/3*11188)^{1/3}*aeff = 0.007207 unit

The TUC is repeated in the y and z directions with periodicity

P_y = 50d = 0.36037 unit 
P_z = 50d = 0.36037 unit

The spheres are assumed to be composed of material with refractive index
m = 1.33 + 0.01i

The wavelength is 0.6328 unit

This problem suffers from slow convergence.  Therefore the error tolerance
has been set to
TOL = 3e-4   rather than the customary 1e-5

The output file w000r000k000.sca reports a mean absorption coefficient
.027882  (fraction of the incident power that is absorbed)

The Mueller matrix elements S_11(M,N) give transmission and
reflection coefficients
T(M,N)=S_11(M,N)*|cos(theta)| for theta < 90
R(M,N)=S_11(M,N)*|cos(theta)| for theta > 90:

T( 0, 0) = 0.82913
T( 1, 0) = 0.083092*cos(48.59)
T(-1, 0) = 0.053209*cos(48.59)
T( 0, 1) = 0.031349*cos(48.59)
T( 0,-1) = 0.032172*cos(48.59)
----------------------------
T(total) = 0.96130

R( 0, 0) = 0.0019379
R( 1, 0) = 0.0011621*cos(48.59)
R(-1, 0) = 0.0029532*cos(48.59)
R( 0, 1) = 0.0096722*cos(48.59)
R( 0,-1) = 0.00099924*cos(48.59)
-------------------------------
R(total) = 0.011718

The inferred absorption coefficient is

1 - (T+R) = 0.02698

This differs by only -0.00090 from the value 0.02788 returned by
ddscat (from direct calculation of absorption), which is reasonable
given that ddscat was only required to converge to an accuracy 3.e-4
