#PBS -N DeafulatADDA_Barton5_200Points
#PBS -l cput=3:00:00
# The working directory
cd /usersc/kh13490/Project/Modelling-Light-Driven-Micro-Machines/ADDA/Linux/SimpleParticleTracking/VaryingBeamSpatialIncrement/DefaultADDA_Barton5_200Points
# Script goes below here
python CalculateDDAForces.py
