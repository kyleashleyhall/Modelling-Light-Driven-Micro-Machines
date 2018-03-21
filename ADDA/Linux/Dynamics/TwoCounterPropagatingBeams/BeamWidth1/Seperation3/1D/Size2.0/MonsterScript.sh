#PBS -N CounterBeamDynamicsSep3Size2
# The working directory
cd /home/mv15885/Project/Size2.0
# Script goes below here
python CalculateCorrectionFactor.py
python 2BeamDynamics.py

