#Preliminary Imports
import os
import subprocess

#Preliminary Variables
x=0
z=5
Initial_y=-5
Step_y=1
Final_y=5

y=Initial_y #Set the value of y to the initial value

while (y<Final_y):
  callString="../adda -size 2 -lambda 1 -prop 0 0 1 -beam barton5 1 "+str(x)+" "+str(y)+" "+str(z)+" -store_beam -store_dip_pol -store_int_field"
  print callString
  subprocess.call(callString,shell=True)
  y=y+Step_y
