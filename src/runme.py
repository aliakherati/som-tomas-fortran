#! /usr/bin/python
# This python script run SOM-TOMAS with various inputs

import os
import numpy as np
import time 

startTime = time.time()

precursor = 'nonorg_alpha-pinene'
regime = 'lownox'
VWL = 0 # [0 or 1] the switch for On/Off vapor wall loss
PWL = 0 # [0 or 1] the switch for On/Off particle wall loss

ke = 0.075 # loss rate constant
kw0= 0.0 # loss rate constant

OH_conc = 1.5E06 # [Molecules/cm3]
nucrate = 0.0 # [cm-3 s-1] nucleation rate
dilt = 0.0 # dilution rate
boxvol = 7000000.0 # teflon [cm3]

endtime = 4.0 # [hours] hours of run time (make even hours)
alpha = 0.001 # accommodation coefficient
Dbk = 1.0e-10 # particle-phase diffusion coefficient [m2/s]
kc = 0.0 # first-order loss rate of species in the particle phase [1/s]
storg = 0.025 # [N/m] surface tension

pres = 101325.0 # [Pa] pressure
temp = 298.0 # [K] temperature
rh = 0.2 # relative humidity
No = 50000.0 # [# cm-3] background number conc.
Dpm = 0.100 # [microns] background median diameter
sigma = 1.1 # background sigma

precemiss = 1 # [0 or 1] if there is any precursor
ninpprec = 1 # [integer] number of input precursors
precname = 'GENVOC' # precursor names
preclenname = 6 # length of the precursor's name (how many letters+spaces)
ippmprec = 0.05 # [ppm] emission of the precursors

oxyspemiss = 0 # [0 or 1] if there is any emission for oxygenated species
ninpoxysp = 1 # [integer] number of inputs for oxygenated species
inpc = 12 # [integer] number of carbons for oxygenated species
inpo = 1 # [integr] number of oxygen for oxygenated species
ippmoxy = 0.05 # [ppm] emission of the oxygenated species

nsomprec = 1 # number of som precursors for parameterizations
somprecname = 'GENSOMG' # AR1SOMG AR2SOMG' # som precursor parameterization's name
dlvp = 1.91 # dLVP values for species above

if os.path.exists('input'):
   os.system('rm input')
   
f = open('input','w')
rname = 'test_%s_%s_%04.1f_%08.1f_%05.3f_%7.1E_%5.3f'%(precursor,regime,endtime,No,Dpm,OH_conc,ippmprec)
#rname = precursor+'_'+regime+'_'+tend+'_'+no+'_'+dp+'_'+oh+'_'+ppm
f.write('%s\n'%rname)
f.write('%1i\n'%VWL)
f.write('%1i\n'%PWL)
f.write('%f\n'%ke)
f.write('%f\n'%kw0)
f.write('%f\n'%OH_conc) # ***
f.write('%7.1f\n'%nucrate)
f.write('%7.1f\n'%dilt)
f.write('%15.1f\n'%boxvol)
f.write('%8.2f\n'%endtime) # ***
f.write('%8.2f\n'%alpha)
f.write('%17.15f\n'%Dbk)
f.write('%10.8f\n'%kc)
f.write('%8.5f\n'%storg)
f.write('%12.3f\n'%pres)
f.write('%8.5f\n'%temp)
f.write('%8.5f\n'%rh)
f.write('%8.5f\n'%No) # ***
f.write('%8.5f\n'%Dpm) # ***
f.write('%8.5g\n'%sigma)

f.write('%1i\n'%precemiss)
if precemiss == 1:
   f.write('%3i\n'%ninpprec)
   f.write('%s\n'%precname)
   f.write('%s\n'%preclenname)
   f.write('%s\n'%ippmprec) # ***
   
f.write('%1i\n'%oxyspemiss)
if oxyspemiss == 1:
   f.write('%s\n'%ninpoxysp)
   f.write('%s\n'%inpc)
   f.write('%s\n'%inpo)
   f.write('%s\n'%ippmoxy)
            
f.write('%02i\n'%nsomprec)
f.write('%s\n'%somprecname)
f.write('%6.3f\n'%dlvp)
f.close()

run_time1 = time.time()
#if os.path.exists('%s.out'rname):
#   os.system('rm %s*'%(rname))
os.system('./box.exe < input > ../outputs/'+rname+'.out')
#os.system('rm input')
run_time2 = time.time()

mins = (run_time2 - run_time1)/60.
secs = (mins - int(mins))*60.
print('Running took %2.0f mins and %2.0f seconds.'%(mins, secs))
