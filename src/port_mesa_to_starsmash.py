import os
import matplotlib.pyplot as plt
import mesa_pars as ms
import numpy as np
import sys, getopt
import argparse
    
parser = argparse.ArgumentParser()
parser.add_argument("-f", help="-f filename",type=str)
args=parser.parse_args()

profilename = args.f

print ('Work with the file  {0}'.format(profilename))
if not os.path.exists(profilename):
    print('ERROR: this profile file is not found')


p1 = ms.mesa_profile(profilename)
x1=np.array(p1.get('mass'))
x2=np.array(p1.get('logR'))
x3=np.array(p1.get('logRho'))
x4=np.array(p1.get('logP'))
x5=np.array(p1.get('x'))
x6=np.array(p1.get('y'))
x7=np.array(p1.get('z'))
x8=np.array(p1.get('entropy'))
x9=np.array(p1.get('energy'))
x10=np.array(p1.get('logT'))

zipped = list(zip(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10))
np.savetxt('profile.starsmash', zipped, newline='\n')

