#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------
#
# Satyariya Gupta
#
#
# -------------------------------------------------------------------
import numpy as np
# import os
# from random import*
# from itertools import chain
# import math
# import re
import sys

f1=sys.argv[1]
d = np.loadtxt(f1, dtype=float,skiprows=0)
d1 = np.genfromtxt('/home/enpeters4/optimize_new/T950_001_1e-2_used.csv', delimiter=",",dtype=float,skip_header=0)
col1G_exp=d1[:,0]
col2G_exp=d1[:,1]
epsilon_E=col1G_exp			# this is the strain
sigma_E=col2G_exp				# this is the stress
epsilon_E= np.log(1+epsilon_E)
sigma_E= sigma_E*np.exp(epsilon_E)
col1G_sim=d[:,1]
col2G_sim=d[:,2]
epsilon_S=list(np.around(np.array(col1G_sim),5))
sigma_S  =list(np.around(np.array(col2G_sim),3))
epsilon_S=list(map(abs, epsilon_S))
sigma_S  =list(map(abs, sigma_S)) 
# epsi_E=[i * 1e2 for i in epsi_E]
#epsi_Esigma_S= [i * 1e-6 for i in sigma_S]
#print (sigma_S)
#plastic_limit = [0.5e-3,1.e-3,1.5e-3,2.e-3,2.5e-3,3.e-3,3.5e-3,4.e-3,4.5e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22]
plastic_limit = [0.5e-3,1.e-3,1.5e-3,2.e-3,2.5e-3,3.e-3,3.5e-3,4.e-3,4.5e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3,0.01,0.02,0.03,0.04,0.05]
#print(plastic_limit[0])
size=len(plastic_limit)
sigma_intE =np.zeros(size)
sigma_intS =np.zeros(size)
#sigma_E_2 = []
#sigma_E_3 = []
error=0
for i in range(0,size):
	sigma_intE[i] = np.interp(plastic_limit[i], epsilon_E, sigma_E)
	sigma_intS[i] = np.interp(plastic_limit[i], epsilon_S, sigma_S)
	error = error + (sigma_intE[i]-sigma_intS[i])**2
print(np.sqrt(error/size))
# print(sigma_intE,sigma_intS)
# #print sigma_E[1]
# File1=open('Interpolate_Exp_100_stress.txt','w')
# File2=open('Interpolate_Sim_100_stress.txt','w')
# for i in range(0,18):
# 	File1.write(str(plastic_limit[i])+'\t'+str(sigma_intE[i])+'\n')
# 	File2.write(str(plastic_limit[i])+'\t'+str(sigma_intS[i])+'\n')
# File1.close()
# File2.close()



