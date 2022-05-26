#!/usr/bin/env python3

import numpy as np
#import sys
import os,shutil
import time
import re
import Optimization
from optparse import OptionParser
from subprocess import run,check_output
np.set_printoptions(linewidth=np.inf)
#-------------------------------------------------------------------------------------------------#

_ratio = None

#calls Optimization.py 
class optimize(Optimization.Optimization):

  def id(self,x):
    return str(self.map2space(x[:self.dimension])).translate(str.maketrans(' []','___')) 
#===========================================#

#writes output files
  def fitness(self,x):
    time.sleep(np.random.randint(2,10)*0)
    if self.id(x) in self.locations:
      self.curr_locations.append(np.append(x,self.locations[self.id(x)]))
      return self.locations[self.id(x)]
    with open('{}/output_gen{}_{}.log'.format(options.root,self.generation+1,self.id(x)),'a') as file:
      file.write("\n Generation %i "%(self.generation+1))

#2d array. Model parameters but not the ones being optimized. These will change depending on what experiment is being simulated. 
#N simulations for one set of conditions at a time. 
#Folders will be created and named using parameter space coordinates. 
#replace temp_strts in .inp file

#2 test conditions in array (two lists within matrix)
      test_conds = np.array([ [1123.0,0.83769E-08,0.6514E+09,90000.0,1.9e-4,0.75,5.73E-4,5.8E-5,8.2e-5],
                              [1223.0,0.83269E-08,0.712E+09,86000.0,1.9e-4,0.71,5.62E-4,6.7E-5,7.8e-5] 
                           ])                            
      
# List that stores error for each set of test conditions.
      error = [] 

#enumerate function assigns integer with list entities 
#for loop runs once for one set of test conditions 
#General logistics. Setting up files and folders. 
#make dirtectory for each optimation parameter

      for test_cond_run,test_cond in enumerate(test_conds):
        dir_loc = "{}/Generation_{}/Test_{}/test_cond_{}".format(self.root,str(self.generation+1),self.id(x),str(test_cond_run+1))
        os.makedirs(dir_loc,exist_ok=True)
        msg = '\n'.join(['generation {} location {} test_cond {}'.format(self.generation+1,self.id(x),str(test_cond_run+1))])
        print(msg)
        file.write(msg+'\n')

        rep = {} # dictionary representing the parameters
        for i in range(len(x[:self.dimension])): 
          rep["%i_coords"%(i+1)] = self.map2space(x[:self.dimension])[i]                             # parameters are in MPa  
#This loop runs as many times as there are parameters in test conditions. 9 is hard-coded. 
#        for j in range(9):
#Change to len(test_conds[0]).
        for j in range(len(test_conds[0])):
          rep["{}_temp_strts".format(j+1)] = test_cond[j]                                            # parameters are in MPa 
        print("rep contains",rep)
        print("rep.items() contains",rep.items())
        
        if self.job_id != None:
#Copies being made into simulation. Once for every test condition.
#          shutil.copy("%s/%s.inp"%(self.root,self.job_id),'%s/'%dir_loc)       #  *.inp file for Abaqus
          shutil.copy("%s/%s.f"%(self.root,self.umat),'%s/'%dir_loc)
          shutil.copy("postproc.py",'%s/'%dir_loc)
          shutil.copy("Read_euler.f",'%s/'%dir_loc)
          shutil.copy("trial.txt",'%s/'%dir_loc)
          print("%s/%s.inp"%(self.root,self.job_id))
          print("%s/"%dir_loc)

#         with open('{}/{}'.format(self.root,self.mat_file),'r') as file_in:
#           contents = file_in.read()
#         for key,value in list(rep.items()):
#           print("replacing with value",str(value))
#           print("string key",str(key))
#           contents = contents.replace(str(key),str(value) )
#         with open("{}/Section_make.f".format(dir_loc),'w') as f_out:
#           f_out.write(contents) 

        with open('{}/{}'.format(self.root,self.disp_file),'r') as file_in:
          # Certain strings in the input file of the simulation are being
          # replaced with actual parameter values. Examples of such strings : 1_temp_strts etc.
          contents1 = file_in.read()
        for key,value in list(rep.items()):
          print("replacing with value",str(value))
          print("string key",str(key))
          contents1=re.sub(r'\b'+str(key)+r'\b', str(value), contents1)
          #contents1 = contents1.replace(str(key),str(value),1)
        #print(contents1)
        with open("{}/10x10x10_tension_Y_pbc_singleCrystal.inp".format(dir_loc),'w') as f_out:
          # New input file written out with strings replaced by numerical
          # values of parameters.
          f_out.write(contents1) 

#--------- Running damask postProcessing and finding error against "experimental" data ----------- #
#         print('Start checking server load')
#         availability={}
#         cpu_num = 4
#         for i in [1,2,6]:
#           cmd0 = 'ssh compute{:02}.egr.msu.edu -t "nc --recv-only localhost 51003"'.format(i)
#           cores = check_output(cmd0, shell=True, universal_newlines=True)
#           cmd0_1 = 'ssh compute{:02}.egr.msu.edu -t "nc --recv-only localhost 51001"'.format(i)
#           load = check_output(cmd0_1, shell=True, universal_newlines=True)
#           availability[i] = float(cores.split()[0]) - float(load.split()[0])  
#         j = max(availability, key=availability.get)
#         print('compute{:02} is least loaded'.format(j)) 
         
#         cmd1 = 'ssh compute{:02}.egr.msu.edu -t "cd {};export DAMASK_NUM_THREADS={ncpus};./postplot.py"'.format(j,dir_loc,ncpus=cpu_num,jobs=self.job_id)
#         cmd1 = 'cd {};gfortran Section_make.f; ./a.out'.format(dir_loc)
#         print(cmd1)
#         run([cmd1],shell=True)

#NO INPUT??---------------
#string being constructed that launches the job
        cmd2 = 'cd {};abaqus job={} user={}.f'.format(dir_loc,self.job_id,self.umat)
        print(cmd2)
        run([cmd2],shell=True)
        tt=0.0
        while not os.path.exists('{}/{}.sta'.format(dir_loc,self.job_id)):
            time.sleep(1)

        if os.path.exists('{}/{}.sta'.format(dir_loc,self.job_id)):
            # Until .sta file is generated, calculation has not started.
            while not os.path.exists('{}/done.txt'.format(dir_loc)):
#we are waiting for along enough (like 1.5 hrs) for high temperature and high strain rate simulations 
                time.sleep(2)
                tt= tt+2
                if tt > 1800:         
                    break
                
#if 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY' in open('{}/10x10x10_tension_pbc_d_5um.sta'.format(dir_loc)).read():
            if os.path.exists('{}/done.txt'.format(dir_loc)):
                cmd3 = 'cd {};abaqus python postproc.py'.format(dir_loc)
                print(cmd3)
                run([cmd3],shell=True)
                cmd4 = 'cd {}; python3 {}/Interpolate_{}.py postProc/S22_LE22.txt'.format(dir_loc,self.root,test_cond_run)
                print(cmd4)
                run([cmd4],shell=True)
                fitness_value = check_output(cmd4, shell=True, universal_newlines=True)
                print("error fitness_value  ",fitness_value)
                error.append(float(fitness_value))
#        elif 'Abaqus/Analysis exited with errors' in open('{}/10x10x10_tension_pbc_d_5um.log'.format(dir_loc)).read():
            else:
                error.append(float(1.e9))
#------------------------------------------------------------------------------------------------- #
      error = np.array(error)
      avg_error = np.mean(error)
#       if np.count_nonzero(np.isnan(error)) <= 0.5*len(error):
#         avg_error = np.nanmean(error)
#       else: 
#         avg_error = np.mean(error)

      print("+++++++++++++++++++++++++++++++++ current fitness and points +++++++++++++++++++++++++++")
      print(float(avg_error))
      print(x[:self.dimension])
      print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
      file.write("\n +++++++++++++++++++++++++++++++++ current fitness and points +++++++++++++++++++++++++++\n")
      file.write("\n {}".format(str(avg_error)))
      file.write("\n {}".format(x[:self.dimension]))
      file.write("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
      self.curr_locations.append(np.append(x,float(avg_error)))
      self.locations[self.id(x)] = float(avg_error)
    
    return float(avg_error)
#    return np.linalg.norm(np.dot(x,x))


  def info_fitness(self,disp_file,job_id,umat):
    self.disp_file  = disp_file
    self.job_id     = job_id
    self.umat       = umat
#    self.exp_file  = exp_file

#-------------------------------- main program starts here ----------------------------------------- #
#adds more options to command

parser = OptionParser()
parser.add_option(      '--mat',
                  dest = 'mat_file',
                  type = 'string', metavar = 'string',
                  help = ' material.config file (including suffix)')
parser.add_option(      '--inp',
                  dest = 'abaqusinput_file',
                  type = 'string', metavar = 'string',
                  help = ' abaqus_input file (including suffix)')
parser.add_option(      '--job_id',
                  dest = 'job_id',
                  type = 'string', metavar = 'string',
                  help = ' job id for abaqus input file')
parser.add_option(      '--root',
                  dest = 'root',
                  type = 'string', metavar = 'string',
                  help = ' desired root of this process ')
parser.add_option('--restart',  action="store_true",
                  dest="restart",
                  help="restart optimization")
parser.add_option(      '--points',
                  dest = 'points_data',
                  type = 'string', metavar = 'string',
                  help = 'points for next generation ')
parser.add_option(      '--pBest',
                  dest = 'posPBest_data',
                  type = 'string', metavar = 'string',
                  help = 'posPBest for current generation ')
parser.add_option(      '--gBest',
                  dest = 'posGBest_data',
                  type = 'string', metavar = 'string',
                  help = 'posGBest for current generation ')
parser.add_option(      '--velo',
                  dest = 'velocities_data',
                  type = 'string', metavar = 'string',
                  help = 'velocities for current generation ')

#--------------Change these default values?-----------------
#making the default values and let them show
parser.set_defaults(disp_file = 'abaqusfile.inp',
		#job_id = 'Testjob',
                    job_id = '10x10x10_tension_Y_pbc_singleCrystal',
                    umat = 'umat_gamma_gammaP_May2020'
                   )

(options,filenames) = parser.parse_args()

if not options.mat_file or not os.path.exists(options.mat_file):
  print("Suitable format material config (file containing parameters) is not supplied ")

#if not os.path.exists('%s.load'%options.job_id):
#  parser.error('No job_file selected for optimization')

options.root = os.path.dirname(os.path.realpath(__file__)) if options.root == None else options.root

tick = time.time()

#make choices for bounds for paramters in optimize matrix
#bounds array is made up of parameters being optimized.
    
#optimization occurs here
#pso = particle swarm optimization   
theOptimizer = optimize( method = 'pso',
                         bounds = np.array([[1.e+7,5.e+7],
                                            [1.e-18,1.e-17],                    
                                            [8.e-19,1.5e-18],
                                            [1.e-6,5.e-6]
                                           ]),
                         tolerance = 0.050,
                         root      = options.root,
                         rigid     = True,
                         )

theOptimizer.info_fitness(options.disp_file, options.job_id, options.umat)
theOptimizer.optimize(verbose = False)
tock = time.time()
print("Time for simulation",(tock - tick))
print(theOptimizer.cost())
print(theOptimizer.best())
with open("{}/output.log".format(options.root),'a') as file:
  file.write("\nTime for simulation {}".format(tock - tick))
  file.write("\n {}".format(theOptimizer.cost()))
  file.write("\n {}".format(theOptimizer.best()))
