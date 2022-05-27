import os, sys
import numpy as np
import math as mh 
from odbAccess import *
from abaqusConstants import *

print('Starting')
name = '10x10x10_tension_Y_pbc_singleCrystal'       #odb name here
path = (os.getcwd())
odbName = '%s/%s.odb'%(path,name)    
odb= openOdb(odbName,readOnly=True)


myAssembly = odb.rootAssembly.instances['PART-1-1']  # name of the part to be written

newpath = 'postProc' #Output folder name 
if not os.path.exists(newpath):
    os.makedirs(newpath)

steps1=odb.steps['Step-1'].frames  # Write the step name here.


## Name of the files of the output to be written
frf1 = open("%s/%s/S11_LE11.txt"% (path,newpath), "w")
frf2 = open("%s/%s/S22_LE22.txt"% (path,newpath), "w")
frf3 = open("%s/%s/S33_LE33.txt"% (path,newpath), "w")   
frf4 = open("%s/%s/accstrain_LE22.txt"% (path,newpath),"w")
s_eqv1= open("%s/%s/eqv_stress_strain.txt"% (path,newpath), "w")

########################################################################################################################
########################################################################################################################
#    Calculating Equivalent strain from components of Log strain[LE11,LE22,LE33,LE12]
########################################################################################################################

def eqv_strain(Str_dis0,Str_dis1,Str_dis2,Str_dis3):
          straineqv=0.
          strainM=np.zeros((3,3))
          strainM[0,0]=Str_dis0
          strainM[1,1]=Str_dis1
          strainM[2,2]=Str_dis2
          strainM[0,1]=strainM[1,0]=Str_dis3
          strainM[0,2]=strainM[2,0]=0.0
          strainM[1,2]=strainM[2,1]=0.0
          eh=0.
          I2=0.
          for i in range(3):
              eh = eh +  1./3.*strainM[i,i]

          for i in range(3):
              for j in range(3):
                     #Invariant
                     if i==j:
                            I2 = I2 + (strainM[i,j]-eh)*(strainM[i,j]-eh)
                     else:
                            I2 = I2 + strainM[i,j]*strainM[i,j]
          straineqv = np.sqrt(2./3.*I2)
          return straineqv

########################################################################################################################
########################################################################################################################
########################################################################################################################
#    Calculating Equivalent Stress from stress components [s11,s22,s33,s12]
########################################################################################################################
    
def eqv_stress(Str_hom0,Str_hom1,Str_hom2,Str_hom3):
            sigma_E= 0.
            sigmaM=np.zeros((3,3))
            sigmaM[0,0]=Str_hom0
            sigmaM[1,1]=Str_hom1
            sigmaM[2,2]=Str_hom2
            sigmaM[0,1]=sigmaM[1,0]=Str_hom3
            sigmaM[0,2]=sigmaM[2,0]=0.0
            sigmaM[1,2]=sigmaM[2,1]=0.0
            Triax_xyz=0.
            Lode_xyz=0.
            p=0. 
            I2=0.   
            for i in range(3):
                p = p +  1./3.*sigmaM[i,i]
            for i in range(3):
              for j in range(3):
                 if i==j:
                    I2 = I2 + (sigmaM[i,j]-p)*(sigmaM[i,j]-p)
                 else:
                    I2 = I2 + sigmaM[i,j]*sigmaM[i,j]
      
  
            sigma_E = np.sqrt(3./2.*I2)
            return  sigma_E  

########################################################################################################################       
for cframe in range(0,len(steps1),1): 
        currentframe1= steps1[cframe]
        sumrf = [] 
        label_rf_x = []
        label_rf_y = []
        label_rf_z = []
        log_strain_components = []           # 6 components of strains e11,e22,e33,e12,e13,e23 for each element
        strain_eqv = 0.
        stress_eqv=0.   
        stress_components = []                # 6 components of stress s11,s22,s33,s12,s13,s23 for each element
        stress_mises = []                     
        volume=[]
        accgamma= []
        accgamma_sum=0.
        volume_sum=0.                              
        s11=0.
        s22=0.
        s33=0.
        s12=0.
        s13=0.
        s23=0.
        e11=0.
        e22=0.
        e33=0.
        e12=0.
        e13=0.
        e23=0.
        smises=0.
        lode=0.    
        i=0            

#################################################################################################################
#      extracting stress and strain and volume for all the elements
#################################################################################################################  
        feildvalues_stress  = currentframe1.fieldOutputs['S'] 
        feildvalues_strain  = currentframe1.fieldOutputs['LE'] 
        feildvalues_volume = currentframe1.fieldOutputs['IVOL'] 
        feildvalues_accgamma = currentframe1.fieldOutputs['SDV56']

        vol_set= feildvalues_volume.values
        S_set= feildvalues_stress.values
        LE_set= feildvalues_strain.values
        acc_ga_set= feildvalues_accgamma.values
        for p1 in LE_set:
          log_strain_components.append(p1.data) 
        for p in vol_set:
          volume.append(p.data) 
        for g in acc_ga_set:
          accgamma.append(g.data)
        for v in S_set:
          stress_components.append(v.data)
        for io in range(len(myAssembly.elements)): 
           s11=s11+stress_components[io][0]*volume[io]
           s22=s22+stress_components[io][1]*volume[io]
           s33=s33+stress_components[io][2]*volume[io]
           s12=s12+stress_components[io][3]*volume[io]
           volume_sum=volume_sum+volume[io] 
           e11=e11+log_strain_components[io][0]*volume[io]
           e22=e22+log_strain_components[io][1]*volume[io]
           e33=e33+log_strain_components[io][2]*volume[io]
           e12=e12+log_strain_components[io][3]*volume[io]
           accgamma_sum=accgamma_sum+accgamma[io]*volume[io]
        stress_eqv=eqv_stress((s11/volume_sum),(s22/volume_sum),(s33/volume_sum),(s12/volume_sum))
        strain_eqv=eqv_strain((e11/volume_sum),(e22/volume_sum),(e33/volume_sum),(e12/volume_sum))

  
  
        frf1.write( ' %.3f '%currentframe1.frameValue  + '   ' +'%.5f'%(e11/volume_sum) + '   ' +'%.5f'%(s11/volume_sum) + "\n"  )  
        frf2.write(' %.3f '%currentframe1.frameValue  + '   ' +'%.5f'%(e22/volume_sum) + '   ' +'%.5f'%(s22/volume_sum) + "\n"  )  
        frf3.write( ' %.3f '%currentframe1.frameValue  + '   ' +'%.5f'%(e33/volume_sum) + '   ' +'%.5f'%(s33/volume_sum) + "\n"  ) 
        frf4.write( ' %.3f '%currentframe1.frameValue  + '   ' +'%.5f'%(e22/volume_sum) + '   ' +'%.5f'%(accgamma_sum/volume_sum) + "\n"  )
        s_eqv1.write( ' %.3f '%currentframe1.frameValue  + '   ' + '%.5f'%(strain_eqv) +'   ' + '%.5f'%(stress_eqv) + "\n"  )



s_eqv1.close()
frf2.close()
frf1.close()
frf3.close()
frf4.close()


print('Sucessfully completed')















