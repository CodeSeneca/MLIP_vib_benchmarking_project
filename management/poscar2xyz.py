#!/usr/bin/env python3.6
import sys
import os
import numpy as np

#print '''
#This script converts a CONTCAR file to a xyz 
#trajectory file. All will be done automatically.
#'''

#
#   Calculate the center of mass (COM) of a molecule
#   Optional flag: move the molecule to the COM if desired
#
# Dictionary of all elements matched with their atomic masses.
elements_dict = {'H' : 1.008,'He' : 4.003, 'Li' : 6.941, 'Be' : 9.012,\
                 'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,\
                 'F' : 18.998, 'Ne' : 20.180, 'Na' : 22.990, 'Mg' : 24.305,\
                 'Al' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,\
                 'Cl' : 35.453, 'Ar' : 39.948, 'K' : 39.098, 'Ca' : 40.078,\
                 'Sc' : 44.956, 'Ti' : 47.867, 'V' : 50.942, 'Cr' : 51.996,\
                 'Mn' : 54.938, 'Fe' : 55.845, 'Co' : 58.933, 'Ni' : 58.693,\
                 'Cu' : 63.546, 'Zn' : 65.38, 'Ga' : 69.723, 'Ge' : 72.631,\
                 'As' : 74.922, 'Se' : 78.971, 'Br' : 79.904, 'Kr' : 84.798,\
                 'Rb' : 84.468, 'Sr' : 87.62, 'Y' : 88.906, 'Zr' : 91.224,\
                 'Nb' : 92.906, 'Mo' : 95.95, 'Tc' : 98.907, 'Ru' : 101.07,\
                 'Rh' : 102.906, 'Pd' : 106.42, 'Ag' : 107.868, 'Cd' : 112.414,\
                 'In' : 114.818, 'Sn' : 118.711, 'Sb' : 121.760, 'Te' : 126.7,\
                 'I' : 126.904, 'Xe' : 131.294, 'Cs' : 132.905, 'Ba' : 137.328,\
                 'La' : 138.905, 'Ce' : 140.116, 'Pr' : 140.908, 'Nd' : 144.243,\
                 'Pm' : 144.913, 'Sm' : 150.36, 'Eu' : 151.964, 'GD' : 157.25,\
                 'Tb' : 158.925, 'Dy': 162.500, 'Ho' : 164.930, 'Er' : 167.259,\
                 'Tm' : 168.934, 'Yb' : 173.055, 'Lu' : 174.967, 'Hf' : 178.49,\
                 'Ta' : 180.948, 'W' : 183.84, 'Re' : 186.207, 'Os' : 190.23,\
                 'Ir' : 192.217, 'Pt' : 195.085, 'Au' : 196.967, 'Hg' : 200.592,\
                 'Tl' : 204.383, 'Pb' : 207.2, 'Bi' : 208.980, 'Po' : 208.982,\
                 'At' : 209.987, 'Rn' : 222.081, 'Fr' : 223.020, 'Ra' : 226.025,\
                 'Ac' : 227.028, 'Th' : 232.038, 'Pa' : 231.036, 'U' : 238.029,\
                 'Np' : 237, 'Pu' : 244, 'Am' : 243, 'Cm' : 247, 'Bk' : 247,\
                 'Ct' : 251, 'Es' : 252, 'Fm' : 257, 'Md' : 258, 'No' : 259,\
                 'Lr' : 262, 'Rf' : 261, 'Db' : 262, 'Sg' : 266, 'Bh' : 264,\
                 'Hs' : 269, 'Mt' : 268, 'Ds' : 271, 'Rg' : 272, 'Cn' : 285,\
                 'Nh' : 284, 'Fl' : 289, 'Mc' : 288, 'Lv' : 292, 'Ts' : 294,\
                 'Og' : 294}



def COM(xyz,natoms,names,elements_dict,move):
   com=np.zeros(3)
   mass=0.0
   for i in range(natoms):
      for j in range(3):
         com[j]=com[j]+elements_dict[names[i]]*xyz[i][j]
      mass=mass+elements_dict[names[i]]

   for i in range(3):
      com[i]=com[i]/mass
   if (move):
      for i in range(natoms):
         for j in range(3):
            xyz[i][j]=xyz[i][j]-com[j]
   return com




#  Determine number of lines in file 
num_lines = sum(1 for line in open('POSCAR'))

xdat_in = open("POSCAR","r")

with xdat_in as infile:
   line = infile.readline().rstrip("\n")
   line = infile.readline().rstrip("\n")
   surf_scale = float(line)
   line = infile.readline().rstrip("\n")
   a_read = line.rstrip().split()[0:3]
   line = infile.readline().rstrip("\n")
   b_read = line.rstrip().split()[0:3]
   line = infile.readline().rstrip("\n")
   c_read = line.rstrip().split()[0:3]
 
   a_vec=np.zeros(3)
   b_vec=np.zeros(3)
   c_vec=np.zeros(3)

   for i in range(3):
      a_vec[i]=float(a_read[i])
      b_vec[i]=float(b_read[i])
      c_vec[i]=float(c_read[i])
   # read in the element ordering
   line = infile.readline().rstrip("\n")
   elements = line.rstrip().split()
   nelem = len(elements)
   # read in the number of elements 
   line = infile.readline().rstrip("\n")
   elem_num0 = line.rstrip().split()
   
   elem_num=[] 
   for i in range(nelem):
      elem_num.append(int(elem_num0[i]))
   # Determine the atomic symbols for the xyz trajectory file 
   names=[] 
   for i in range(nelem):
      for j in range(elem_num[i]):
         names.append(elements[i])

   natoms = np.sum(elem_num)   
   # Determine the number of frames 
   
   xyz = np.zeros((natoms,3))
   line = infile.readline().rstrip("\n")
   line_split = line.rstrip().split()
   if line_split[0] == "Selective" or line_split[0] == "selective":
      line = infile.readline().rstrip("\n")
      line_split = line.rstrip().split()     
   #  Determine if cartesian or direct coordinates are used 
   coord_direct = False
   if line_split[0] == "Cartesian":
      coord_direct = False
   elif line_split[0] == "Direct":
      coord_direct = True
  
   for j in range(natoms):
      line = infile.readline().rstrip("\n")
      pos = line.rstrip().split()
      for k in range(3):
         xyz[j][k]=float(pos[k])

# Center the geometry at the origin (COM = 0.0 , 0.0 , 0.0)


#com_mol1 = COM(xyz,natoms,names,elements_dict,move=True)

# Write out the trajectory to file trajectory.xyz

#surf_scale = 1.0
xyz_act=np.zeros(3)
original_stdout=sys.stdout

with open("output.xyz","w") as f:
   sys.stdout = f
   for i in range(1):
      print(natoms)
      print("comment " + str(i))
      for j in range(natoms):
         if coord_direct == True: 
            xyz_act[0]=(xyz[j][0]*a_vec[0]+xyz[j][1]*b_vec[0]+\
                      xyz[j][2]*c_vec[0])*surf_scale
            xyz_act[1]=(xyz[j][0]*a_vec[1]+xyz[j][1]*b_vec[1]+\
                      xyz[j][2]*c_vec[1])*surf_scale
            xyz_act[2]=(xyz[j][0]*a_vec[2]+xyz[j][1]*b_vec[2]+\
                      xyz[j][2]*c_vec[2])*surf_scale
            print(names[j] + "     " +  str(xyz_act[0]) + "    " + str(xyz_act[1]) +\
                 "   " + str(xyz_act[2])) 
         elif coord_direct == False:
            print(names[j] + "     " +  str(xyz[j][0]*surf_scale) + "    " +\
                    str(xyz[j][1]*surf_scale) + "   " + str(xyz[j][2]*surf_scale))
 
   sys.stdout=original_stdout

print("POSCAR file conversion successful. Output stored to output.xyz.")          
