#!/usr/bin/env python3
#
#   rotate_mol: Rotate a molecule in an xyz file 
#     by arbitrary angles in the moments of intertia 
#     reference orientation
# 
import sys     
import os   
import re
import numpy as np
from numpy import linalg as LA

print('''
This script rotates a molecule given as xyz file within its
moments of intertia rotation frame. 
You can either just let it rotate to the standard rotation
or give the Euler angles explicitly for a new rotation.

Usage: rotate_mol name.xyz [prec=precession nut=nutation rot=rotation]
''')


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


#
#   Calculate Euler rotation matrix for given angles (precision, nutation, rotation)
#
def EulerMatrix(pre,nut,rot):
   Emat=np.zeros((3,3))
   Emat[0][0] = np.cos(rot)*np.cos(pre)-np.sin(rot)*np.cos(nut)*np.sin(pre)
   Emat[0][1] = np.cos(rot)*np.sin(pre)+np.sin(rot)*np.cos(nut)*np.cos(pre)
   Emat[0][2] = np.sin(rot)*np.sin(nut)
   Emat[1][0] = -np.sin(rot)*np.cos(pre)-np.cos(rot)*np.cos(nut)*np.sin(pre)
   Emat[1][1] = -np.sin(rot)*np.sin(pre)+np.cos(rot)*np.cos(nut)*np.cos(pre)
   Emat[1][2] = np.cos(rot)*np.sin(nut)
   Emat[2][0] = np.sin(nut)*np.sin(pre)
   Emat[2][1] = -np.sin(nut)*np.cos(pre)
   Emat[2][2] = np.cos(nut)
   return Emat


#
#    Class molecule: for the current molecuar structure
#
class Molecule:
#
#    The general constructor
#
   def __init__(self, IL_prec, IL_nut, IL_rot):
      self.prec = IL_prec
      self.nut = IL_nut
      self.rot = IL_rot
      self.natoms=0
#
#    Read in the structure of the molecule
#
   def ReadMolecule(self,filename):
#
#    The attributes for structure and element names
#
      self.xyz = np.zeros((1,3))
      self.names=[]
      with open(filename) as infile:
         line = infile.readline().rstrip("\n")
         self.natoms=int(line)
         self.xyz = np.zeros((self.natoms,3))
         line = infile.readline().rstrip("\n")
         for i in range(self.natoms):
            line = infile.readline().rstrip("\n")
            self.names.append(line.rstrip().split()[0])
            xyz_read = line.rstrip().split()[1:4]
            for j in range(3):
               self.xyz[i][j]=float(xyz_read[j])
#
#    Print the current structure of the molecule
# 
   def PrintMolecule(self,filename):
      print(self.natoms,"\n", file=open(filename, 'w'))
      for i in range(self.natoms):
         print(self.names[i]," ",self.xyz[i][0]," ",self.xyz[i][1],self.xyz[i][2],file=open(filename,"a"))

#
#   Calculate the center of mass (COM) of a molecule
#   Optional flag: move the molecule to the COM if desired
#
   def COM_Molecule(self):
      self.com=np.zeros(3)
      self.mass=0.0
      for i in range(self.natoms):
         for j in range(3):
            self.com[j]=self.com[j]+elements_dict[self.names[i]]*self.xyz[i][j]
         self.mass=self.mass+elements_dict[self.names[i]]

      for i in range(3):
         self.com[i]=self.com[i]/self.mass
      if (True):
         for i in range(self.natoms):
            for j in range(3):
               self.xyz[i][j]=self.xyz[i][j]-self.com[j]


#
#   Calculate the moments of inertia tensor of the molecule 
#
   def MIT_Molecule(self):

      self.MIT=np.zeros((3,3))
      for i in range(self.natoms):
   # the diagonal elements
         self.MIT[0][0]=self.MIT[0][0]+elements_dict[self.names[i]]*\
                  (self.xyz[i][1]**2+self.xyz[i][2]**2)
         self.MIT[1][1]=self.MIT[1][1]+elements_dict[self.names[i]]*\
                  (self.xyz[i][0]**2+self.xyz[i][2]**2)
         self.MIT[2][2]=self.MIT[2][2]+elements_dict[self.names[i]]*\
                  (self.xyz[i][0]**2+self.xyz[i][1]**2)
   # the off-diagonal elements
         self.MIT[0][1]=self.MIT[0][1]-elements_dict[self.names[i]]*\
                  (self.xyz[i][0]*self.xyz[i][1])
         self.MIT[0][2]=self.MIT[0][2]-elements_dict[self.names[i]]*\
                  (self.xyz[i][0]*self.xyz[i][2])
         self.MIT[1][2]=self.MIT[1][2]-elements_dict[self.names[i]]*\
                  (self.xyz[i][1]*self.xyz[i][2])
# the mirrored off-diagonal elements
         self.MIT[1][0]=self.MIT[0][1]
         self.MIT[2][0]=self.MIT[0][2]
         self.MIT[2][1]=self.MIT[1][2]

#
#   Rotate a molecule with respect to the given rotation matrix
#
   def Rot_Molecule(self,matrix):
      xyz_rot = np.zeros((self.natoms,3))
      for i in range(self.natoms):
         coord = np.zeros(3)
         coord_new = np.zeros(3)
         for j in range(3):
            coord[j]=self.xyz[i][j]
            coord_new=np.dot(matrix,coord)
         for j in range(3):
            xyz_rot[i][j]=coord_new[j]
      self.xyz = xyz_rot
#  End of Class molecule


# Define class object of the new molecule first
# Filename is first command line argument
filename=sys.argv[1]

# Default values for Euler rotation angles
prec=0
nut=0
rot=0

# Read in the three angles (optional)
for arg in sys.argv:
   if re.search("=",arg):
      param,actval=arg.split("=")
      if param == "-prec":
         prec=float(actval)
      if param == "-nut":
         nut=float(actval)
      if param == "-rot":
         rot=float(actval)

#   Define the molecule object (constructor)
mol_act=Molecule(prec,nut,rot)

print("Read in file ",filename,"...")

if os.path.isfile(filename):
#   Read in the coordinates and elements of the molecule
   mol_act.ReadMolecule(filename)
else:
   print("The file ",filename," does not exist!")
   sys.exit(1)

print("Requested Euler rotation angles:")
print(" * precession: ",prec,"°")
print(" * nutation: ",nut,"°")
print(" * rotation: ",rot,"°")

print("\nPerform the rotation...")
#   Calculate center of mass of the molecule
mol_act.COM_Molecule()

#   Calculate the moments of intertia tensor of the molecule
mol_act.MIT_Molecule()

# solve the eigenvalue problem to get the principial axis of inertia
eigvals,eigvecs = LA.eig(mol_act.MIT)

# Transpose eigenvector matrix for vector transformation
eigvecs=eigvecs.transpose()

# Now rotate all position vectors of the molecule along the prinicipial axis of inertia 
mol_act.Rot_Molecule(eigvecs)
#
#   Now rotate the molecule to the desired orientation
#
Emat=np.zeros((3,3))
Emat=EulerMatrix(mol_act.prec,mol_act.nut,mol_act.rot)

mol_act.Rot_Molecule(Emat)

#
#   Print out the xyz file to struc_rotated.xyz
#

mol_act.PrintMolecule("struc_rotated.xyz")
print("\nRotated molecule written to file 'struc_rotated.xyz'\n")
