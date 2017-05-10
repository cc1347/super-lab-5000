GlowScript 2.4 VPython
import numpy as np
from visual import *

atom, co_x, co_y, co_z = np.loadtxt("data.txt", unpack= true) e

N_atom = len(atom)

for i in range (1, N_atom):

  if atom(i) in ["O","o"]:
      colour = "green"
  else if atom(i)  in ["Mg","mg"]:
      colour = "red"
  else
      colour = "blue"
      
   po_x = co_x(i)*4.2
   po_y = co_y(i)*(len(atoms)/8.0)
   po_z = co_z(i)*4.2
  

  sphere(pos = vector(po_x,po_y,po_z), radius = 0.3, color = color.colour)
