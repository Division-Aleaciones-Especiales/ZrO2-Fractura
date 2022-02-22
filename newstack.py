from ase.io.vasp import read_vasp
import Tools as tl
from ase.build import bulk
from ase.visualize import view
import numpy as np
import pdb

zr1010 = read_vasp('Zr1010_1x1.vasp')
zr1010 = tl.make_symmetric(zr1010, natoms=1)
zro2001Oterm = read_vasp('ZrO2001Oterm_2x2_0.vasp')
tl.get_adsite(zr1010, site='top', face='bottom')
tl.get_adsite(zro2001Oterm, site='top', face='top')

zr1010adstruc = tl.make_adstruc(zr1010, '1010ad', theface='bottom', thesite='top')
zro2001Otermadstruc = tl.make_adstruc(zro2001Oterm, 'rzosad', theface='top', thesite='top')

# make the stak

#1: scale the slabs:
from ase import Atoms
atoms1 = zro2001Oterm.copy()
atoms2 = zr1010.copy() # Atoms([atom for atom in zr1010])

cell = zro2001Oterm.cell.copy()

height1 = tl.get_slab_height(atoms1)
height2 = tl.get_slab_height(atoms2)

cellheight = np.array([0, 0, height1 + height2+2+10])

cell[2] = cellheight

auxcell = atoms1.cell.copy()
auxcell[2] = cellheight
atoms1.set_cell(auxcell, scale_atoms=False)
atoms1.center(axis=2)

auxcell2 = atoms2.cell.copy()
auxcell2[2] = cellheight
atoms2.set_cell(auxcell2, scale_atoms = False)
atoms2 = Atoms( atoms2.get_chemical_symbols(), scaled_positions=atoms2.get_scaled_positions(), cell= cell)
atoms2.rotate(0,'z')
atoms2.center(axis=2)

tl.get_adsite(atoms1, site='top', face='top')
tl.get_adsite(atoms2, site='top', face='bottom')

scaledadstruc1 = tl.make_adstruc(atoms1, 'atom1adstr')
scaledadstruc2 = tl.make_adstruc(atoms2, 'atom2adstr', theface='bottom', thesite='top')

scaledadsite1 = atoms1.info['adatom']['top']['top']
scaledadsite2 = atoms2.info['adatom']['bottom']['top']

translation = scaledadsite1-scaledadsite2 + np.array([0,0,2])

stacked = atoms2.copy() 
print(stacked.cell)
stacked.translate(translation)

stack = atoms1.copy()
stack.extend(stacked)
stack.center(axis=2)
stack.wrap(pbc = (1,1,0))
view(stack)


