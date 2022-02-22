from ase.io.vasp import read_vasp
import Tools as tl
from ase.build import bulk
from ase.visualize import view
import numpy as np

zr1010 = read_vasp('Zr1010_1x1.vasp')
zr1010 = tl.make_symmetric(zr1010, natoms=1)
zro2001Oterm = read_vasp('ZrO2001Oterm_2x2_0.vasp')

tl.get_adsite(zr1010, site='top', face='bottom')
tl.get_adsite(zro2001Oterm, site='top', face='top')

zr1010adstruc = tl.make_adstruc(zr1010, '1010ad', theface='bottom', thesite='top')
zro2001Otermadstruc = tl.make_adstruc(zro2001Oterm, 'rzosad', theface='top', thesite='top')

view(zr1010adstruc)
view(zro2001Otermadstruc)

# make the stak

#1: scale the slabs:

atoms1 = zro2001Oterm.copy()
atoms2 = zr1010.copy()

cell = zro2001Oterm.cell.copy()

height1 = tl.get_slab_height(zro2001Oterm)
height2 = tl.get_slab_height(zr1010)

cellheight = np.array([0, 0, height1 + height2+4])

cell[2] = cellheight
zro2001Oterm.set_cell(cell, scale_atoms=False)

auxcell = zr1010.cell.copy()
auxcell[2] = cellheight
zr1010.set_cell(auxcell, scale_atoms=False)

zr1010.pbc=False
zr1010.set_cell(np.eye(3), scale_atoms=True)
zr1010.set_cell(cell, scale_atoms=True)

scaled_adsite1 = np.multiply( np.divide(zro2001Oterm.info['adatom']['top']['top'], zro2001Oterm.cell.lengths()), cell.lengths() )
scaled_adsite2 = np.multiply( np.divide(zr1010.info['adatom']['bottom']['top'], zr1010.cell.lengths()), cell.lengths() )

atoms1.info['adatom']['top']['top'] = scaled_adsite1
atoms2.info['adatom']['bottom']['top'] = scaled_adsite2

scaled_adstruc1 = tl.make_adstruc(atoms1, 'atom1adstr')
scaled_adstruc2 = tl.make_adstruc(atoms2, 'atom2adstr', theface='bottom', thesite='top')

view(scaled_adstruc1)
view(scaled_adstruc2)
