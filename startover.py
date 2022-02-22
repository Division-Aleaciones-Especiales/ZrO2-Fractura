from ase.io.vasp import read_vasp
import Tools as tl
from ase.build import bulk
from ase.visualize import view

zr1010 = read_vasp('Zr1010_2x2.vasp')
zro20001Oterm = read_vasp('ZrO2001Oterm_2x2_0.vasp')


view(zr1010)
view(zro20001Oterm)
