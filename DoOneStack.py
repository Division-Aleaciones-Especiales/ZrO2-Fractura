import pdb
import subprocess
import sys
try:
    sys.path.insert(0,'/data/git/ase/')
except FileNotFoundError as E:
    pass
try:
    import ase
except ModuleNotFoundError as E:
    subprocess.check_call([sys.executable, '-m','pip','install','git+https://gitlab.com/mdforti/ase.git'])
    import ase

from ase.io.vasp import read_vasp
from ase.build import make_supercell
from ase.build import bulk
from ase.build.surfaces_with_termination import surfaces_with_termination, atom_index_in_top, atom_index_in_bottom
from ase.geometry import get_layers
from ase.io.vasp import write_vasp
import matplotlib.pyplot as plt
plt.rc('figure', figsize=(25,10))
import numpy as np
from ase.build.surfaces_with_termination import (atom_index_in_bottom, atom_index_in_top)
from Tools import get_adsite, stack, remove_bottom_atom

zrhcpa = 3.2313 # Zr.cell[0][0]
zrhcpc = 5.1479 # Zr.cell[-1][-1]
zro2 = read_vasp('Structure/POSCAR')

orthozrhcp = bulk('Zr',crystalstructure='hcp',a=zrhcpa, c = zrhcpc,orthorhombic=True)
orthozrhcp.write('OrthoZrHCP.vasp', direct=True, wrap=True, sort=True)
zrhcp = bulk('Zr',crystalstructure='hcp',a=zrhcpa, c = zrhcpc,orthorhombic=False)
zrhcp.write('OrthoZrHCP.vasp', direct=True, wrap=True, sort=True)
multiplicities = { '1x1': [[1,0,0],[0,1,0],[0,0,1]], '2x2':[[2,0,0],[0,2,0],[0,0,1]] }

oxideterm = {'Oterm': 'O', 'Zrterm': 'Zr'}

thisoxideterm = 'Oterm'
oxidesite = 'top'

oxidename ='ZrO20001'+thisoxideterm

ZrO2001 = {thisoxideterm: {
            '1x1': surfaces_with_termination(
                zro2, [0,0,1],layers=4, vacuum=15., termination=oxideterm[thisoxideterm],verbose=True, symmetric=True
                ) } 
            }

write_1x1_oxide_surfaces = [
        write_vasp(f'ZrO2001{thisoxideterm}_1x1_{i}.vasp', thisone, sort=True, direct=True) for i, thisone in enumerate(ZrO2001[thisoxideterm]['1x1'])
        ]

ZrO2001[thisoxideterm]['2x2'] = [
        make_supercell(theoxidesurf,multiplicities['2x2'],tol=1e-10, wrap=True) 
        for i, theoxidesurf in enumerate(ZrO2001[thisoxideterm]['1x1']) 
        ]

translate = [this.translate([this[0].x]) for this in ZrO2001[thisoxideterm]['2x2']]
center = [this.center() for this in ZrO2001[thisoxideterm]['2x2']]
findtheoxidesites  = [get_adsite(theintf,site=oxidesite, face='top') for theintf in ZrO2001[thisoxideterm]['2x2']]
writetheoxides = [this.write(f'ZrO2{thisoxideterm}_2x2_{i}.vasp',sort=True, direct=True, format='vasp') for i,this in enumerate(ZrO2001[thisoxideterm]['2x2'])]


metalindex='0001'
metalmulti ='2x2'
metalname = f'Zr{metalindex}_{metalmulti}'
metalfilename = f'{metalname}.vasp'
metalsite = 'top'

Zr = {metalindex: { '1x1': surfaces_with_termination(orthozrhcp, (0,0,1), 6, vacuum=10, termination='Zr')[0] }}
Zr[metalindex]['1x1'] = remove_bottom_atom(Zr[metalindex]['1x1'])
Zr[metalindex]['1x1'] = remove_bottom_atom(Zr[metalindex]['1x1'])
Zr[metalindex][metalmulti] = make_supercell(Zr[metalindex]['1x1'], multiplicities[metalmulti],)
Zr[metalindex][metalmulti].write(metalfilename,sort=True,wrap=True, direct=True,vasp5=True)

adsites = get_adsite(Zr[metalindex][ '2x2' ],site='top',face='bottom')

def poscarname(thei): 
    return f'{oxidename}_{oxidesite}_2x2_{thei}_{metalname}_{metalsite}_2x2.vasp'

thestacks = [
        stack(oxidesurf, Zr[metalindex][metalmulti], 'top','top',3.0, cell = oxidesurf.cell.copy())  
        for oxidesurf in  ZrO2001[thisoxideterm]['2x2']
        ]

poscarnames = [poscarname(thisi) for thisi in range(len(thestacks))]

writethestack = [ thistack.write(poscarnames[i],format='vasp',wrap=True,direct=True,sort=True) for i, thistack in enumerate(thestacks)]

import subprocess
cmd_opts = ' '.join(poscarnames)
print(cmd_opts)
subprocess.run(['vesta']+poscarnames, shell=False)
