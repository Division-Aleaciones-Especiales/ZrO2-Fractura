import pdb
import os, sys, subprocess
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
from ase.visualize import view
zrhcpa = 3.2313 # Zr.cell[0][0]
zrhcpc = 5.1479 # Zr.cell[-1][-1]
zro2 = read_vasp('Structure/POSCAR')

vaspopts={'format':'vasp','sort':True,'wrap':True, 'direct':True,'vasp5':True}

orthozrhcp = bulk('Zr',crystalstructure='hcp',a=zrhcpa, c = zrhcpc,orthorhombic=True)
orthozrhcp.write('OrthoZrHCP.vasp', direct=True, wrap=True, sort=True)
zrhcp = bulk('Zr',crystalstructure='hcp',a=zrhcpa, c = zrhcpc,orthorhombic=False)
zrhcp.write('OrthoZrHCP.vasp', **vaspopts)
multiplicities = { '1x1': [[1,0,0],[0,1,0],[0,0,1]], '2x2':[[2,0,0],[0,2,0],[0,0,1]] }

metalindex={'0001': ( 0,0,1 ), '11m20':( 1,0,0 ), '1010': (0,1,0)}
metalmulti ='1x1'
metalsite = 'top'
thisindex='11m20'
metalname = f'Zr{thisindex}_{metalmulti}'
metalfilename = f'{metalname}.vasp'

def poscarname(adistance, anoxidesite, anoxidename): 
    return f'{anoxidename}_{anoxidesite}_2x2_{metalname}_{metalsite}_2x2_{adistance}.vasp'

Zr = {thisindex: { '1x1': surfaces_with_termination(orthozrhcp, metalindex[thisindex], 6, vacuum=10, termination='Zr')[1] }}
Zr[thisindex]['1x1'] = remove_bottom_atom(Zr[thisindex]['1x1'])
Zr[thisindex]['1x1'] = remove_bottom_atom(Zr[thisindex]['1x1'])
Zr[thisindex][metalmulti] = make_supercell(Zr[thisindex]['1x1'], multiplicities[metalmulti],)
Zr[thisindex][metalmulti].write(metalfilename, **vaspopts)
adsites = get_adsite(Zr[thisindex][metalmulti],site='top',face='bottom')

oxideterm = {'Oterm': 'O', 'Zrterm': 'Zr'}
oxidesites = ['top']#, 'hollow', 'bridge']
thisoxideterm = 'Oterm'

thisoxidename ='ZrO2001'+thisoxideterm

distances = np.array([2.5, 2.75, 3.0, 3.25, 3.5, 4.0, 5.0, 6.0])

for thisoxidesite in oxidesites:

    ZrO2001 = {thisoxideterm: {
                '1x1': surfaces_with_termination(
                    zro2, [0,0,1],layers=4, vacuum=15., termination=oxideterm[thisoxideterm],verbose=True, symmetric=True
                    ) } 
                }

    for i in range(3):
        ZrO2001[thisoxideterm]['1x1'] = [ remove_bottom_atom( thisoxideslab ) for thisoxideslab in ZrO2001[thisoxideterm]['1x1'] ]

    write_1x1_oxide_surfaces = [
            thisone.write(f'ZrO2001{thisoxideterm}_1x1_{i}.vasp', **vaspopts) for i, thisone in enumerate(ZrO2001[thisoxideterm]['1x1'])
            ]

    ZrO2001[thisoxideterm]['2x2'] = [
            make_supercell(theoxidesurf,multiplicities['2x2'],tol=1e-10, wrap=True) 
            #for i, theoxidesurf in enumerate(ZrO2001[thisoxideterm]['1x1']) 
            for theoxidesurf in ZrO2001[thisoxideterm]['1x1']
            ]

    write_2x2_oxide_surfaces = [
            thisone.write(f'ZrO2001{thisoxideterm}_2x2_{i}.vasp', **vaspopts) for i, thisone in enumerate(ZrO2001[thisoxideterm]['2x2'])
            ]

    center = [this.center(axis=2) for this in ZrO2001[thisoxideterm]['2x2']]
    findtheoxidesites  = [get_adsite(theintf,site=thisoxidesite, face='top') for theintf in ZrO2001[thisoxideterm]['2x2']]
    writetheoxides = [this.write(f'ZrO2{thisoxideterm}_2x2_{i}.vasp', **vaspopts) for i,this in enumerate(ZrO2001[thisoxideterm]['2x2'])]

    stackcell=ZrO2001[thisoxideterm]['2x2'][0].cell.copy()
    thestacks = []
    thestacks = [
            stack(ZrO2001[thisoxideterm]['2x2'][0].copy(), Zr[thisindex][metalmulti].copy(), thisoxidesite ,'top',d, cell = stackcell)   
            for d in distances
            ]

    poscarnames = [poscarname(thisdistance, thisoxidesite, thisoxidename) for thisdistance in distances]

    writethestack = [ thistack.write(poscarnames[i],format='vasp',wrap=True,direct=True,sort=True) for i, thistack in enumerate(thestacks)]

    cmd_opts = ' '.join(poscarnames)

    #os.popen('vesta '+cmd_opts)
