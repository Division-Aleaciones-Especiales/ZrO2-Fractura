import pdb
import os
from ase import Atoms
from ase.geometry import get_layers
import numpy as np
import sys
sys.path.insert(0, '/data/git/ase/')
from ase.build.surfaces_with_termination import atom_index_in_top, atom_index_in_bottom
from ase.build import sort

vasp_write_options = {'format':'vasp', 'direct':True, 'wrap': True}
def check_has_adatom(theatoms):
    if hasattr(theatoms, 'info'):
        if isinstance(theatoms.info, dict):
            if 'adatom' in theatoms.info.keys():
    #            print(theatoms.info['adatom'])
                pass
            else:
                theatoms.info.update({'adatom':{}})
    else:
        raise ValueError('no info property')


def get_top_positions(atoms, top_layer):
    return atoms.positions[top_layer[0]].copy()

def get_bridge_position(atoms, top_layer):
    assert(len(top_layer)>=2)
    top_layer = np.array(top_layer).astype(int)
    first_x = atoms.positions[top_layer,0].argsort()
    return  atoms.positions[top_layer[first_x[:2]], :].mean(axis=0)

def get_adsite(atoms, site = None, face='top'):
    """ 
    site = ['top', 'hollow', 'bridge']
    face= ['top', 'bottom']
    """
    check_has_adatom(atoms)
    layer, hs = get_layers(atoms, [0, 0, 1])
    if face == 'top':
        layer = atom_index_in_top(layer)
    elif face == 'bottom':
        layer = atom_index_in_bottom(layer)

    info = atoms.info

    if face not in info['adatom']:
        info['adatom'].update({face:{}})

    if site == 'hollow':
        info['adatom'][face].update({'hollow': atoms.positions[layer].mean(axis=0)})
    elif site == 'top':
        info['adatom'][face].update({'top': get_top_positions(atoms, layer)})
    elif site == 'bridge':
        info['adatom'][face].update({'bridge':get_bridge_position(atoms, layer)})
    elif site == None:
        info['adatom'][face].update( {'hollow': atoms.positions[layer].mean(axis=0)})
        info['adatom'][face].update({'top': get_top_positions(atoms, layer)})
        info['adatom'][face].update({'bridge':get_bridge_position(atoms, layer)})

    return info['adatom']

def make_symmetric(theatoms, natoms=3):
    symmetricatoms = theatoms.copy()
    for i in range(natoms):
        symmetricatoms = remove_bottom_atom(symmetricatoms) 
    return symmetricatoms

def make_adstruc(theatoms, name, theface='top', thesite='top', d=2):
    from ase.io.vasp import write_vasp
    if theface == 'bottom':
        d = -d 

    ad_pos = [theatoms.info['adatom'][theface][thesite]+[0,0,d]]
    adatom_inface_insite = Atoms('H', positions=ad_pos, pbc=True, cell=theatoms.cell.copy())
    adstruc = theatoms.copy()
    adstruc.extend(adatom_inface_insite)
    write_vasp(f'{name}_{theface}face_ad{thesite}.vasp', adstruc,direct=True, sort=True, vasp5=True, wrap=True)
    return adstruc

def get_slab_height(theatoms):
    return theatoms.positions[:,2].max() - theatoms.positions[:,2].min()

def get_scaled_site(thesite, thecell):
    return np.array( [ad / np.linalg.norm(v) for ad, v in zip(thesite, thecell.array)] )

def remove_bottom_atom(theatoms):
    atoms = theatoms.copy()
    layers, hs = get_layers(atoms,(0,0,1))
    atoms_in_bottom = atom_index_in_bottom(layers)
    atoms.pop(atoms_in_bottom[-1])
    return atoms

def stack(in_atoms1, in_atoms2, tagadsite1, tagadsite2, distance, mix=0.5, cell = None, return_parts=False):
    d = np.array([0,0,distance])

    height1 = get_slab_height(in_atoms1)
    height2 = get_slab_height(in_atoms2)

    if cell is None:
        cell = in_atoms1.cell.copy() # + mix*(atoms2.cell.copy() - atoms1.cell.copy())
    cell[2] = np.array([0,0,height1+height2+2*distance])

    atoms1 = sort(in_atoms1.copy())
    atoms1 = in_atoms1.copy()
    atoms1.cell[2] = cell[2]
    atoms1.center(axis=2)

    atoms2 = sort(in_atoms2.copy())
    atoms2.cell[2] = cell[2].copy()
    atoms2.write('atoms2.xyz', format='xyz')
    atoms2.set_cell(np.eye(3), scale_atoms = True)
    atoms2.set_cell(cell.copy(), scale_atoms = True)
    atoms2.write('rescaled_atoms2.xyz', format='xyz')

    adsite1 = get_adsite(atoms1, site=tagadsite1, face='top')['top'][tagadsite1]
    adsite2 = get_adsite(atoms2, site=tagadsite2, face='bottom')['bottom'][tagadsite2]

    translation = adsite1 - adsite2 + d
    print(translation)

    atoms2.translate(translation)
    atoms2.write('translated_rescaled_atoms2.xyz', format='xyz')

    thestack = atoms1.copy()
    thestack.extend(atoms2)

    if return_parts:
        return thestack, atoms1, atoms2
    else:
        return thestack


# Plotting tools
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms

def plotviews(atoms_object,
              rotation1='-90x',
              rotation2='-90x, 45y',
              rotation3='0x'):
  fig = plt.figure()
  ax1 = fig.add_axes((0.1,0.1, 0.3,0.8))
  ax2 = fig.add_axes((0.3,0.1, 0.3,0.8))
  ax3 = fig.add_axes((0.6,0.1, 0.3,0.8))
  plot_atoms(atoms_object,rotation=rotation1, ax=ax1)
  plot_atoms(atoms_object,rotation=rotation2, ax=ax2)
  plot_atoms(atoms_object,rotation=rotation3,ax=ax3)
  return ax1, ax2, ax3
