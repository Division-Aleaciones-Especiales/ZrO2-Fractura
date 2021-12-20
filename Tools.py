import pdb
import os
from ase import Atoms
from ase.geometry import get_layers
import numpy as np
import sys
sys.path.insert(0, '/data/git/ase/')
from ase.build.surfaces_with_termination import atom_index_in_top, atom_index_in_bottom

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


def get_bridge_position(theatoms, thelayer):
    positions = theatoms.positions[thelayer]
    theys = np.unique(positions[:, 1])
    equaly = [positions[positions[:, 1] == they].mean(axis=0) for they in theys]
    print(equaly)
    return positions[equaly] 


def get_top_positions(atoms, top_layer):
    return atoms.positions[top_layer[0]]


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
    elif site == None:
        info['adatom'][face].update( {'hollow': atoms.positions[layer].mean(axis=0)})
        info['adatom'][face].update({'top': get_top_positions(atoms, layer)})

    return info['adatom']


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

def stack(in_atoms1, in_atoms2, adsite1, adsite2, distance, mix=0.5, cell = None):

    d = np.array([0,0,distance])

    height1 = get_slab_height(in_atoms1)
    height2 = get_slab_height(in_atoms2)
    if cell is None:
        cell = in_atoms1.cell.copy() # + mix*(atoms2.cell.copy() - atoms1.cell.copy())
    cell[2] = np.array([0,0,height1+height2+2*distance])

    atoms1 = in_atoms1.copy()
    atoms1.center()
    correction = (atoms1.get_positions() - in_atoms1.get_positions()).mean(axis=0)

    atoms2 = in_atoms2.copy()
    atoms2.cell[2] = cell[2].copy()
    atoms2.set_cell(np.eye(3), scale_atoms = True)
    atoms2.set_cell(cell.copy(), scale_atoms = True)

    try:
        atoms1.info['adatom']['top'][adsite1]+=correction
        adsite1 = atoms1.info['adatom']['top'][adsite1]
    except KeyError:
        print('in_atoms1 do not have adatom . falling back to max heght')
        adsite1 = atoms1.get_positions().max(axis=0)
    try:
        adsite2 = atoms2.info['adatom']['bottom'][adsite2]
    except KeyError:
        print('in_atoms2 do not have adatom . falling back to min height')
        adsite2 = atoms1.get_positions().min(axis=0)
    translation = adsite1 - adsite2 + d

    atoms2.translate(translation)

    thestack = atoms1.copy()
    thestack.cell /= np.linalg.norm(thestack.cell)

    thestack.extend(atoms2)

    newcell = in_atoms1.cell.copy()

    newcell[2]=[0,0,height1+height2+2*distance]

    thestack.set_cell(newcell)

    return thestack







