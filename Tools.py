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
                print(theatoms.info['adatom'])
    else:
        raise ValueError('no info property')


def get_bridge_position(theatoms, thelayer):
    positions = theatoms.positions[thelayer]
    theys = np.unique(positions[:, 1])
    equaly = [positions[positions[:, 1] == they].mean(axis=0) for they in theys]
    print(equaly)
    return positions[equaly] 


def get_top_positions(atoms, top_layer):
    return atoms.positions[top_layer]


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






