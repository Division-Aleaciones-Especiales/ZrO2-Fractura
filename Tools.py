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
    ad_pos = [theatoms.info['adatom'][theface][thesite]+[0,0,d]]
    adatom_inface_insite = Atoms('H', positions=ad_pos, pbc=True, cell=theatoms.cell.copy())
    adstruc = theatoms.copy()
    adstruc.extend(adatom_inface_insite)
    write_vasp(f'{name}_{theface}face_ad{thesite}.vasp', adstruc,direct=True, sort=True, vasp5=True, wrap=True)
    return adstruc

def get_slab_height(theatoms):
    return theatoms.positions[:,2].max() - theatoms.positions[:,2].min()

def stack(in_atoms1, in_atoms2, adsite1, adsite2, distance, mix=0.5, cell = None):

    atoms1 = in_atoms1.copy()
    atoms2 = in_atoms2.copy()

    if cell is None:
        cell = atoms1.cell.copy() + mix*(atoms2.cell.copy() - atoms1.cell.copy())
    height1 = get_slab_height(in_atoms1)
    height2 = get_slab_height(in_atoms2)
    cell[2] = np.array([0,0,height1+height2+2*distance])

    atoms1.set_cell(np.array([ atoms1.cell[0].copy(), atoms1.cell[1].copy(), cell[2] ]))
    atoms2.set_cell(np.array([ atoms2.cell[0].copy(), atoms2.cell[1].copy(), cell[2] ]))

    _ = [theatoms.set_cell([[1,0,0],[0,1,0],[0,0,1]], scale_atoms=True) for theatoms in [atoms1, atoms2]]


    t1 = atoms1.positions.min(axis=0)
    atoms1.translate(-t1)
    
    t2 = atoms2.positions.min(axis=0)
    atoms2.translate(-t2+[0,0,( height1+distance)/cell[2][2]])

    atoms1.extend(atoms2)

    atoms1.set_cell(cell, scale_atoms = True)
    atoms1.pbc = True

    return atoms1, atoms2







