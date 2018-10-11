import os, os.path
import numpy
import json
from ase.io import read
import ase.db
from gpaw import GPAW, PW, FermiDirac
from ase.parallel import parprint, paropen
from numpy import polyfit
from ase.io.trajectory import Trajectory

import ase.optimize
from ase.constraints import UnitCellFilter, StrainFilter

from ase.data import covalent_radii

from scipy.interpolate import interp1d

def best_d(d, E):
    sp = interp1d(d, E, kind="cubic")
    xx = numpy.linspace(d[0], d[-1], 100)
    yy = sp(xx)
    x_opt = xx[numpy.argmin(yy)]
    y_opt = numpy.min(yy)
    return x_opt, y_opt



def get_thick(atom_row):
    pos = atom_row.positions[:, -1]
    diff = covalent_radii[atom_row.numbers]
    zmax = numpy.max(pos + diff) - numpy.min(pos - diff)
    return zmax

# Relax ingle
def relax(atoms, name="",
          base_dir="./"):
    # If pbc in z direction, use vdW for relaxation as default!
    if atoms.pbc[-1]:           # Do not compare with np.bool_ !
        use_vdW = True
    else:
        use_vdW = False

    
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    param_file = os.path.join(curr_dir, "../parameters.json")
        
    gpw_file = os.path.join(base_dir, "gs.gpw")
    if os.path.exists(gpw_file):
        parprint("Relaxation already done, will use gpw directly!")
        return 0
    if os.path.exists(param_file):
        params = json.load(open(param_file, "r"))
    else:
        raise FileNotFoundError("no parameter file!")

    # vdw = None
    # Change xc to vdW
    # if use_vdW:
        # params["relax"].pop("poissonsolver", None)  # Remove the dipole correction since pbc_z=True
        # params["relax"]["xc"] = {"name": "vdW-DF2",
                                 # "backend": "libvdwxc"}
        # params["relax"]["parallel"] = dict(augment_grids=True)
        # params["gs"].pop("poissonsolver", None)  # Remove the dipole correction since pbc_z=True
        
        
    # calculation asign
    calc = GPAW(**params["relax"])
    
    # atoms.set_calculator(calc)
    traj_filename = os.path.join(base_dir,
                                 "{}_relax.traj".format(name))
    log_filename = os.path.join(base_dir,
                                 "{}_relax.log".format(name))
    
    # Get the min(E) by quadratic fitting
    # c0 = atoms.cell[-1][-1]
    atoms.set_calculator(calc)
    c0 = get_thick(atoms)
    if c0 < 2:
        c0 = 3.3
    else:
        c0 = c0 + 1
       # print(c0)
    n_loops = 30
    c_range = numpy.linspace(0.75 * c0, 1.5 * c0, n_loops)  # run 20 runs
    parprint(c_range)
    E = []
    for c in c_range:
        atoms.cell[-1][-1] = c
        e = atoms.get_potential_energy()  # I'll trust the self-consistent result
        E.append(e)
    fd = paropen(os.path.join(base_dir, "E_vdw.txt"),
                 mode="w")
    [parprint(c, e, file=fd) for c, e in zip(c_range, E)]
    fd.close()

    d, e = best_d(c_range, numpy.array(E))
    parprint("Best distance: {} with energy {}".format(d, e))

    # Calculate the ground state 
    atoms.cell[-1][-1] = d     # Reduce the cell-c
    atoms.set_calculator(calc)  # make sure
    calc.set(**params["gs"])    # Now use PBE only
    
    atoms.get_potential_energy()
    calc.write(gpw_file)
    
    
    
    
