import os, os.path
import numpy
import json
from ase.io import read
import ase.db
from gpaw import GPAW, PW, FermiDirac
from ase.parallel import parprint

import time
import numpy as np
import ase.optimize
import gpaw.mpi as mpi
from ase.constraints import UnitCellFilter

from gpaw.xc.exx import EXX
from gpaw.xc.tools import vxc
from gpaw.spinorbit import get_spinorbit_eigenvalues as get_soc_eigs

def runhse(base_dir):
    hse(base_dir)
    mpi.world.barrier()
    hse_spinorbit(base_dir)

def get_kpts_size(atoms, density):
    """trying to get a reasonable monkhorst size which hits high
    symmetry points
    """
    from gpaw.kpt_descriptor import kpts2sizeandoffsets as k2so
    size, offset = k2so(atoms=atoms, density=density)
    size[2] = 1
    for i in range(2):
        if size[i] % 6 != 0:
            size[i] = 6 * (size[i] // 6 + 1)

    kpts = {'size': size, 'gamma': True}
    return kpts

# calculate band gap by GGA
def hse(base_dir="./",
        kdens=6.0):
    emptybands = 20
    convbands = 10
    # If pbc in z direction, use vdW for relaxation as default!
    # if atoms.pbc[-1]:           # Do not compare with np.bool_ !
        # use_vdW = True
    # else:
        # use_vdW = False
    
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    param_file = os.path.join(curr_dir, "../parameters.json")
        
    gpw_file = os.path.join(base_dir, "gs.gpw")
    hse_file = os.path.join(base_dir, "hse.gpw")
    hse_nowfs_file = os.path.join(base_dir, "hse_nowfs.gpw")
    hse_eigen_file = os.path.join(base_dir, "hse_eigenvalues.npz")
    if not os.path.exists(gpw_file):
        parprint("No ground state calculation? Exit...")
        return 0
    if os.path.exists(param_file):
        params = json.load(open(param_file, "r"))
    else:
        raise FileNotFoundError("no parameter file!")

    if not os.path.exists(hse_file):
        calc = GPAW(gpw_file)       # reload the calculation
        atoms = calc.get_atoms()
        kpts = get_kpts_size(atoms, kdens)
        calc.set(nbands=-emptybands,
                 fixdensity=True,
                 kpts=kpts,
                 convergence={'bands': -convbands})
        calc.get_potential_energy()
        calc.write(hse_file, 'all')
        calc.write(hse_nowfs_file)  # no wavefunction
    
    mpi.world.barrier()
    time.sleep(10)  # is this needed?
    calc = GPAW(hse_file, txt=None)
    ns = calc.get_number_of_spins()
    nk = len(calc.get_ibz_k_points())
    nb = calc.get_number_of_bands()
    vxc_pbe_skn = vxc(calc, 'PBE')
    vxc_pbe_nsk = numpy.ascontiguousarray(vxc_pbe_skn.transpose(2, 0, 1))
    vxc_pbe_nsk = calc.wfs.bd.collect(vxc_pbe_nsk, broadcast=True)
    vxc_pbe_skn = vxc_pbe_nsk.transpose(1, 2, 0)[:, :, :-convbands]
    e_pbe_skn = np.zeros((ns, nk, nb))
    for s in range(ns):
        for k in range(nk):
            e_pbe_skn[s, k, :] = calc.get_eigenvalues(spin=s, kpt=k)

    e_pbe_skn = e_pbe_skn[:, :, :-convbands]
    hse_calc = EXX(hse_file, xc='HSE06', bands=[0, nb - convbands])
    hse_calc.calculate()
    vxc_hse_skn = hse_calc.get_eigenvalue_contributions()
    e_hse_skn = e_pbe_skn - vxc_pbe_skn + vxc_hse_skn
    ranks = [0]
    if mpi.world.rank in ranks:
        dct = dict(vxc_hse_skn=vxc_hse_skn,
                   e_pbe_skn=e_pbe_skn,
                   vxc_pbe_skn=vxc_pbe_skn,
                   e_hse_skn=e_hse_skn)
        with open(hse_eigen_file, 'wb') as f:
                numpy.savez(f, **dct)
    parprint("Single HSE06 finished!")
    
def hse_spinorbit(base_dir="./"):
    hse_file = os.path.join(base_dir, "hse.gpw")
    hse_nowfs_file = os.path.join(base_dir, "hse_nowfs.gpw")
    hse_eigen_file = os.path.join(base_dir, "hse_eigenvalues.npz")
    hse_eigen_soc_file = os.path.join(base_dir, 'hse_eigenvalues_soc.npz')
    if not os.path.isfile(hse_eigen_file):
        return
    if not os.path.isfile(hse_eigen_file):
        return

    ranks = [0]
    comm = mpi.world.new_communicator(ranks)
    if mpi.world.rank in ranks:
        calc = GPAW(hse_nowfs_file, communicator=comm, txt=None)
        with open(hse_eigen_file, 'rb') as fd:
            dct = dict(np.load(fd))

        e_skn = dct.get('e_hse_skn')
        dct = {}
        e_mk, s_kvm = get_soc_eigs(calc, gw_kn=e_skn, return_spin=True,
                                   bands=np.arange(e_skn.shape[2]))
        dct['e_hse_mk'] = e_mk
        dct['s_hse_mk'] = s_kvm[:, 2, :].transpose()
        with open(hse_eigen_soc_file, 'wb') as fd:
            np.savez(fd, **dct)
    parprint("SOC HSE06 finished!")
    
