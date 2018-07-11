import os, os.path
import json
from ase.io import read
import ase.db
from gpaw import GPAW, PW, FermiDirac
from gpaw.response.df import DielectricFunction
from ase.parallel import parprint

def excited(base_dir="./"):
    gs_gpw = os.path.join(base_dir, "gs.gpw")
    es_gpw = os.path.join(base_dir, "es.gpw")
    if os.path.exists(es_gpw):
        parprint("Excited states calculated, will use gpw file directly!")
        return
    
    if not os.path.exists(gs_gpw):
        raise FileNotFoundError("Ground state not calculated!")

    curr_dir = os.path.dirname(os.path.abspath(__file__))
    param_file = os.path.join(curr_dir, "../parameters.json")
    if os.path.exists(param_file):
        params = json.load(open(param_file, "r"))
    else:
        raise FileNotFoundError("no parameter file!")

    calc = GPAW(gs_gpw)
    calc.set(**params["es"])
    calc.get_potential_energy()
    calc.diagonalize_full_hamiltonian(nbands=60)
    calc.write(es_gpw, mode="all")  # full matrix

def polarizability(base_dir="./", mode="df"):
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    es_gpw = os.path.join(base_dir, "es.gpw")
    param_file = os.path.join(curr_dir, "../parameters.json")

    if os.path.exists(param_file):
        params = json.load(open(param_file, "r"))
    else:
        raise FileNotFoundError("no parameter file!")
    
    if not os.path.exists(es_gpw):
        raise FileNotFoundError("Excited state not calculated!")
    
    if mode not in ("df", "tetra"):
        raise ValueError("Mode should be df or tetra")
    
    data_file = os.path.join(base_dir,
                             "polarizability_{}.npz".format(mode))

    if os.path.exists(data_file):
        parprint("Polarizability file exists!")
        return 0
    
    df = DielectricFunction(calc=es_gpw,
                            **params[mode])
    alpha0x, alphax = df.get_polarizability(q_c=[0, 0, 0],
                                            direction="x",
                                            pbc=[True, True, False],
                                            filename=None)
    alpha0y, alphay = df.get_polarizability(q_c=[0, 0, 0],
                                            direction="y",
                                            pbc=[True, True, False],
                                            filename=None)
    alpha0z, alphaz = df.get_polarizability(q_c=[0, 0, 0],
                                            direction="z",
                                            pbc=[True, True, False],
                                            filename=None)
    
    freq = df.get_frequencies()
    data = dict(frequencies=freq,
                alpha_x=alphax,
                alpha_y=alphay,
                alpha_z=alphaz,
                alpha_x0=alpha0x,
                alpha_y0=alpha0y,
                alpha_z0=alpha0z,)
    from ase.parallel import world
    import numpy
    if world.rank == 0:
        numpy.savez_compressed(data_file, **data)
