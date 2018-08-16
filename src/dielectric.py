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

def permittivity(base_dir="./", mode="df"):
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    es_gpw = os.path.join(base_dir, "es.gpw")
    param_file = os.path.join(curr_dir, "../parameters.json")

    if os.path.exists(param_file):
        params = json.load(open(param_file, "r"))
    else:
        raise FileNotFoundError("no parameter file!")

    # No reason to use 2D truncation for dielectric
    params["df"].pop("truncation", None)
    params["tetra"].pop("truncation", None)
    
    if not os.path.exists(es_gpw):
        raise FileNotFoundError("Excited state not calculated!")
    
    if mode not in ("df", "tetra"):
        raise ValueError("Mode should be df or tetra")
    
    data_file = os.path.join(base_dir,
                             "eps_{}.npz".format(mode))

    if os.path.exists(data_file):
        parprint("Polarizability file exists!")
        return 0
    
    df = DielectricFunction(calc=es_gpw,
                            **params[mode])
    eps0x, epsx = df.get_dielectric_function(q_c=[0, 0, 0],
                                            direction="x",
                                            filename=None)
    eps0y, epsy = df.get_dielectric_function(q_c=[0, 0, 0],
                                            direction="y",
                                            filename=None)
    eps0z, epsz = df.get_dielectric_function(q_c=[0, 0, 0],
                                            direction="z",
                                            filename=None)
    
    freq = df.get_frequencies()
    data = dict(frequencies=freq,
                eps_x=epsx,
                eps_y=epsy,
                eps_z=epsz,
                eps_x0=eps0x,
                eps_y0=eps0y,
                eps_z0=eps0z,)
    from ase.parallel import world
    import numpy
    if world.rank == 0:
        numpy.savez_compressed(data_file, **data)
