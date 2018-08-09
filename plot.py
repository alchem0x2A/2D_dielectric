import os, os.path
import sys
import matplotlib.pyplot as plt
import numpy

root_mono  = "/cluster/scratch/ttian/2D/"
root_bulk = "/cluster/scratch/ttian/2D-bulk/"


def plot(name, bulk=False, kind=None):
    
    if bulk:
        root = root_bulk
        pol_file = "eps_df.npz"
    else:
        root = root_mono
        pol_file = "polarizability_df.npz"
    dir_candidates = os.listdir(root)
    if kind is None:
        search = name
    else:
        search = "{}-{}".format(name, kind)
    print(search)
    for directory in dir_candidates:
        if kind is None:
            if search != directory.split("-")[0]:
                continue
        else:
            if search != directory:
                continue
        print("Found {}".format(directory))
        try:
            f = numpy.load(os.path.join(root, directory, pol_file))
        except Exception:
            return
        freq = f["frequencies"]

        if bulk:
            eps_x = f["eps_x"]
            eps_z = f["eps_z"]
            plt.figure()
            plt.plot(freq, eps_x.real, label="xx")
            plt.plot(freq, eps_z.real, label="zz")
            plt.xlabel("$\\hbar \\Omega$ (eV)")
            plt.ylabel("$\\varepsilon$")
            plt.title(directory)
            
            plt.figure()
            plt.plot(freq, eps_x.imag, label="xx")
            plt.plot(freq, eps_z.imag, label="zz")
            plt.xlabel("$\\hbar \\Omega$ (eV)")
            plt.ylabel("$\\varepsilon$")
            plt.title(directory)
        else:
            alpha_x = f["alpha_x"]
            alpha_z = f["alpha_z"]
            plt.figure()
            plt.plot(freq, alpha_x.real, label="xx")
            plt.plot(freq, alpha_z.real, label="zz")
            plt.xlabel("$\\hbar \\Omega$ (eV)")
            plt.ylabel("$\\alpha$ ($\\AA$)")
            plt.title(directory)

            plt.figure()
            plt.plot(freq, alpha_x.imag, label="xx")
            plt.plot(freq, alpha_z.imag, label="zz")
            plt.xlabel("$\\hbar \\Omega$ (eV)")
            plt.ylabel("$\\alpha$ ($\\AA$)")
            plt.title(directory)
    plt.show()
    return 0

if __name__ == "__main__":
    assert len(sys.argv) >= 2
    name = sys.argv[1]
    
    if len(sys.argv) == 2:
        is_bulk = False         # default mono
        kind = None
    if len(sys.argv) == 3:      # bulk or kind?
        opt = sys.argv[2]
        if opt not in ("mono", "bulk"):  # then monolayer with kind
            try:
                kind = sys.argv[2]
            except IndexError:
                kind = None
            is_bulk = False
        else:
            is_bulk = opt == "bulk"
            kind = None
    if len(sys.argv) == 4:
        try:
            kind = sys.argv[2]
        except IndexError:
            kind = None
        opt = sys.argv[3]
        if opt in ("mono", "bulk"):
            is_bulk = opt == "bulk"
        else:
            is_bulk = False
            
    print(is_bulk, kind)
    plot(name, kind=kind, bulk=is_bulk)
    
        
