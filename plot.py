import os, os.path
import sys
import matplotlib.pyplot as plt
import numpy

root = "/cluster/scratch/ttian/2D/"
pol_file = "polarizability_df.npz"

def plot(name, kind=None):
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
        alpha_x = f["alpha_x"]
        alpha_z = f["alpha_z"]
        plt.figure()
        plt.plot(freq, alpha_x.real, label="xx")
        plt.plot(freq, alpha_z.real, label="zz")
        plt.xlabel("$\\hbar \\Omega$ (eV)")
        plt.ylabel("$\\alpha$ ($\\AA$)")
        plt.title(directory)
    plt.show()
    return 0

if __name__ == "__main__":
    assert len(sys.argv) >= 2
    name = sys.argv[1]
    try:
        kind = sys.argv[2]
    except IndexError:
        kind = None
    plot(name, kind)
    
        
