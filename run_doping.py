import sys
import os, os.path
# May need this for the path issue
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from src.build import get_structure
from src.relax import relax
from src.polarizability import excited, polarizability
import shutil
from ase.parallel import paropen, parprint, world, rank, broadcast

def main(formula, kind, fermi_shift=0,
         root="/cluster/scratch/ttian/2D/",
         clean=False):
    # candidates = {}
    # if rank == 0:
    # candidates = get_structure(formula)
    # candidates = broadcast(candidates, root=0)
    # parprint(rank, candidates)
    name = "{}-{}".format(formula, kind)
    base_dir = os.path.join(root, name)

    # Directory manipulation
    if rank == 0:
        if clean:
            shutil.rmtree(base_dir, ignore_errors=True)
                
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)
        
    world.barrier()
    if clean:
        return                  # on all ranks

    # On all ranks
        # Relaxation and gs
        # relax(mol, name=name,
              # base_dir=base_dir)
        # parprint("Relaxation for {} finished!".format(name))

        # excited(base_dir=base_dir)
        # parprint("Excitation for {} finished!".format(name))

    polarizability(base_dir=base_dir, mode="df",
                   fermi_shift=fermi_shift)  # add fermilevel
    parprint("Polarizability for {} with Fermi Shift {} finished!".format(name, fermi_shift))

        # polarizability(base_dir=base_dir, mode="tetra")  # 
        # parprint("Polarizability using tetra {} finished!".format(name))
    return 0

if __name__ == "__main__":
    if len(sys.argv) < 4:
        raise ValueError("[formula] [type] [fermi shift]")
    elif len(sys.argv) == 4:
        formula = sys.argv[1]
        kind = sys.argv[2]
        fermi_shift = float(sys.argv[3])
        main(formula, kind, fermi_shift)
    else:
        raise ValueError("Parameter ill defined!")
