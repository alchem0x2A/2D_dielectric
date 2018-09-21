import sys
import os, os.path
# May need this for the path issue
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from src.build import get_structure
from src.relax import relax
from src.dielectric import excited, permittivity
import shutil
from ase.parallel import paropen, parprint, world, rank, broadcast

def main(formula, root="/cluster/scratch/ttian/2D-bulk/",
         clean=False):
    # candidates = {}
    # if rank == 0:
    # If MX2 then use a larger c
    if any([s in formula for s in ("I", "O", "S", "Se", "Te")]):
        c = 6.0
    else:
        c = 3.0

    candidates = get_structure(formula, c=c)
    # candidates = broadcast(candidates, root=0)
    parprint(rank, candidates)

    # Directory manipulation
    if rank == 0:
        for name in candidates:
            base_dir = os.path.join(root, name)
            if clean:
                shutil.rmtree(base_dir, ignore_errors=True)
                
            if not os.path.exists(base_dir):
                os.makedirs(base_dir)
        
    world.barrier()
    if clean:
        return                  # on all ranks

    # On all ranks
    for name in candidates:
        base_dir = os.path.join(root, name)
        mol = candidates[name]
        # Relaxation and gs
        relax(mol, name=name,
              base_dir=base_dir)
        parprint("Relaxation for {} finished!".format(name))

        excited(base_dir=base_dir)
        parprint("Excitation for {} finished!".format(name))

        permittivity(base_dir=base_dir, mode="df")
        parprint("Permittivity for {} finished!".format(name))

        # polarizability(base_dir=base_dir, mode="tetra")  # 
        # parprint("Polarizability using tetra {} finished!".format(name))
    return 0

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise ValueError("Only 1 parameter is needed!")
    elif len(sys.argv) == 2:
        formula = sys.argv[1]
        main(formula)
    elif len(sys.argv) == 3:
        formula = sys.argv[1]
        if sys.argv[2] == "clean":
            main(formula, clean=True)
        else:
            main(formula)
    else:
        raise ValueError("Parameter ill defined!")
