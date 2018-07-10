import ase.db
from ase.atoms import Atoms
from ase.parallel import parprint, world, broadcast
import os, os.path

def get_structure(formula, base_dir="./", **filters):
    db_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../c2db.db")
    # Serial version, only on rank 0
    candidates = {}
    # if world.rank == 0:
    db = ase.db.connect(db_file)
    res = db.select(formula=formula, **filters)
    for mol in res:
        symbol = mol.formula
        pos = mol.positions
        cell = mol.cell
        pbc = mol.pbc
        name = "{}-{}".format(symbol,
                              mol.prototype)
            # name = os.path.join(os.path.abspath(base_dir),
            # "{}-{}.traj".format(symbol, mol.prototype))
        atoms = Atoms(symbol, positions=pos,
                      cell=cell, pbc=pbc)
        candidates[name] = atoms
    return candidates
        # atoms.write(name)

if __name__ == "__main__":
    import os, sys
    if len(sys.argv) < 2:
        raise ValueError("not enough parameters")
    else:
        formula = sys.argv[1]
        print(get_structure(formula))
