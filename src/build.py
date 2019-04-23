import ase.db
from ase.atoms import Atoms
from ase.parallel import parprint, world, broadcast
import os, os.path

def get_structure(formula, prototype=None, base_dir="./", c=None, **filters):
    db_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../c2db.db")
    # Serial version, only on rank 0
    candidates = {}
    # if world.rank == 0:
    db = ase.db.connect(db_file)
    if prototype is None:
        res = db.select(formula=formula, **filters)
    else:
        res = db.select(formula=formula, prototype=prototype)
    for mol in res:
        symbol = mol.formula
        pos = mol.positions
        cell = mol.cell
        pbc = mol.pbc
        # Change distance
        if (c is not None) and (isinstance(c, (float, int))):
            cell.setflags(write=True)
            cell[-1][-1] = c
            pbc = (True, True, True)  # Use full periodic
        name = "{}-{}".format(symbol,
                              mol.prototype)
            # name = os.path.join(os.path.abspath(base_dir),
            # "{}-{}.traj".format(symbol, mol.prototype))
        atoms = Atoms(symbol, positions=pos,
                      cell=cell, pbc=pbc)
        candidates[name] = atoms

        if (c is not None) and (isinstance(c, (float, int))):
            atoms.center()      # center the atoms, although not really needed
    return candidates
        # atoms.write(name)

if __name__ == "__main__":
    import os, sys
    if len(sys.argv) < 2:
        raise ValueError("not enough parameters")
    else:
        formula = sys.argv[1]
        print(get_structure(formula))
