import os, os.path
import sys
from subprocess import run
import ase.db


sub_string = ("bsub -n 24 -W 24:00 -R \"rusage[mem=2048]\" "
              "-J \"{0}-{1}\" "
              "\"mpirun gpaw-python run_bulk.py {0} {1}\"")

def sub_job(name, proto):
    run(sub_string.format(name, proto), shell=True)
    return True

default_values = dict(bulk_calculated=False,
                      bulk_L=-1,
                      bulk_eps_x=-1,
                      bulk_eps_y=-1,
                      bulk_eps_z=-1,
                      bulk_gap=-1,
                      bulk_gap_hse=-1)

def main(dryrun=True):
    db_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "./c2db.db")
    if os.path.exists(db_file) is not True:
        print("Please download the c2db file into {}".format(os.path.dirname(db_file)))
        return False
    db = ase.db.connect(db_file)
    new_db = ase.db.connect(os.path.join("bulk.db"))
    candidates = db.select(selection="gap_gw>0.5")
    params = []
    for mol in candidates:
        res = list(new_db.select(unique_id=mol.unique_id))  # exists?
        if len(res) == 0:
            db_id = new_db.write(mol.toatoms(), mol.key_value_pairs)
        else:
            db_id = res[0].id
        # Update with default, don't cover
        new_db.update(db_id, delete_keys=["data"], **default_values)
        params.append((mol.formula, mol.prototype))
        print((mol.formula, mol.prototype))
    print(len(params))
    if dryrun is not True:
        for n, p in params:
            sub_job(n, p)
        return True
    else:
        return True

if __name__ == "__main__":
    choice = input("Do you need to dryrun? [Y] N")
    if choice.lower()[0] == "y":
        main(dryrun=True)
    else:
        main(dryrun=False)
