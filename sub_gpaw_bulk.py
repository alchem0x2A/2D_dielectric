import os, os.path
import sys
from subprocess import run
import ase.db
from check_status import check

sub_string = ("bsub -n 24 -W 48:00 -R \"rusage[mem=2048]\" "
              "-J \"{0}-{1}\" "
              "-o {0}-{1}-%J "
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
    candidates = db.select(selection="gap_gw>0.05")
    # print(len(list(candidates)))
    params = []
    not_run = []
    for mol in candidates:
        res = list(new_db.select(formula=mol.formula,
                                 prototype=mol.prototype))  # exists?
        if len(res) == 0:                                   # has not been run yet
            db_id = new_db.write(mol.toatoms(), mol.key_value_pairs)
        else:
            db_id = res[0].id
        # Update with default, don't cover
        if hasattr(new_db[db_id], "data"):
            new_db.update(db_id, delete_keys=["data"])
        # Key for bulk not stored yet, otherwise do nothing
        if any(key not in new_db.get(db_id).key_value_pairs \
               for key in default_values):
            new_db.update(db_id, **default_values)
        res = check(mol.formula, mol.prototype, data=False)
        if res["step"] != "full":  # not converged calculations
            not_run.append((mol.formula, mol.prototype))

        # All data
        params.append((mol.formula, mol.prototype))
        # print((mol.formula, mol.prototype))
    
    print(len(params))
    print("Not finished calculations:")
    print(not_run)
    if dryrun is not True:
        for n, p in not_run:
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
