import numpy
import ase.db
import os, os.path
from gpaw import GPAW


# Check which step for calculation
def check(name, proto,
          root_dir = "/cluster/scratch/ttian/2D-bulk",
          data=True):
    # print(name, proto)
    base_dir = os.path.join(root_dir, "{}-{}".format(name, proto))
    gs_file = os.path.join(base_dir, "gs.gpw")
    es_file = os.path.join(base_dir, "es.gpw")
    eps_file = os.path.join(base_dir, "eps_df.npz")
    res = {"step":None}
    if os.path.exists(base_dir) is not True:
        res["step"] = None
    elif os.path.exists(gs_file) is not True:
        res["step"] =  "failed"
    elif os.path.exists(es_file) is not True:
        res["step"] = "gs"
    elif os.path.exists(eps_file) is not True:
        res["step"] = "eps"
    else:
        res["step"] = "full"
        if data is True:
            f = numpy.load(eps_file)
            if "L" not in f:
                calc = GPAW(gs_file)
                L = calc.get_atoms().cell[-1][-1]
                numpy.savez(eps_file, **{k:v for k, v in f.items()},
                            L=L)
            else:
                L = float(f["L"])  # original format is scalar
            epsx = f["eps_x"][0].real
            epsy = f["eps_y"][0].real
            epsz = f["eps_z"][0].real

            res["data"] = (epsx, epsy, epsz, L)
    return res


# Update to atoms db
def update_data(db, db_id, data):
    epsx, epsy, epsz, L = data
    try:
        db.update(db_id,
                  bulk_calculated=True,
                  bulk_L=L,
                  bulk_eps_x=epsx,
                  bulk_eps_y=epsy,
                  bulk_eps_z=epsz)
    except Exception:
        return False
    return True

def main(db_file="./bulk.db"):
    if os.path.exists(db_file) is not True:
        print("database not valid yet!")
        return False
    db = ase.db.connect(db_file)
    old_db = ase.db.connect("c2db.db")
    stats = {"not started": {"count": 0, "names": []},
             "failed": {"count": 0, "names": []},
             "gs": {"count": 0, "names": []},
             "es": {"count": 0, "names": []},
             "success": {"count": 0, "names": []}}
    candidates = old_db.select(selection="gap_gw>0.5")  # Maybe more robust?
    # print(len(list(candidates)))
    for mol in candidates:
        # print(mol.formula, mol.prototype)
        r = list(db.select(formula=mol.formula,
                           prototype=mol.prototype))  # exists?
        if len(r) == 0:
            continue
        else:
            db_id = r[0].id
        res = check(mol.formula, mol.prototype)
        name = "{}-{}".format(mol.formula, mol.prototype)
        if res["step"] is None:
            stats["not started"]["count"] += 1
            stats["not started"]["names"].append(name)
        elif res["step"] == "failed":
            stats["failed"]["count"] += 1
            stats["failed"]["names"].append(name)
        elif res["step"] == "gs":
            stats["gs"]["count"] += 1
            stats["gs"]["names"].append(name)
        elif res["step"] == "es":
            stats["es"]["count"] += 1
            stats["es"]["names"].append(name)
        elif res["step"] == "full":
            stats["success"]["count"] += 1
            stats["success"]["names"].append(name)
            # update
            # print(res["data"])
            print(update_data(db, db_id, res["data"]))
        else:
            pass                # ???
    
    for key in ("not started", "failed", "gs", "es", "success"):
        print("{0}, Counts: {1}".format(key, stats[key]["count"]))
    print(stats["failed"]["names"])
    return    

if __name__ == "__main__":
    main()
