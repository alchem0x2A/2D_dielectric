import sys
import os, os.path
import subprocess
import numpy

if __name__ == "__main__":
    formula = sys.argv[1]
    kind = sys.argv[2]
    program = "bsub"
    params = ["-n", "4", "-W", "04:00", "-R", "\"rusage[mem=4096]\""]
    run_string = "\"mpirun -n 4 gpaw-python run_doping.py {0} {1} {2:.2f}\""
    for fermi_shift in numpy.linspace(-1.0, 1.0, 41):
        call_cmd = [program] + params + [run_string.format(formula, kind, fermi_shift)]
        print(" ".join(call_cmd))
    choice = input("Should we proceed? [y/n]")
    if choice in ("y", "Y"):
        for fermi_shift in numpy.linspace(-1.0, 1.0, 41):
            call_cmd = [program] + params + [run_string.format(formula, kind, fermi_shift)]
            proc = subprocess.run(" ".join(call_cmd),
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
            print(proc.stdout.decode("utf8"))
    
        
	
