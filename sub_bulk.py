import os
import sys
from subprocess import run

sub_string = "bsub -n 24 -W 24:00 -R \"rusage[mem=2048]\" -J {0} \"mpirun -n 24 gpaw-python run_bulk.py {0}\""
if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise ValueError("At least provide one material name!")
    else:
        run(sub_string.format(sys.argv[1]), shell=True)
