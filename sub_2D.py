import os
import sys
from subprocess import run

sub_string = "bsub -n 24 -W 24:00 -R \"rusage[mem=2048]\" -J {0}-{1} \"mpirun -n 24 gpaw-python run.py {0} {1}\""
if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise ValueError("At least provide one material name!")
    else:
        run(sub_string.format(sys.argv[1], sys.argv[2]), shell=True)
