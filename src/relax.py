import os, os.path
import numpy
import json
from ase.io import read
import ase.db
from gpaw import GPAW, PW, FermiDirac
from ase.parallel import parprint

import time
import numpy as np
import ase.optimize
from ase.constraints import UnitCellFilter


class QuasiNewton:
    def __init__(self, atoms, logfile=None, trajectory=None):
        self.atoms = atoms
        self._logfile = logfile
        self.trajectory = trajectory

    def run(self, fmax, smax, smask=None, emin=-np.inf):
        self.smax = smax
        self.smask = smask
        self.emin = emin
        uf = UnitCellFilter(self.atoms, mask=smask)

        self.opt = ase.optimize.BFGS(uf,
                                     logfile=self._logfile,
                                     trajectory=self.trajectory)
        self.opt.log = self.log
        self.opt.converged = self.converged
        self.force_consistent = self.opt.force_consistent
        self.step0 = self.opt.step
        self.opt.step = self.step
        self.opt.run(fmax)

    def step(self, f):
        m = self.atoms.get_magnetic_moments()
        self.step0(f)
        self.atoms.set_initial_magnetic_moments(m)

    def converged(self, forces):
        if self.atoms.get_potential_energy() < self.emin:
            return True
        if (forces[:-3]**2).sum(axis=1).max() > self.opt.fmax**2:
            return False
        stress = self.atoms.get_stress() * self.smask
        return abs(stress).max() < self.smax

    @property
    def nsteps(self):
        return self.opt.nsteps

    @property
    def logfile(self):
        return self.opt.logfile

    def log(self, forces):
        fmax = (forces[:-3]**2).sum(axis=1).max()**0.5
        stress = self.atoms.get_stress() * self.smask
        smax = abs(stress).max()
        e = self.atoms.get_potential_energy(
            force_consistent=self.force_consistent)
        m = self.atoms.get_magnetic_moment()
        ms = self.atoms.get_magnetic_moments()
        T = time.localtime()
        if self.logfile is not None:
            name = self.__class__.__name__
            if self.nsteps == 0:
                self.logfile.write(
                    '%s  %4s %8s %15s  %12s %12s %8s %8s\n' %
                    (' ' * len(name), 'Step',
                     'Time', 'Energy', 'fmax', 'smax', 'totmm', 'maxmm'))
                if self.force_consistent:
                    self.logfile.write(
                        '*Force-consistent energies used in optimization.\n')
            self.logfile.write(
                '%s:  %3d %02d:%02d:%02d %15.6f%1s %12.4f %12.4f '
                '%8.4f %8.4f\n' %
                (name, self.nsteps, T[3], T[4], T[5], e,
                 {1: '*', 0: ''}[self.force_consistent],
                 fmax, smax, m, abs(ms).max()))
            self.logfile.flush()



# Relax ingle
def relax(atoms, name="",
          base_dir="./",
          smax=2e-4):
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    param_file = os.path.join(curr_dir, "../parameters.json")
    gpw_file = os.path.join(base_dir, "gs.gpw")
    if os.path.exists(gpw_file):
        parprint("Relaxation already done, will use gpw directly!")
        return 0
    if os.path.exists(param_file):
        params = json.load(open(param_file, "r"))
    else:
        raise FileNotFoundError("no parameter file!")
    
    # calculation asign
    calc = GPAW(**params["relax"])
    atoms.set_calculator(calc)
    traj_filename = os.path.join(base_dir,
                                 "{}_relax.traj".format(name))
    log_filename = os.path.join(base_dir,
                                 "{}_relax.log".format(name))
    opt = QuasiNewton(atoms,
                      trajectory=traj_filename,
                      logfile=log_filename)
    mask = [1, 1, 0, 0, 0, 0]   # only relax for xy direction
    opt.run(fmax=0.01, smax=smax, smask=mask)
    
    # Calculate the ground state 
    calc.set(**params["gs"])
    atoms.get_potential_energy()
    calc.write(gpw_file)
    
    
    
    
