# Class and functions relating to Metropolis cooling

from .lattice import Lattice

class MetropolisRun(Lattice):
    def __init__(self, spin_lattice, site_lattice=None, energy_function=None):
        self.initial_lattice = lattice.spins
        self.initial_lattice = site_lattice.spins
        mean_magnetisation = []
        mean_monopoles = []
        if energy_function is None:
            energy = 0
        else:
            energy = energy_function(self.initial_lattice, site_lattice=self.initial_site_lattice)
        
    def cool(self, temps, N, energy_change_function):
        initial_lattice = self.initial_lattice
        
    def _determine_flip():
        return True
        
    def save_run(self, filename):
        pass
        
        