# Class and functions relating to Metropolis cooling

from .lattice import Lattice

class MetropolisRun(Lattice):
    """A lattice object to which the Metropolis algorithm can be applied. Child class of Lattice.
    
    --- Definitions ---
    Boltzmann factor: exp(- E/hbar T) where E is the energy, hbar is the reduced Planck's constant and T is the temperature
    Metropolis algorithm: an algorithm used to compute thermodynamic properties of objects using their micro-states.
    A macroscopic physical system can have many microstates (a specific combination of miscroscopic properties of the system's constituents).
    Its thermodynamic properties depend on the energetic structure of its microstates and the time spent in each one, which is governed by the Boltzmann distribution at a certain temperature T.
    The microstates with their associated probabilities, and the derived macroscopic properties, are collectively known as the canonical ensemble.
    The Metropolis algorithm computes thermodynamic properties as follows:
    - For each microscopic component of the system, propose a change of state
    - Calculate the energy change resulting from that change of state
    - If the energy change is negative (favourable), apply the change
    - If the energy change is positive (unfavourable), apply the change with a probability of exp(- delta E/hbar T), i.e. a Boltzmann factor
    - If this process is repeated a large number of times the thermodynamic properties converge to those of the canonical ensemble.
    Magnetic moment: A vector giving the axis orientation of a magnetic dipole (north-south pole pair). The magnitude is the magnetic dipole strength.
    Magnetisation: Overall magnetic moment of a sample per unit volume. A vector quantity.
    Magnetic monopoles: individual magnetic north poles or south poles. These cannot exist in free space but can exist in certain exotic magnetic materials at low temperatures.
    
    For this project, the only micro-properties of interest are the spins of the ions sitting on lattice sites.
    A spin typically has a magnetic moment.
    """
    def __init__(self, Lattice, energy_function=None, energy_change_function=None, initial_site_lattice=None, **kwargs):
        """Initialises the test lattice sample.
        There may or may not be a site_lattice involved. Some systems require this paradigm to be able to calculate the energy or the energy change.
        Note that
        - There must be either an energy function or an energy change function
        - If there is no energy function, the initial energy will be assumed to be zero
        - If the energy change function requires a site lattice, there must be either an energy function which provides an initial site lattice
        an initial site lattice which is passed manually. If both are present the result of the energy change function will take precedent.
        - If a site lattice is needed but not provided, the error handling depends on the implementation of the energy_change_function but is likely to be a TypeError relating to a forbidden operation on a NoneType object. Ideally the energy_change_function implementation would raise a more informative error if the site_lattice is None and is not supposed to be.
        
        :Lattice: a Lattice object representing the spins
        :energy_function: a function which takes a lattice and any necessary parameters, and returns a generator which yields the energy, followed by a site lattice if applicable
        :energy_change_function: a function which takes a lattice and index location. It may also take a site lattice and any necessary parameters. It returns the change in energy due to changing one spin. This can drastically reduce the time complexity for systems with short range interations.
        :initial_site_lattice: an initial site lattice to pass, in case an energy function is not provided but a site lattice is required
        :kwargs: any parameters required for the energy functions should be passed as keyword arguments
        :return: None"""
        
        self.basis = Lattice.basis
        self.side = Lattice.side
        self.n = Lattice.basis.shape[0] # number of spins in the basis
        self.D = Lattice.basis.shape[1] # dimensionality of space
        self.spin_lattice = Lattice.spin_lattice
        self.initial_spin_lattice = self.spin_lattice.copy()
        
        self.temp = None # the last temperature at which the lattice was equilibrated. None if the lattice has not been equilibrated.
        
        self.energies = []
        
        self.energy_function = energy_function
        self.energy_change_function = energy_change_function
        
        self.params = kwargs

        if energy_function is None and energy_change_function is None:
            raise ValueError('There must be either an energy_function or an energy_change_function')
            
        elif energy_function is None:
            self.initial_energy = 0
            
        else:
            energy_gen = energy_function(self.initial_spin_lattice, self.params)
            self.initial_energy = next(energy_gen)
            try:
                self.initial_site_lattice = next(energy_gen)
            except StopIteration:
                self.initial_site_lattice = initial_site_lattice
        
        self.energy = self.initial_energy
        self.params['site_lattice'] = self.initial_site_lattice.copy()
        
        
    def cool(self, temps, N, lattice_filename=None, seed=None):
        """Cool a lattice using the Metropolis algorithm.
        To find the ground state of a lattice, decrease the temperature gradually, whilst equilibrating the lattice at each step using the Metropolis algorithm
        
        --- Definitions ---
        Ground state: lowest energy microstate of a system
        
        :temps: list(float), temperatures through which to cool the lattice
        :N: int, number of Metropolis algorithm iterations at each temperature
        :lattice_filename: str, an extensionless filename to save the lattices from the run in. If None, the lattices will not be saved. One might want to save the lattices to calculate other properties like the magnetisation. The lattices will be saved in the form of a text file, where each line contains a lattice, a new iteration is denoted by / and a new temperature step is denoted by //. See example_lattice_file.txt for an example.
        :seed: int, seed for random number generator
        :return: None"""
        
        # prepare lattice text file
        
        if lattice_filename is not None:
            file = open(lattice_filename + '.txt', 'w')
        
        import numpy as np

        for T in temps:
            
            if lattice_filename is not None: # add delimiters to text file: new temp step
                file.write('//')
                file.write('\n')
            
            energies  = []
            
            for i in range(0, N):
                
                if lattice_filename is not None: # add new line to text file: new iteration
                    file.write('\n')
                
                E = self.energy
                
                # Index over lattice of spins

                lattice_iter = np.nditer(self.spin_lattice, op_flags = ['readwrite'], flags = ['multi_index'])
                
                for spin in lattice_iter:
                    
                    spin_loc = np.array(lattice_iter.multi_index)
                    
                    energy_change_gen = self.energy_change_function(self.spin_lattice, spin_loc, self.params)
                    energy_change = next(energy_change_gen)
                    flip_bool = _determine_flip(energy_change, T)
                    
                    if flip_bool:
                        spin[...] = -spin[...]
                        E += energy_change
                        try:
                            self.params['site_lattice'] = next(energy_change_gen)
                        except StopIteration:
                            pass
                        
                    else:
                        pass
                    
                energies.append(E)
                
                if lattice_filename is not None:
                    file.write(np.array2string(self.spin_lattice, separator=','))
                    file.write('\n')
                
            self.energy = np.mean(energies)
            self.energies.append(self.energy)
            
            self.temp = T
            
        if lattice_filename is not None:
            file.close()
            
    def reset(self):
        """Resets lattice to the initial state and wipes run data
        :return: None"""
        
        self.spin_lattice = self.initial_spin_lattice
        self.temp = None
        self.params['site_lattice'] = self.initial_site_lattice
        self.energy = self.initial_energy
        self.energies = []
        
def _determine_flip(energy_change, T, seed=None):
    """Determine whether a spin should be flipped.
    :energy_change: float, change in energy due to flipping a spin
    :seed: int, seed for random number generator
    :return: bool, whether the spin should be flipped"""
    
    import numpy as np
    
    rng = np.random.default_rng(seed=seed)

    a = rng.random()

    if energy_change <= 0:
        return True

    elif a < np.exp(-energy_change/T):
        return True

    else:
        return False
        