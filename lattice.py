# lattice class and functions

class Lattice:
    """A lattice of Ising spins.
    
    --- Definitions ---
    Ising spin: a microscopic magnetic moment constrained to point along a defined axis
    Lattice: a repeating structure in space, i.e. one which is invariant under certain spatial translations (known as lattice vectors).
    A lattice has a unit cell which is the repeating element.
    The unit cell can contain several Ising spins.
    The locations of the Ising spins in the unit cell relative to the origin are collectively called the basis.
    
    --- Lattice representation concept ---
    A numpy array containing the values 1 and -1 (representing spin axis up and spin axis down) represents the states of Ising spins at different locations in the lattice.
    So that the array can represent a general lattice in n dimensions, it must be n+1 dimensional:
    n dimensions for the spatial position of the unit cell which the spin belongs to,
    and one additional to label the position of the spin in its unit cell.
    The spatial coordinates of a spin can be obtained from a mapping of the labels onto the unit cell basis, along with the unit cell's spatial position.
    In future I use the term 'index location' to refer to the array index location of a spin in this representation (as opposed to its actual location in space)."""
    
    def __init__(self, basis, side, seed=None):
        """Initialises a lattice with random spins.
        Dimensionality is inferred from the basis.
        :basis: np.ndarray, an array of n-dim numpy row vectors denoting the locations of spins in a unit cell (note: the coordinates are given as fractions of the unit cell side length)
        :side: int, the side length of the lattice (in number of unit cells)
        :seed: int, for random number generator
        :return: np.ndarray(int), numpy array lattice representation containing random spins"""
        
        import numpy as np
        
        self.basis = basis
        self.side = side
        self.n = basis.shape[0] # number of spins in the basis
        self.D = basis.shape[1] # dimensionality of space
        
        self.spin_lattice = createLattice(self.n, self.D, side, seed=seed)
        
        
    
def createLattice(n, D, side, seed=None):
    """Creates the representation of a D-dimensional lattice with n random spins per unit cell.
    :D: int, dimensionality
    :n: int, number of spins per unit cell
    :side: int, side length of lattice (in number of unit cells)
    :seed: int, for random number generator
    :return: np.ndarray(int), numpy array lattice representation containing random spins
    >>> createLattice(1, 2, 2, seed=42)
    np.array([[[1, -1], [1, 1]]])
    >>> createLattice(2, 3, 1, seed=42)
    np.array([[[[1, -1]]]])"""
    
    import numpy as np
    
    rng = np.random.default_rng(seed=seed)

    array_dimensions = np.append(
        side * np.ones(D, dtype=int), 
        n
    )

    lattice = np.array(
        2 * np.rint(rng.random(array_dimensions)) - 1,
        dtype=int
    )

    return lattice

def save_lattice(spin_lattice, filename):
    """Saves the lattice as a numpy array string in a .txt file.
    The lattice could be recovered by reading the file and passing its contents to np.array().
    :spin_lattice: np.ndarray(int) lattice representation as numpy array
    :filename: str, extensionless filename as a string
    :return: None"""
    
    with open('filename' + '.txt', 'x') as f:
        f.write(
            np.array2string(spin_lattice, seperator=',')
        )
        
        f.close()
        
def get_mag(filename):
    """Finds the magnetisations of a series of lattices saved as part of an execution of MetropolisRun.cool(). See the documentation for this method for more details of the file structure.
    Takes the mean magnetisation across iterations
    :filename: str, name of the file in which the lattice data is saved
    :return: list(np.ndarray(float)), list of mean magnetisations at each temperature"""
    
    mean_mags = []
    
    directions = np.array([
        [1,-1,-1],
        [-1,1,-1],
        [1,1,1],
        [-1,1,-1],
        [-1,1,-1],
        [1,1,1],
        [-1,-1,1],
        [1,1,1],
        [-1,-1,-1],
        [-1,1,1],
        [1,1,-1],
        [1,-1,1],
        [1,1,-1],
        [-1,1,1],
        [1,1,-1],
        [-1,1,1]
    ])
    
    # if a line has [: make an empty list
    # if a subsequent line has [: make another empty list
    # if a subsequent line has numbers: append them to the innermost list
    # if a subsequent line has ]: append the innermost list to the next innermost list
    
    with open(filename, 'r') as f:
        
        next(f) # first line will be //, get rid of it
        
        mags = []
        
        while True:
            try:
                lattice_str = next(f)
                
                if lattice_str != '//':
                    lattice = np.array(list(
                        lattice_str
                    ))
                    
                else:
                    pass
            
            except StopIteration:
                break
            
    
    side = iceLattice.shape[0]
    v = side ** 3

    latticeIter = nditer(iceLattice,op_flags = ['readwrite'], flags = ['multi_index'])
    M = array([0,0,0])
    for spin in latticeIter:
        indexLocation = latticeIter.multi_index
        basisIndex = indexLocation[3]
        M = M + (spin/sqrt(3)) * directions[(basisIndex,...)]
    M = M/v
    return(M)