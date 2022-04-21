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

def save_lattice(spin_lattice, filename, overwrite=False):
    """Saves the lattice as a numpy array string in a .txt file.
    Protocol:
    - A lattice takes up one line
    - first numbers are the shape
    - the delimiter is /
    - then come the array elements
    - For example '2 2/1 2 3 4' would be np.array([[1, 2], [3, 4]])
    :spin_lattice: np.ndarray(int) lattice representation as numpy array
    :filename: str, extensionless filename as a string
    :overwrite: whether to overwrite if the file already exists
    :return: None"""
    
    if overwrite:
        flag = 'w'
    else:
        flag = 'x'
    
    file_str = lattice2str(spin_lattice)
    
    filename += '.txt'
    
    with open(filename, flag) as f:
        
        f.write(file_str)
        f.close()
        
def lattice2str(lattice):
    """Converts a lattice to a line in a file according to convention described in save_lattice.
    :lattice: np.ndarray(int) lattice representation as numpy array
    :return: str, lattice as file string"""
    
    import numpy as np
    
    shape_str = str(
        lattice.shape
    ).strip('()').replace(',', '')
    
    arr_str = np.array_str(
        np.ravel(lattice)
    ).strip('[]')
    
    return shape_str + '/' + arr_str

        
def read_lattice(line):
    """Reads a line of a file containing a lattice.
    :line: str, line of a file in format described under save_lattice
    :return: np.ndarray(int), lattice representation of the saved lattice"""
    
    import numpy as np
    
    split_file_str = line.split('/')
    
    shape = tuple(
        [int(i) for i in split_file_str[0].split()]
    )
    
    array = np.fromstring(
        split_file_str[1],
        dtype=int,
        sep=' '
    ).reshape(shape)
    
    return array

# 18 April 2022: This function is yet to be implemented, needs tests writing
# Two options for implementation:
# - pass as an additional argument to MetropolisRun, magnetisation is calculated during the run (reduced time/space complexity)
# - calculate using saved lattices from a previous run (more flexible, saves time during the run if magnetisation is not desired)
        
def get_mag_pyrochlore(lattice):
    """Finds the magnetisation of a given lattice.
    See MetropolisRun documentation for definitions.
    :lattice: np.ndarray(int), lattice representation of spins
    :return: np.ndarray(float), magnetisation of the lattice"""

    # Gives the directions in space of the pyrochlore spins
    
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
    

           