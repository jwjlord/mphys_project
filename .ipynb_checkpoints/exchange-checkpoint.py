# exchange energy functions

def pyrochlore_exchange_energy(lattice, params):
    """Calculate the exchange energy of the pyrochlore lattice.
    
    --- Definitions ---
    Pyrochlore lattice: lattice of corner sharing tetrahedra
    Diamond lattice: structure exhibited by carbon atoms in diamond form.
    It happens that the centres of tetrahedra in a pyrochlore lattice form a diamond lattice.
    Exchange energy: an energy arising from interaction between nearby spins, depending on their relative alignment, caused by quantum mechanical effects.
    
    :lattice: np.ndarray(int), a lattice representation of the spins in the pyrochlore lattice
    :params: dict(str: float), must contain the exchange energy constant J as the key with a float as its value
    :yield: float, total exchange energy of the lattice
    :yield: np.ndarray(int), a diamond-like site lattice containing integers from 0 to 4, denoting how many spins point out of each tetrahedron.
    The latter is required because the exchange energy in the pyrochlore lattice depends on how many spins point into or out of each tetrahedron.
    
    >>>> lattice_energy = pyrochlore_exchange_energy(np.array([[[[1, -1, 1, 1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1]]]]), {'J': 1})
    >>>> next(lattice_energy)
    10
    >>>> next(lattice_energy)
    np.array([[[[3, 3, 3, 2, 1, 1, 3, 0]]]])"""
    
    import numpy as np
    
    J = params['J']
    
    side = lattice.shape[0]
    
    # The following key maps sites in the site lattice (which has the diamond structure) to spins in the pyrochlore lattice.
    # Given a pyrochlore lattice of spins, one uses it to deduce the number of spins which point into or out of each tetrahedron
    # Key is 4 x 4 x 8 array where each 4x4 layer contains the index locations of the spins interacting with each site
    # In the 4x4 layers the rows are the spin index locations (3 lattice + basis)
    # Convention adopted: for each site the first two spins point out of the tetrahedron when positive and the second two point into it.
    # Explanation of convention: each tetrahedron has four corners, with a spin sitting on each, but there are two spins per tetrahedron overall.
    
    key = np.array([
        [[0,0,0,0],[0,0,0,1],[-1,-1,0,14],[0,0,0,8]],
        [[0,0,0,2],[0,0,0,3],[0,0,0,9],[-1,0,0,12]],
        [[0,0,0,4],[0,0,0,5],[0,-1,0,10],[0,-1,0,15]],
        [[0,0,0,6],[0,0,0,7],[0,0,-1,11],[0,0,-1,13]],
        [[0,0,0,8],[0,0,0,9],[0,0,0,4],[0,0,0,6]],
        [[0,0,0,10],[0,0,0,11],[0,0,0,2],[0,1,1,0]],
        [[0,0,0,12],[0,0,0,13],[0,0,0,5],[1,0,1,1]],
        [[0,0,0,14],[0,0,0,15],[0,0,0,7],[1,0,0,3]]
                ])
    
    E = 0
    
    site_lattice = np.zeros((side,side,side,8), dtype='int') # initialise a site lattice to iterate over
    
    site_lattice_out = site_lattice.copy() # and one to return
    
    lattice_iter = np.nditer(site_lattice, op_flags = ['readwrite'], flags = ['multi_index'])
    
    for site in lattice_iter:
        
        site_index = np.array(lattice_iter.multi_index) # Get index location of site in question

        spin_loc = site_to_lattice(site_index, key)
        
        # Periodic boundary conditions
        
        spin_loc_per = [
            impose_pbc(loc, side)
            for loc in spin_loc
        ]
        
        # Get the values of the spins
            
        spins = [
            lattice[tuple(loc)]
            for loc in spin_loc_per
        ]
        
        # Populate the lattice of sites
        
        out = 0.5 * (spins[0] + spins[1] - spins[2] - spins[3]) + 2
        
        site_lattice_out[tuple(site_index)] = out
        
        # Tally contributions to exchange energy

        E += pyrochlore_site_energy(out, J)
        
        # note: there is a trade-off here between flexibility in the pyrochlore_site_energy function and time complexity.
        # in projects like these it's quite useful to be able to change how you calculate the energy
            
    yield E
    yield site_lattice_out
    

def pyrochlore_exchange_energy_change(lattice, spin_loc, params):
    """Find the change in the exchange energy of the lattice due to changing one spin.
    :lattice: np.ndarray(int), a lattice representation of the spins (on tetrahedron corners) in the pyrochlore lattice
    :spin_loc: np.ndarray(int), index loction of the spin which is to be changed
    :params: dict(str: float, str: np.ndarray(int)), must contain the exchange energy constant J as the key with a float as its value, and site_lattice, a lattice representation of the sites (tetrahedron centres) in the pyrochlore lattice
    :yield: float, energy change when the spin is flipped
    :yield: np.ndarray(int), updated lattice represetation of the sites
    
    >>>> lattice = np.array([[[[1, -1, 1, 1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1]]]])
    >>>> site_lattice = np.array([[[[3, 3, 3, 2, 1, 1, 3, 0]]]])
    >>>> energy_change = pyrochlore_exchange_energy_change(lattice, np.array([0, 0, 0, 1]), {'J': 1, 'site_lattice': site_lattice})
    >>>> next(energy_change)
    2
    >>>> next(energy_change)
    np.array([[[[4, 3, 3, 2, 1, 1, 2, 0]]]])"""
    
    # Preliminaries
    
    import numpy as np
    
    J = params['J']
    site_lattice = params['site_lattice']
    
    side = lattice.shape[0]
    
    spin = lattice[tuple(spin_loc)]
    
    # this key is for the inverse mapping associated with the key we saw earlier
    # it shows the two sites which are associated with each spin
    
    key = np.array([
        [[0,0,0,0],[0,-1,-1,5]],
        [[0,0,0,0],[-1,0,-1,6]],
        [[0,0,0,1],[0,0,0,5]],
        [[0,0,0,1],[-1,0,0,7]],
        [[0,0,0,2],[0,0,0,4]],
        [[0,0,0,2],[0,0,0,6]],
        [[0,0,0,3],[0,0,0,4]],
        [[0,0,0,3],[0,0,0,7]],
        [[0,0,0,4],[0,0,0,0]],
        [[0,0,0,4],[0,0,0,1]],
        [[0,0,0,5],[0,1,0,2]],
        [[0,0,0,5],[0,0,1,3]],
        [[0,0,0,6],[1,0,0,1]],
        [[0,0,0,6],[0,0,1,3]],
        [[0,0,0,7],[1,1,0,0]],
        [[0,0,0,7],[0,1,0,2]],
    ])
    
    site_loc = site_to_lattice(spin_loc, key)
    
    # impose periodic boundary conditions
    
    site_loc_per = [
        impose_pbc(loc, side)
        for loc in site_loc
    ]
    
    # Get values associates with each site
    
    sites = [
        site_lattice[tuple(loc)]
        for loc in site_loc_per
    ]
    
    # Get initial energy for the sites
    
    e0 = sum([
        pyrochlore_site_energy(out, J)
        for out in sites
    ])
        
    # Effect of changing the spin on the sites
    # Convention: positive spin points out of first site and into second site
    
    site_lattice_out = site_lattice.copy()
    convention = [-1, 1]
    
    if spin == 1:
        sites_changed = [item[0] + item[1] for item in zip(sites, convention)]
    elif spin == -1:
        sites_changed = [item[0] - item[1] for item in zip(sites, convention)]
    else:
        raise ValueError('Lattice contains invalid spin value')
        
    # what the site lattice will become if the spin is changed
        
    for item in zip(site_loc_per, sites_changed):
        site_lattice_out[tuple(item[0])] = item[1]
        
    # energy after the change:
        
    e1 = sum([
        pyrochlore_site_energy(out, J) for out in sites_changed
    ])
        
    energy_change = e1 - e0
    
    yield energy_change
    yield site_lattice_out

def site_to_lattice(index_loc, key):
    """Finds the index locations of spins in the lattice associated with a site with given index location.
    :index_loc: np.ndarray(int), index location of site of interest in lattice of sites
    :key: np.ndarray(int), mapping from spin lattice to site lattice as described previously
    :return: list(np.ndarray(int)), list of lattice index locations associated with that site
    >>>> site_to_lattice([1, 1, 0], np.array([[[0, 0, 0], [0, 0, 1], [0, -1, 1], [-1, 0, 0]]]))
    [np.array([1, 1, 0]), np.array([1, 1, 1]), np.array([1, 0, 1]), np.array([0, 1, 0])]"""
    
    import numpy as np
    
    n = key.shape[1] # number of spins associated with a site
    D = index_loc.shape[0] - 1 # dimensionality
    
    lattice_loc = index_loc[0:D] # lattice location of unit cell
    
    basis_loc = index_loc[D] # basis label (proxy for location within unit cell)
    
    spin_loc = [
        np.append(lattice_loc, 0) + key[(basis_loc, i, ...)] 
        for i in range(0, n)
    ]
    
    return spin_loc

def impose_pbc(loc, side):
    """Ensures that spatial locations are compliant with periodic boundary conditions.
    :loc: np.ndarray(int), a spatial location (each dimension is measured in unit cell side length)
    :side: the side length of a lattice in number of unit cells
    :return: np.ndarray(int), spatial location which is compliant with periodic boundary conditions
    >>>> impose_pbc(np.array([3, 2, -1, 1]), 3)
    np.array([0, 2, 2, 1])"""
    
    import numpy as np
    
    D = loc.shape[0] - 1 # dimensionality
    
    for i in range(0, D):
        
        if loc[i] == side: # spins past the end of the lattice considered to be at the beginning
            loc[i] = 0
            
        elif loc[i] == -1: # spins before the beginning of the lattice considered to be at the end
            loc[i] = side - 1
        
        else:
            pass
        
    return loc

def pyrochlore_site_energy(out, J):
    """Return the exchange energy of a site in the pyrochlore lattice.
    :out: int, number of spins pointing out of the site (tetrahedron).
    :J: exchange energy constant
    :return: float, site energy
    >>>> pyrochlore_site_energy(2, 1)
    0"""
    
    if out == 4 or out == 0:
        E = 4 * J
    elif out == 3 or out == 1:
        E = J
    elif out == 2:
        E = 0
    else:
        raise ValueError('The site is not valid for the pyrochlore lattice')
        
    return E