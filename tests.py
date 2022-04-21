########## Tests file ###########

# TODO: some features of metropolis_run are not tested
# visualisations module has no unit tests

from .lattice import *
from .metropolis import *
from .coulomb import *
from .exchange import *
import numpy as np

# ----- lattice ------

def test_create_lattice():
    basis =  np.array([[0, 0], [0.5, 0.5]])
    side = 1
    test_lattice = Lattice(basis, side, seed=42)
    assert np.all(test_lattice.basis == basis)
    assert test_lattice.side == side
    assert np.all(test_lattice.spin_lattice == np.array([[[1, -1]]]))
    
    # read a line of a lattice file
    
    with open('test_read_file.txt') as f:
        assert np.all(read_lattice(f.read()) == np.array([[[[-1, 1,-1,-1,-1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1]]]]))
    
# ----- metropolis ------

def test_run_metropolis():
    ice_basis = np.array([[0,0,0],[0,0,0],[0,0.5,0.5],[0,0.5,0.5],[0.5,0,0.5],[0.5,0,0.5],[0.5,0.5,0],[0.5,0.5,0],[0.25,0.25,0.25],[0.25,0.25,0.25],[0.25, 0.75, 0.75],[0.25, 0.75, 0.75],[0.75, 0.25, 0.75],[0.75, 0.25, 0.75],[0.75, 0.75, 0.25],[0.75, 0.75, 0.25]])
    side = 1
    test_lattice = Lattice(ice_basis, side, seed=42)
    test_run = MetropolisRun(test_lattice, energy_function=pyrochlore_exchange_energy, energy_change_function=pyrochlore_exchange_energy_change, J=1)
    
    # test __init__():
    
    assert test_run.n == 16 and test_run.D == 3
    
    assert test_run.initial_energy == 10
    assert np.all(test_run.initial_site_lattice == np.array([[[[3, 3, 3, 2, 1, 1, 3, 0]]]]))
    
    assert test_run.params['J'] == 1
    assert np.all(test_run.params['site_lattice']) == np.all(np.array([[[[3, 3, 3, 2, 1, 1, 3, 0]]]]))
    
    test_run.cool(np.linspace(5, 0.01, 100), 1)
    
    assert test_run.temp == 0.01
    assert test_run.energy == 0
    assert np.all(test_run.params['site_lattice']) == np.all(np.array([[[[2, 2, 2, 2, 2, 2, 2, 2]]]]))
    assert len(test_run.energies) == 100
    
    test_run.reset()
    
    assert test_run.energy == 10
    
    test_run.cool(np.linspace(5, 0.01, 100), 1, lattice_filename='test_cool')
    
    assert test_run.temp == 0.01
    assert test_run.energy == 0
    assert np.all(test_run.params['site_lattice']) == np.all(np.array([[[[2, 2, 2, 2, 2, 2, 2, 2]]]]))
    assert len(test_run.energies) == 100
    
    
#    assert _determine_flip() == 0
    
# ------ exchange ------

def test_exchange():
    
    # pyrochlore_exchange_energy
    
    ice_basis = np.array([[0,0,0],[0,0,0],[0,0.5,0.5],[0,0.5,0.5],[0.5,0,0.5],[0.5,0,0.5],[0.5,0.5,0],[0.5,0.5,0],[0.25,0.25,0.25],[0.25,0.25,0.25],[0.25, 0.75, 0.75],[0.25, 0.75, 0.75],[0.75, 0.25, 0.75],[0.75, 0.25, 0.75],[0.75, 0.75, 0.25],[0.75, 0.75, 0.25]])
    side = 1
    test_lattice = Lattice(ice_basis, side, seed=42)
    test_lattice_spins = test_lattice.spin_lattice
    J = 1
    test_site_lattice = np.array([[[[3, 3, 3, 2, 1, 1, 3, 0]]]])
    lattice_energy = pyrochlore_exchange_energy(test_lattice_spins, {'J': J})
    assert next(lattice_energy) == 10
    assert np.all(next(lattice_energy) == test_site_lattice)
    
    # pyrochlore_exchange_energy_change
    
    energy_change = pyrochlore_exchange_energy_change(test_lattice_spins, np.array([0, 0, 0, 1]), {'J': J, 'site_lattice': test_site_lattice})
    assert next(energy_change) == 2
    assert np.all(next(energy_change) == np.array([[[[4, 3, 3, 2, 1, 1, 2, 0]]]]))
    
    # site_to_lattice
    
    lattice_sites = site_to_lattice(
        np.array([1, 1, 0]), 
        np.array([[[0, 0, 0], [0, 0, 1], [0, -1, 1], [-1, 0, 0]]])
    )
    
    lattice_sites_output = [
        np.array([1, 1, 0]), 
        np.array([1, 1, 1]), 
        np.array([1, 0, 1]), 
        np.array([0, 1, 0])]
    
    for element in [*zip(lattice_sites, lattice_sites_output)]:
        assert np.all(element[0] == element[1])
        
    # impose_pbc
    
    assert np.all(
        impose_pbc(np.array([3, 2, -1, 1]), 3) == np.array([0, 2, 2, 1])
    )
        
    # pyrochlore_site_energy
        
    assert pyrochlore_site_energy(2, 1) == 0
    
# ------ coulomb ------

def test_coulomb():
    assert coulomb_energy(site_lattice) == 0
    assert coulomb_energy_change(site_lattice, lattice, spin_index) == 0
    