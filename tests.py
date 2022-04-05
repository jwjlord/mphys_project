########## Tests file ###########

from .lattice import *
from .metropolis import *
from .coulomb import *
from .exchange import *

# ----- lattice ------

def test_create_lattice():
    test_lattice = Lattice(basis, side)
    assert test_lattice.basis == basis
    assert test_lattice.side == side
    assert test_lattice.spins == 0
    
# ----- metropolis ------

def test_run_metropolis():
    assert MetropolisRun.cool() is None
    assert _determine_flip() == 0
    assert mean_magnetisation == None
    assert mean_monopoles == None
    
# ------ exchange ------

def test_exchange():
    assert pyrochlore_exchange_energy(lattice, site_lattice) == 0
    assert pyrochlore_exchange_energy_change(lattice, site_lattice, spin_index) == 0
    
# ------ coulomb ------

def test_coulomb():
    assert coulomb_energy(site_lattice) == 0
    assert coulomb_energy_change(site_lattice, lattice, spin_index) == 0
    