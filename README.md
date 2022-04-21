# mphys_project

Some extracts from the code I wrote for my MPhys project.

## What the project is about

An exotic class of magnetic materials known as spin ice (so called due to structural and thermodynamic similarities to water ice).

All magnetic materials derive their magnetism from the magnetic moments of one or more of their microscopic constituents (which could include atoms, ions or molecules). These are referred to loosely as 'spins'. They have a single directional axis (like a bar magnet with a north pole and a south pole).

Spin ice derives its magnetisation from lanthanide ions which are arranged in a pyrochlore lattice (a repeating geometric structure made of corner sharing tetrahedra where the spins sit on the corners). The lanthanide ions have a magnetic moment, so they are hereafter referred to as the spins. Other (non-magnetic and therefore irrelevant) species exist in the structure.

Spin ice is unusual because below about 10 Kelvin (-268 degrees Celsius) the crystal environment (the combined effects of the arrangement of other nearby ions in the crystal) constrains the spins to point along certain axes. Each tetrahedron takes on its lowest energy state when two corner spins point into the tetrahedron and two point out of it. The number of ways in which this arrangement can be achieved increases exponentially with the lattice size. This means that a sample of spin ice in its ground state (lowest energy state) at low temperature may have a variety of magnetic moments, or none at all. This is in contrast with ordinary magnetic materials such as iron, where there is only one possible magnetisation (albeit in any chosen direction) at relatively low temperatures.

When the temperature is above absolue zero (-278 degrees Celsius), a material may exist not just in its ground state but also in excited states (which have energy higher than the ground state). In spin ice, these excitations are produced when a spin flips direction, producing a pair of neighbouring tetrahedra which have more or less than two spins pointing in and out. It turns out that these 'defects' behave as magnetic monopoles (isolated north and south poles), a rare and surprising phenomenon in physics (recall that ordinary bar magnets always have both a north and a south pole).

The goal of the project was to simulate the structure and basic properties of spin ice, and get some results on the behaviour of the monopoles which could be compared with some experimental work.

I encourage the reader to take a look at the mphys project report to learn more about spin ice and the nature of the simulation. There is also a paper which describes the Metropolis algorithm, a type of Monte Carlo simulation, which is used in the project. Some explanations of the methodology are also provided in the documentation.

## Contents

- lattice module: a virtual representation of a spin ice lattice and some associated methods
- metropolis algorithm module: fundamentally contains a 'spin ice sample' object which can be virtually cooled and measured.
- exchange energy module: contains functions which calculate the exchange energy of a spin ice lattice, given its state (arrangement of spins)
- magnetic coulomb energy module (no content yet)
- magnetic properties module (no content yet)
- visualisations module: currently contains a simple animation of a pyrochlore lattice. See the jupyter notebook for some example usage.

## Some design notes

I tried to make the layers of the project loosely coupled. This means that the lattice and metropolis algorithm modules can in principle be used for any magnetic lattice with Ising spins (spins which can point along only one fixed axis), or indeed for other types of Ising systems (such as beta-brass). The included exchange energy functions are specific to spin ice, but in the project I used other energy functions as well. The visualisations module contains some general functions for making 3D animations and some specific to visualising the pyrochlore lattice.

All the code in the project is based on notebooks which were used to produce the results shown in the report, but has since been refactored and polished.