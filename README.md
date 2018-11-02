# amaunet
This is an implementation of a object-oriented framework for symbolically computing a high-order perturbative expansion of the equation of state for interacting fermions in any dimension. It forms the analytic basis for a series of publications and additional upcoming manuscripts, including:

* Andrew C. Loheac and Joaquin E. Drut, Third-order perturbative lattice and complex Langevin analyses of the finite-temperature equation of state of non-relativistic fermions in one dimension. Physical Review D **95**, 033602 (2017).

* Andrew C. Loheac, Jens Braun, and Joaquin E. Drut, Polarized fermions in one dimension: density and polarization from complex Langevin calculations, perturbation theory, and the virial expansion. Physical Review D **98**, 054507 (2018).

The first article contains detailed information on the formalism used behind this code. This software does not fully evaluate the equation of state, but rather computes a symbolic representation of the contracted path-integral representation of the pressure equation of state, which results in an analytic expression in terms of fully-connected Feynman diagrams. The resulting expression can then be numerically evaluated, from which the density and polarization equations of state can be derived -- this is done in a separate code. This package is designed in a modular way such that similar expansions can be computed easily with minimal modifications to the code.

_Please note that the official repository is kept on another server, so this portfolio may not contain the latest updates._
