from openfermion.hamiltonians import MolecularData
from openfermionpsi4 import run_psi4
from openfermion.utils import uccsd_generator,uccsd_convert_amplitude_format
from forestopenfermion import exponentiate,qubitop_to_pyquilpauli
from openfermion.transforms import get_fermion_operator, jordan_wigner
import numpy as np
from pyquil.quil import Program
import pyquil.api as api
from pyquil.gates import *
from forestopenfermion import exponentiate
from grove.pyvqe.vqe import VQE
from scipy.optimize import minimize

def ansatz(amps,n_orbitals=14,n_electrons=10):
    """ Return pyquil program UCC ansatz"""
    amp_singles = amps[:n_orbitals**2]
    amp_doubles = amps[n_orbitals**2:]
    x = amp_singles.reshape(n_orbitals,n_orbitals)
    y = amp_doubles.reshape(n_orbitals,n_orbitals,n_orbitals,n_orbitals)
    p = Program()
    for i in range(n_electrons):
        p += X(i)
    p += exponentiate(jordan_wigner(uccsd_generator(x,y)) / (-1j))
    return p

def expectation(p1,p2,p3,multi=1):
    """ 
    Return UCC expectation value for a specified geometry 
    * First runs a psi4 ccsd calculation to get single and double
      amplitudes to use as ansatz for UCC
    * Generates a Hamiltonian for the specified geometry
    * Obtains expectation value using VQE 
    """
    geometry = [['O',p1], ['H',p2], ['H',p3]]
    molecule = MolecularData(geometry, 
                            basis='sto-3g',
                            multiplicity=multi,
                            description=str(round(rad, 2)) 
                            + "_" + str(round(ang,2))
                            )
    # Run Psi4.
    molecule = run_psi4(molecule,
                    run_ccsd=1,
                    run_fci=1)
    # Print out some results of calculation.
    print('\nRAD: {}, ANG: {}\n'.format(rad,ang))
    print('FCI energy: {} Hartree.'.format(molecule.fci_energy))

    singles_initial = molecule.ccsd_single_amps.flatten()
    doubles_initial = molecule.ccsd_double_amps.flatten()
    amps = np.concatenate((singles_initial,doubles_initial),axis=0)

    print("Compiling the Hamiltonian...")
    hamiltonian=jordan_wigner(get_fermion_operator(molecule.get_molecular_hamiltonian()))
    hamiltonian.compress()
    hamiltonian = qubitop_to_pyquilpauli(hamiltonian)
    print("Hamiltonian complete")

    vqe = VQE(minimizer=minimize, 
              minimizer_kwargs={'method': 'nelder-mead',
                                'options': {'fatol':1.5e-3}}
              )
    result = vqe.expectation(ansatz(amps), hamiltonian, None, qvm)
    print("VQE Expectation Value: {} Hartree".format(result))

    return result

qvm = api.QVMConnection()
# Set grid parameters
N = 12
r_min = 0.6
r_max = 1.5
a_min = 30
a_max = 180

# Get expectation values for a grid of parameters
for rad in np.linspace(r_min,r_max,N):
    for ang in np.linspace(a_min,a_max,N):
        p1 = [0,0,0]
        p2 = [rad,0,0]
        p3 = [rad*np.cos(ang*np.pi/180),rad*np.sin(ang*np.pi/180),0]

        exp = expectation(p1,p2,p3) 
        with open('resources/contour_results.txt','a') as f:
            f.write(str(exp)+"\n")
        print("Result saved.") 
