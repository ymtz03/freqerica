from pyscf import gto
import freqerica

mol = gto.Mole()
mol.basis = "sto-6g"
mol.verbose = 0
mol.spin = 0
mol.atom = [["Li", (0.0, 0.0, 0.0)],["H", (0.0, 0.0, 1.5)]]
mol.symmetry = True
mol.build()

result = freqerica.kernel(mol, norb=5, nelec=2)
