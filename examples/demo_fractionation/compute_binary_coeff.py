import numpy as np 

def mass_atom_kg(m_atom):
    m_atom_kg = m_atom * 1.66053906660e-27  # Convert amu to kg
    return m_atom_kg

def mass_molecule_2atoms_kg(m_atom1, m_atom2):
    m_molecule = m_atom1 + m_atom2 # molecule mass in amu
    m_mol_kg = m_molecule * 1.66053906660e-27  # Convert to kg
    return m_mol_kg

def factor_Genda_Ikoma_2008_approx(m1, m2, m3, m4):
    factor = np.sqrt( ((m1 + m2)/(m1 * m2)) / ((m3 + m4) / (m3 * m4)) ) # Equation C3 from Genda and Ikoma 2008
    return factor

print("---------------------------------")
print ('Test to retrieve the mass of H2, HD and D2')

# Test to retrieve the Appendix results from Genda and Ikoma 2008
m_h2 = mass_molecule_2atoms_kg(1.00784, 1.00784)  # H2
m_hd = mass_molecule_2atoms_kg(1.00784, 2.014102)  # HD
m_d2 = mass_molecule_2atoms_kg(2.014102, 2.014102)  # D2
print(f"Mass of H2: {m_h2:.5e} kg")
print(f"Mass of HD: {m_hd:.5e} kg")
print(f"Mass of D2: {m_d2:.5e} kg")
factor = factor_Genda_Ikoma_2008_approx(m_h2, m_hd, m_h2, m_d2)
print(f"Factor between the 2 binary diffusion coefficients for H2,HD and H2,D2: {factor:.5e}")

print("---------------------------------")
print ('Results factor for H,N and H,S')
print("---------------------------------")

# Single atoms
m_H = mass_atom_kg(1.00784)  # H    
m_N = mass_atom_kg(14.00674)  # N
m_S = mass_atom_kg(32.065)  # S
m_O = mass_atom_kg(15.999)  # O
m_C = mass_atom_kg(12.011)  # C
print(f"Mass of H: {m_H:.5e} kg")
print(f"Mass of N: {m_N:.5e} kg")       
print(f"Mass of S: {m_S:.5e} kg")
print(f"Mass of O: {m_O:.5e} kg")
print(f"Mass of C: {m_C:.5e} kg")

factor_HN_HH = factor_Genda_Ikoma_2008_approx(m_H, m_N, m_H, m_H)
factor_HS_HH = factor_Genda_Ikoma_2008_approx(m_H, m_S, m_H, m_H)
print(f"Factor between the 2 binary diffusion coefficients for H,N and H,H: {factor_HN_HH:.5e}")   
print(f"Factor between the 2 binary diffusion coefficients for H,S and H,H: {factor_HS_HH:.5e}")

print("---------------------------------")

factor_HN_HO = factor_Genda_Ikoma_2008_approx(m_H, m_N, m_H, m_O)
factor_HS_HO = factor_Genda_Ikoma_2008_approx(m_H, m_S, m_H, m_O)
print(f"Factor between the 2 binary diffusion coefficients for H,N and H,O: {factor_HN_HO:.5e}")
print(f"Factor between the 2 binary diffusion coefficients for H,S and H,O: {factor_HS_HO:.5e}")
print("---------------------------------")

print ('Compute factor for He mixture for isofate')
print("---------------------------------")

m_He = mass_atom_kg(4.002602)  # He
print(f"Mass of He: {m_He:.5e} kg")
factor_HeN_HeO = factor_Genda_Ikoma_2008_approx(m_He, m_N, m_He, m_O)
factor_HeS_HeO = factor_Genda_Ikoma_2008_approx(m_He, m_S, m_He, m_O)

print(f"Factor between the 2 binary diffusion coefficients for He,N and He,O: {factor_HeN_HeO:.5e}")
print(f"Factor between the 2 binary diffusion coefficients for He,S and He,O: {factor_HeS_HeO:.5e}")
print("---------------------------------")


#factor_HC_HO = factor_Genda_Ikoma_2008_approx(m_H, m_C, m_H, m_O)
#factor_HC_HH = factor_Genda_Ikoma_2008_approx(m_H, m_C, m_H, m_H)
#print(f"Factor between the 2 binary diffusion coefficients for H and C with H and O: {factor_HC_HO:.5e}")
#print(f"Factor between the 2 binary diffusion coefficients for H and C with H and H: {factor_HC_HH:.5e}")

