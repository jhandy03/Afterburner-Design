import cantera as ct

ct.suppress_thermo_warnings()
gas = ct.Solution('JetSurf_V1_1.yaml')

gas.TP = 300,500000
gas.set_equivalence_ratio(phi=0.3872,fuel='NC12H26',oxidizer='O2:1,N2:3.76')
gas.equilibrate('HP')
print(gas.T)
gammapre = gas.cp_mass/gas.cv_mass
print(gas.density)


i_fuel = gas.species_index('NC12H26')
i_O2 = gas.species_index('O2')

gas.TP = 1045,330000
gas.set_equivalence_ratio(phi=1.0, fuel='NC12H26', oxidizer='O2:0.1267,N2:0.768364,CO2:0.05037,H2O:0.054566')
gas.equilibrate('HP')
print(gas.T)
gamma = gas.cp_mass/gas.cv_mass
print(gamma)
