import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

ct.suppress_thermo_warnings()
gas = ct.Solution('JetSurf_V1_1.yaml')
gas.TP = 300,500000
gas.set_equivalence_ratio(phi=0.3872,fuel='NC12H26',oxidizer='O2:1,N2:3.76')
gas.equilibrate('HP')

i_fuel = gas.species_index('NC12H26')
i_O2 = gas.species_index('O2')

print(f'Temp = {gas.T} K, and the mole fraction of O2 = {gas.X[i_O2]}')
gas.TP = 1045,330000

phi_values = np.linspace(0.7,1.1,100)

# Variables to store the maximums
max_temp = -1.0
max_phi = -1.0
temp_values = np.ones(len(phi_values))

for i, phi in enumerate(phi_values):
    gas.set_equivalence_ratio(phi=phi, fuel='NC12H26', oxidizer='O2:0.1267,N2:0.768364,CO2:0.05037,H2O:0.054566')
    gas.equilibrate('HP')
    
    current_temp = gas.T 
    temp_values[i] = gas.T
    # print(f'Phi = {phi:.4f}, Temp = {current_temp:.2f} K')
    
    # Check if this is the highest temperature seen so far
    if current_temp > max_temp:
        max_temp = current_temp
        max_phi = phi

print(f'\nMaximum Temperature found: {max_temp:.2f} K at Phi = {max_phi:.4f}')

# # Plot
# plt.figure(figsize=(10, 6))
# plt.plot(phi_values, temp_values, 'b-', linewidth=2)
# plt.xlabel('Equivalence Ratio (φ)')
# plt.ylabel('Temperature (K)')
# plt.title('Temperature vs Equivalence Ratio')
# plt.grid(True, alpha=0.3)
# plt.plot(max_phi, max_temp, 'ro', markersize=8, label=f'Max: {max_temp:.2f} K at φ={max_phi:.4f}')
# plt.legend()
# plt.tight_layout()
# plt.show()
