import cantera as ct
from pathlib import Path

ct.suppress_thermo_warnings()

# Resolve paths relative to this script so loading works from any working directory.
MECH_DIR = Path(__file__).resolve().parent
MECH_FILE = MECH_DIR / 'JetSurf_V1_1.yaml'
ct.add_directory(str(MECH_DIR))
gas = ct.Solution(str(MECH_FILE))

gas.TP = 300,500000
gas.set_equivalence_ratio(phi=0.3872,fuel='NC12H26',oxidizer='O2:1,N2:3.76')
gas.equilibrate('HP')
gammapre = gas.cp_mass/gas.cv_mass

i_fuel = gas.species_index('NC12H26')
i_O2 = gas.species_index('O2')

gas.TP = 1045,330000
gas.set_equivalence_ratio(phi=1.0, fuel='NC12H26', oxidizer='O2:0.1267,N2:0.768364,CO2:0.05037,H2O:0.054566')
gas.equilibrate('HP')

