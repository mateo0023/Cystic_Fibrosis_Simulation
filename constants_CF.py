# General Naming: where_what_ofWhat

A_V_H2O = 1.8e-6  # Molar volume of H2O. ------ m^3 mol^-1 (meters cube / moles)
TEMP = 310.15  # Temperature is assumed to be constant. ------ K (degrees kelvin)

# Permeability Constants
P_PERM_H2O = 3.1e-5  # Permeability of water through the paracellular. ------ m/s
A_PERM_H20 = 2.4e-4  # Permeability of H2O through the apical. ------ m/s
B_PERM_H20 = 2.4e-5  # Permeability of H2O through the basolateral. ------ m/s
P_PERM_NA = 5e-8  # Permeability of sodium ions (Na+) through the paracellular. ------ m/s
P_PERM_K = 7.2e-10  # Permeability of potassium ions (K+) through the P. ------ m/s
P_PERM_CL = 1.2e-8  # Permeability chloride ions (Cl-) through the paracellular. ------ m/s

# Osmomolarity Constants
B_OSM = 279.1  # Osmomolarity on the basolateral membrane. ------ mM
PHI = 0.93
GAMMA = 0.76
B_CONS_NA = 106.4  # Na+ concentration in the basolateral compartment. ------ mM (milli-moles / litter)
B_CONS_K = 4.0  # Potassium concentration in the basolateral compartment. ------ mM (milli-moles / litter)
B_CONS_CL = 91.2  # Chloride-ion concentration in the basolateral compartment. ------ mM (milli-moles / litter)
A_CONS_OI = 2.7  # Other ions concentration in the apical compartment. ------ mM (milli-moles / litter)
C_CONS_OI = 30.4  # Concentration of other ions in the cellular compartment. ------ mM (milli-moles / litter)
B_ACT_NA = B_CONS_NA * GAMMA  # Activity of the ion, see above.
B_ACT_K = B_CONS_K * GAMMA
B_ACT_CL = B_CONS_CL * GAMMA
A_ACT_OI = A_CONS_OI * GAMMA
C_ACT_OI = C_CONS_OI * GAMMA

# Flow
FARADAY = 96485  # Faraday's Constant. ------ C / mol
R = 8.31447  # Gas constant. ------ J / (K mol)
F_RT = FARADAY / (R * TEMP)  # F/(R*T), To save computations.
J_Pump_max = 4.4e-6  # Max flow of the Na-K-ATPase pump. ------ milli-mol m^-2 s^-1
K_Na_In_pump = 0.99  # Look Ref 44 for more info. ------ mM
K_K_in_pump = 9.1  # Look Ref 44 for more info.  -   mM
K_K_ext_pump = 0.11  # Look Ref 44 for more info. ------ mM
K_Na_ext_pump = 24.3  # Look Ref 44 for more info. ------ mM

# REF 22
# COT_K_1 = 1000  # Cotransporter (Cot) off binding constant Na+. ------ mol^-1
# COT_K_2 = 0.7428  # Cot off binding constant Cl. ------ mol^-1
# COT_K_03 = 892.3  # Luminal off binding constant K. ------ mol^-1
# COT_K_I3 = 458.8  # Cytosolic off binding constant K. ------ mol^-1
# COT_K_4 = 46.71  # Cot of binding rate constant for Cl_2. ------ mol^-1
# COT_K_FF = 1824  # Cot translocation rate constant. ------ s^-1
# COT_K_BF = 2724  # Cot translocation rate constant. ------ s^-1
# COT_K_FE = 18781  # Cot translocation rate constant. ------ s^-1
# COT_K_BE = 13259  # ------ s^-1
# COT_K_ON = 1e8  # ------ mol^-1 s^-1

# NKCC2 Cotransporter - from Benjamin-Jonson (1997)
COT_D = 0.4e-6  # Density of the NKCC Cotransporter - mol / m^2
COT_K_Cl = 2.42  # mM
COT_K_Na = 22.38  # mM
COT_K_K = 234.74  # mM
COT_K_f_empty = 37767  # s^-1
COT_K_f_full = 1406  # s^-1
COT_K_b_empty = 13196  # s^-1
COT_K_b_full = 4025  # s^-1
COT_Z = (COT_K_Cl * COT_K_K * COT_K_Na * COT_K_b_empty,
         COT_K_Cl**2 * COT_K_K, COT_K_f_empty,
         COT_K_Cl * COT_K_Na * COT_K_b_empty,
         COT_K_Cl * COT_K_K * COT_K_f_empty,
         COT_K_Na * COT_K_b_empty,
         COT_K_Cl * COT_K_b_empty,
         COT_K_b_empty + COT_K_b_full,
         COT_K_f_empty + COT_K_f_full,
         COT_K_f_full / COT_K_Na,
         COT_K_f_full / COT_K_Cl,
         COT_K_b_full / (COT_K_Cl * COT_K_Na),
         COT_K_f_full / (COT_K_Cl * COT_K_K),
         COT_K_b_full / (COT_K_Cl**2 * COT_K_K),
         COT_K_b_full / (COT_K_Cl * COT_K_K * COT_K_Na),
         (COT_K_b_full + COT_K_f_full) / (COT_K_Cl**2 * COT_K_K * COT_K_Na),
         COT_K_Cl**2 * COT_K_K * COT_K_Na * (COT_K_b_empty + COT_K_f_empty))

# Permeability EQUATIONS
CACC_PERM_MAX = 6.9e-9  # The estimated MAX permeability of the CaCC Channel. ------ m / s
CACC_ATP_AT_HALF_PERM = 3.6e-5  # [ATP] when half receptors occupied (half perm). ------ muM

BK_PERM_MAX = 4.4e-9  # Similar concept, with the same units.
BK_ATP_AT_HALF_PERM = 1.5e-3  # Similar concept, with the same units.

ENAC_PERM_MAX = 2.3e-8  # Similar concept, with the same units.
ENAC_ATP_AT_HALF_PERM = 5.1e-2  # Similar concept, with the same units.

CAKC_PERM_MAX = 2.3e-7  # Similar concept, with the same units.
CAKC_ATP_AT_HALF_PERM = 3.3e-3  # Similar concept, with the same units.

# Voltage (Capacitance)
A_CAPACI = 3.23  # Capacitance of the apical membrane. ------ muF m^-2 (micro-F / meter squared)
B_CAPACI = 33.4  # Capacitance of the basolateral membrane. ------ muF m^-2 (micro-F / meter squared)

# Current uses only Faraday's constant, is shown at FLOW

# ATP, ADP, AMP, etc. NUCLEOTIDE REGULATION
J_atp = 0.0011  # Nucleotide release. ------ nmol min^-1 L^-1 (nano-moles, minutes, liters)
J_adp = 0.0131  # Nucleotide release. ------ nmol min^-1 L^-1
J_amp = 0.0125  # Nucleotide release. ------ nmol min^-1 L^-1
V_1_max = 6.3  # ------ nmol min^-1 L^-1
K_1_m = 16.3  # ------ muM (micro-moles)
V_2_max = 15.5  # ------ nmol min^-1 L^-1
K_2_m = 114.9  # ------ muM
V_3_max = 20  # ------ nmol min^-1 L^-1
K_3_m = 418  # ------ muM
V_4_max = 0.5  # ------ nmol min^-1 L^-1
K_4_m = 2.8  # ------ muM
V_5_max = 10.7  # ------ nmol min^-1 L^-1
K_5_m = 83.9  # ------ muM
V_6_max = 1.7  # ------ nmol min^-1 L^-1
K_6_m = 13  # ------ muM
V_7_max = 6.2  # ------ nmol min^-1 L^-1
K_7_m = 27.2  # ------ muM
V_8_max = 11.9  # ------ nmol min^-1 L^-1
K_8_m = 694.9  # ------ muM
V_9_max = 0.3  # ------ nmol min^-1 L^-1
K_9_m = 17  # ------ muM
V_10_max = 1.2  # ------ nmol min^-1 L^-1
K_10_m = 38.2  # ------ muM
V_F_max = 2.2  # ------ nmol min^-1 L^-1
K_F_atp = 30.4  # ------ muM
K_F_amp = 24.7  # ------ muM
V_B_max = 2.2  # ------ nmol min^-1 L^-1
K_B_adp = 61.8  # ------ muM
V_U1_max = 0.2  # ------ nmol min^-1 L^-1
K_U1_m = 1.2  # ------ muM
V_U2_max = 0.2  # ------ nmol min^-1 L^-1
K_U2_m = 1.2  # ------ muM
K_IN_atp = 28.4  # ------ muM
K_IN_adp = 20.4  # ------ muM
