import pandas as pd
import numpy as np
from itertools import product as prod
from scipy.optimize import least_squares

class Constants:
	"""
	This is a class that will be used to store the values of the constants.
	"""

	def __init__(self, cf, init_vals=None, constants=None):
		"""
		It will create all of the necessary constants for the model to run.

		General Naming: where_what_ofWhat
			* Concentration is "CONS" because it sounds closer.
		
		:arg cf			The cf is a boolean whether to default to the Cystic Fibrosis or regular constants
		:arg init_vals	The initial values of the system (necessary to calculate the parameters).
						The calculation will be done if constants is None and init_vals is set.
		:arg constants	The constants should be a dictionary with the precisse name and value of the
						constants to which one wants to override the default value.

		"""
		self.isCF = cf

		self.data = {}

		if self.isCF:
			# This should have all of the parameters that need to be calculated
			self.parametersMeanVariance = {'P_PERM_NA': (5e-8, 0.08), 'P_PERM_K': (7.2e-10, .15),
											'P_PERM_CL': (1.2e-8, 0.07), 'J_Pump_max': (4.4e-6, 0.09),
											'CACC_PERM_MAX': (6.9e-9, 0.03), 'CACC_ATP_AT_HALF_PERM': (3.6e-5, 0.28),
											'BK_PERM_MAX': (4.4e-9, 0.16), 'BK_ATP_AT_HALF_PERM': (1.5e-3, 0.38),
											'ENAC_PERM_MAX': (2.3e-8, 0.08), 'ENAC_ATP_AT_HALF_PERM': (5.1e-2, 0.07),
											'CAKC_PERM_MAX': (2.3e-7, 0.14), 'CAKC_ATP_AT_HALF_PERM': (3.3e-3, 0.15)}
			self.variable_params = [5e-08, 7.2e-10, 1.2e-08, 4.4e-06, 6.9e-09,
									3.6e-05, 4.4e-09, 0.0015, 2.3e-08, 0.051, 2.3e-07, 0.0033]
			self.indexDict = {"P_PERM_NA": 0, "P_PERM_K": 1,
								"P_PERM_CL": 2,"J_Pump_max": 3,
								"CACC_PERM_MAX": 4, "CACC_ATP_AT_HALF_PERM": 5,
								"BK_PERM_MAX": 6, "BK_ATP_AT_HALF_PERM": 7,
								"ENAC_PERM_MAX": 8, "ENAC_ATP_AT_HALF_PERM": 9,
								"CAKC_PERM_MAX": 10, "CAKC_ATP_AT_HALF_PERM": 11}
		else:
			# This should have all of the parameters that need to be calculated
			self.paremetersMeanVariance = {'P_PERM_NA': (1e-8, 0.06), 'P_PERM_K': (3.4e-10, 0.24),
											'P_PERM_CL': (3.9e-8, 0.05), 'J_Pump_max': (1.5e-6, 0.06),
											'CACC_PERM_MAX': (4.8e-8, 0.25), 'CACC_ATP_AT_HALF_PERM': (8.1e-5, 0.37),
											'BK_PERM_MAX': (7.4e-10, 0.22), 'BK_ATP_AT_HALF_PERM': (1.3e-3, 0.43),
											'ENAC_PERM_MAX': (8.3e-9, 0.07), 'ENAC_ATP_AT_HALF_PERM': (7.2e-2, 0.23), 'ENAC_ADO_AT_HALF_PERM': (3.1e-1, 0.28),
											'CAKC_PERM_MAX': (5.7e-8, 0.08), 'CAKC_ATP_AT_HALF_PERM': (2.1e-4, 0.26),
											'CFTR_PERM_MAX': (9.4e-8, 0.29), 'CFTR_ATP_AT_HALF_PERM': (6.8e-5, 0.34), 'CFTR_ADO_AT_HALF_PERM': (1.9e-1, 0.53)}
			self.variable_params = [1e-08, 3.4e-10, 3.9e-08, 1.5e-06, 4.8e-08, 8.1e-05, 7.4e-10, 0.0013,
									8.3e-09, 0.072, 0.31, 5.7e-08, 0.00021, 9.4e-08, 6.8e-05, 0.19]
			self.indexDict = {"P_PERM_NA": 0, "P_PERM_K": 1,
							  "P_PERM_CL": 2, "J_Pump_max": 3,
							  "CACC_PERM_MAX": 4, "CACC_ATP_AT_HALF_PERM": 5,
							  "BK_PERM_MAX": 6, "BK_ATP_AT_HALF_PERM": 7,
							  "ENAC_PERM_MAX": 8, "ENAC_ATP_AT_HALF_PERM": 9, "ENAC_ADO_AT_HALF_PERM": 10,
							  "CAKC_PERM_MAX": 11, "CAKC_ATP_AT_HALF_PERM": 12,
							  "CFTR_PERM_MAX": 13, "CFTR_ATP_AT_HALF_PERM": 14, "CFTR_ADO_AT_HALF_PERM": 15}

		# Create lambda functions for the constants that are calculated based on others.
		self.calcCELL_H = lambda: (self.data['CELL_VOL'] / 2) ** (1 / 3) * 2e-6
		self.calcB_ACT_NA = lambda: self.data['B_CONS_NA'] * self.data['GAMMA']
		self.calcB_ACT_K = lambda: self.data['B_CONS_K'] * self.data['GAMMA']
		self.calcB_ACT_CL = lambda: self.data['B_CONS_CL'] * self.data['GAMMA']
		self.calcA_ACT_OI = lambda: (self.data['A_CONS_IO'] + self.data['A_CONS_IO']
										* self.data['GAMMA'])
		self.calcC_ACT_OI = lambda: (self.data['C_CONS_IO'] + self.data['C_CONS_NCL']
										* self.data['GAMMA'])
		self.calcF_RT = lambda: self.data['FARADAY'] / (self.data['R'] * self.data['TEMP'])
		self.calcCOT_Z = lambda: (self.data['COT_K_Cl'] * self.data['COT_K_K'] * self.data['COT_K_Na'] * self.data['COT_K_b_empty'],
							self.data['COT_K_Cl'] ** 2 * self.data['COT_K_K'], self.data['COT_K_f_empty'],
							self.data['COT_K_Cl'] * self.data['COT_K_Na'] * self.data['COT_K_b_empty'],
							self.data['COT_K_Cl'] * self.data['COT_K_K'] * self.data['COT_K_f_empty'],
							self.data['COT_K_Na'] * self.data['COT_K_b_empty'],
							self.data['COT_K_Cl'] * self.data['COT_K_b_empty'],
							self.data['COT_K_b_empty'] + self.data['COT_K_b_full'],
							self.data['COT_K_f_empty'] + self.data['COT_K_f_full'],
							self.data['COT_K_f_full'] / self.data['COT_K_Na'],
							self.data['COT_K_f_full'] / self.data['COT_K_Cl'],
							self.data['COT_K_b_full'] / (self.data['COT_K_Cl'] * self.data['COT_K_Na']),
							self.data['COT_K_f_full'] / (self.data['COT_K_Cl'] * self.data['COT_K_K']),
							self.data['COT_K_b_full'] / (self.data['COT_K_Cl'] ** 2 * self.data['COT_K_K']),
							self.data['COT_K_b_full'] / (self.data['COT_K_Cl'] * self.data['COT_K_K']
														* self.data['COT_K_Na']),
							(self.data['COT_K_b_full'] + self.data['COT_K_f_full']) /
							(self.data['COT_K_Cl'] ** 2 * self.data['COT_K_K']
							* self.data['COT_K_Na']),
							self.data['COT_K_Cl'] ** 2 * self.data['COT_K_K'] * self.data['COT_K_Na']
							* (self.data['COT_K_b_empty'] + self.data['COT_K_f_empty']))

		# Craete a dictionary to store the names of the variables that are calculated,
		# which constants are needed to calculate them and the lambda to calculate them.
		self.calculated = { 'CELL_H': (('CELL_VOL'), self.calcCELL_H),
							'B_ACT_NA': (('B_CONS_NA', 'GAMMA'), self.calcB_ACT_NA),
							'B_ACT_K': (('B_CONS_K', 'GAMMA'), self.calcB_ACT_K),
							'B_ACT_CL': (('B_CONS_CL', 'GAMMA'), self.calcB_ACT_CL),
							'A_ACT_OI': (('A_CONS_IO', 'A_CONS_IO', 'GAMMA'), self.calcA_ACT_OI),
							'C_ACT_OI': (('C_CONS_IO', 'C_CONS_NCL', 'GAMMA'), self.calcC_ACT_OI),
							'F_RT': (('FARADAY', 'R', 'TEMP'), self.calcF_RT)}

		# All of the constants used for the Benjamin-Johnson model.
		# In case that COT_Z needs to be re-calculated.
		self.cotransporter_benj = ('COT_D', 'COT_K_Cl', 'COT_K_Na', 'COT_K_K', 'COT_K_f_empty',
									'COT_K_f_full', 'COT_K_b_empty', 'COT_K_b_full')

		if cf:
			# Permeability Constants
			self.data['P_PERM_NA'] = 5e-8  # Permeability of sodium ions (Na+) through the paracellular. ------ m/s
			self.data['P_PERM_K'] = 7.2e-10	 # Permeability of potassium ions (K+) through the P. ------ m/s
			self.data['P_PERM_CL'] = 1.2e-8	 # Permeability chloride ions (Cl-) through the paracellular. ------ m/s

			# Osmomolarity Constants
			self.data['A_CONS_NCL'] = 48.7	# Apical concentration on non-Chloride anions ------ mM
			self.data[
				'A_CONS_IO'] = 2.7	# Concentration of impermeable osmolytes in the apical compartment. ------ mM (milli-moles / litter)
			self.data['C_CONS_NCL'] = 69.1	# Cellular concentration on non-Chloride anions ------ mM
			self.data[
				'C_CONS_IO'] = 30.4	 # Concentration of impermeable osmolytes in the cellular compartment. ------ mM (milli-moles / litter)

			# Flow
			self.data['J_Pump_max'] = 4.4e-6  # Max flow of the Na-K-ATPase pump. ------ milli-mol m^-2 s^-1

			# Permeability EQUATIONS
			self.data['CACC_PERM_MAX'] = 6.9e-9	 # The estimated MAX permeability of the CaCC Channel. ------ m / s
			self.data['CACC_ATP_AT_HALF_PERM'] = 3.6e-5	 # [ATP] when half receptors occupied (half perm). ------ muM

			self.data['BK_PERM_MAX'] = 4.4e-9  # Similar concept, with the same units.
			self.data['BK_ATP_AT_HALF_PERM'] = 1.5e-3  # Similar concept, with the same units.

			self.data['ENAC_PERM_MAX'] = 2.3e-8	 # Similar concept, with the same units.
			self.data['ENAC_ATP_AT_HALF_PERM'] = 5.1e-2	 # Similar concept, with the same units.

			self.data['CAKC_PERM_MAX'] = 2.3e-7	 # Similar concept, with the same units.
			self.data['CAKC_ATP_AT_HALF_PERM'] = 3.3e-3	 # Similar concept, with the same units.
		else:
			# Permeability Constants
			self.data['P_PERM_NA'] = 1e-8  # Permeability of sodium ions (Na+) through the paracellular. ------ m/s
			self.data['P_PERM_K'] = 3.4e-10	 # Permeability of potassium ions (K+) through the P. ------ m/s
			self.data['P_PERM_CL'] = 3.9e-8	 # Permeability chloride ions (Cl-) through the paracellular. ------ m/s

			# Osmomolarity Constants
			self.data['A_CONS_NCL'] = 40.8	# Apical concentration on non-Chloride anions ------ mM
			self.data['A_CONS_IO'] = 8	# Concentration of impermeable osmolytes in the apical compartment. ------ mM (milli-moles / litter)
			self.data['C_CONS_NCL'] = 77.6	# Cellular concentration on non-Chloride anions ------ mM
			self.data['C_CONS_IO'] = 29	 # Concentration of impermeable osmolytes in the cellular compartment. ------ mM (milli-moles / litter)

			# Flow
			self.data['J_Pump_max'] = 1.5e-6  # Max flow of the Na-K-ATPase pump. ------ milli-mol m^-2 s^-1

			# Permeability EQUATIONS
			self.data['CACC_PERM_MAX'] = 4.8e-9	 # The estimated MAX permeability of the CaCC Channel. ------ m / s
			self.data['CACC_ATP_AT_HALF_PERM'] = 8.1e-5	 # [ATP] when half receptors occupied (half perm). ------ muM

			self.data['BK_PERM_MAX'] = 7.4e-9  # Similar concept, with the same units.
			self.data['BK_ATP_AT_HALF_PERM'] = 1.3e-3  # Similar concept, with the same units.

			self.data['ENAC_PERM_MAX'] = 8.3e-9	 # Similar concept, with the same units.
			self.data['ENAC_ATP_AT_HALF_PERM'] = 7.2e-2	 # Similar concept, with the same units.
			self.data['ENAC_ADO_AT_HALF_PERM'] = 3.1e-1	 # Similar concept, with the same units.

			self.data['CAKC_PERM_MAX'] = 2.3e-7	 # Similar concept, with the same units.
			self.data['CAKC_ATP_AT_HALF_PERM'] = 3.3e-3	 # Similar concept, with the same units.

			self.data['CFTR_PERM_MAX'] = 9.4e-8	 # Similar concept, with the same units.
			self.data['CFTR_ATP_AT_HALF_PERM'] = 6.8e-5	 # Similar concept, with the same units.
			self.data['CFTR_ADO_AT_HALF_PERM'] = 1.e-1 # Similar concept, with the same units.

		self.data['A_V_H2O'] = 1.8e-6  # Molar volume of H2O. ------ m^3 mol^-1 (meters cube / moles)
		self.data['TEMP'] = 310.15	# Temperature is assumed to be constant. ------ K (degrees kelvin)

		self.data['CELL_VOL'] = 1450	# Cell volume of water, used for init cell height. ------ (mu m)^3 (micro-meters cubed)
		self.data['CELL_H'] = self.calcCELL_H()	 # Cell shape assumed to be around 2NxNxN. ------ m (meters)

		# Permeability Constants
		self.data['P_PERM_H2O'] = 3.1e-5  # Permeability of water through the paracellular. ------ m/s
		self.data['A_PERM_H20'] = 2.4e-4  # Permeability of H2O through the apical. ------ m/s
		self.data['B_PERM_H20'] = 2.4e-5  # Permeability of H2O through the basolateral. ------ m/s

		# Osmomolarity & Concentration Constants
		self.data['B_OSM'] = 279.1	# Osmomolarity on the basolateral membrane. ------ mM
		self.data['PHI'] = 0.93
		self.data['GAMMA'] = 0.76
		self.data['B_CONS_NA'] = 106.4  # Na+ concentration in the basolateral compartment. ------ mM (milli-moles / litter)
		self.data['B_CONS_K'] = 4.0  # Potassium concentration in the basolateral compartment. ------ mM (milli-moles / litter)
		self.data['B_CONS_CL'] = 91.2	 # Chloride-ion concentration in the basolateral compartment. ------ mM (milli-moles / litter)
		self.data['B_ACT_NA'] = self.calcB_ACT_NA()	 # Activity of the ion, see above.
		self.data['B_ACT_K'] = self.calcB_ACT_K()
		self.data['B_ACT_CL'] = self.calcB_ACT_CL()
		self.data['A_ACT_OI'] = self.calcC_ACT_OI()	 # Activity of the non-Chloride ions and the impermeable osmolytes.
		self.data['C_ACT_OI'] = self.calcC_ACT_OI()	 # Activity of the non-Chloride ions and the impermeable osmolytes.

		# Flow
		self.data['FARADAY'] = 96485  # Faraday's Constant. ------ C / mol
		self.data['R'] = 8.31447  # Gas constant. ------ J / (K mol)
		self.data['F_RT'] = self.calcF_RT()	 # F/(R*T), To save computations.
		self.data['K_Na_In_pump'] = 0.99  # Look Ref 44 for more info. ------ mM
		self.data['K_K_in_pump'] = 9.1	# Look Ref 44 for more info.  -	  mM
		self.data['K_K_ext_pump'] = 0.11  # Look Ref 44 for more info. ------ mM
		self.data['K_Na_ext_pump'] = 24.3  # Look Ref 44 for more info. ------ mM

		# REF 22
		# self.data['COT_K_1'] = 1000  # Cotransporter (Cot) off binding constant Na+. ------ mol^-1
		# self.data['COT_K_2'] = 0.7428	 # Cot off binding constant Cl. ------ mol^-1
		# self.data['COT_K_03'] = 892.3	 # Luminal off binding constant K. ------ mol^-1
		# self.data['COT_K_I3'] = 458.8	 # Cytosolic off binding constant K. ------ mol^-1
		# self.data['COT_K_4'] = 46.71	# Cot of binding rate constant for Cl_2. ------ mol^-1
		# self.data['COT_K_FF'] = 1824	# Cot translocation rate constant. ------ s^-1
		# self.data['COT_K_BF'] = 2724	# Cot translocation rate constant. ------ s^-1
		# self.data['COT_K_FE'] = 18781	 # Cot translocation rate constant. ------ s^-1
		# self.data['COT_K_BE'] = 13259	 # ------ s^-1
		# self.data['COT_K_ON'] = 1e8  # ------ mol^-1 s^-1

		# NKCC2 Cotransporter - from Benjamin-Jonson (1997)
		self.data['COT_D'] = 0.4e-6	 # Density of the NKCC Cotransporter - mol / m^2
		self.data['COT_K_Cl'] = 2.42  # mM
		self.data['COT_K_Na'] = 22.38  # mM
		self.data['COT_K_K'] = 234.74  # mM
		self.data['COT_K_f_empty'] = 37767	# s^-1
		self.data['COT_K_f_full'] = 1406  # s^-1
		self.data['COT_K_b_empty'] = 13196	# s^-1
		self.data['COT_K_b_full'] = 4025  # s^-1
		self.data['COT_Z'] = self.calcCOT_Z()

		# Voltage (Capacitance)
		self.data['A_CAPACI'] = 3.23  # Capacitance of the apical membrane. ------ muF m^-2 (micro-F / meter squared)
		self.data['B_CAPACI'] = 33.4	# Capacitance of the basolateral membrane. ------ muF m^-2 (micro-F / meter squared)

		# Current uses only Faraday's constant, is shown at FLOW

		# ATP, ADP, AMP, etc. NUCLEOTIDE REGULATION
		self.data['J_atp'] = 0.0011	 # Nucleotide release. ------ nmol min^-1 L^-1 (nano-moles, minutes, liters)
		self.data['J_adp'] = 0.0131	 # Nucleotide release. ------ nmol min^-1 L^-1
		self.data['J_amp'] = 0.0125	 # Nucleotide release. ------ nmol min^-1 L^-1
		self.data['V_1_max'] = 6.3	# ------ nmol min^-1 L^-1
		self.data['K_1_m'] = 16.3  # ------ muM (micro-moles)
		self.data['V_2_max'] = 15.5	 # ------ nmol min^-1 L^-1
		self.data['K_2_m'] = 114.9	# ------ muM
		self.data['V_3_max'] = 20  # ------ nmol min^-1 L^-1
		self.data['K_3_m'] = 418  # ------ muM
		self.data['V_4_max'] = 0.5	# ------ nmol min^-1 L^-1
		self.data['K_4_m'] = 2.8  # ------ muM
		self.data['V_5_max'] = 10.7	 # ------ nmol min^-1 L^-1
		self.data['K_5_m'] = 83.9  # ------ muM
		self.data['V_6_max'] = 1.7	# ------ nmol min^-1 L^-1
		self.data['K_6_m'] = 13	 # ------ muM
		self.data['V_7_max'] = 6.2	# ------ nmol min^-1 L^-1
		self.data['K_7_m'] = 27.2  # ------ muM
		self.data['V_8_max'] = 11.9	 # ------ nmol min^-1 L^-1
		self.data['K_8_m'] = 694.9	# ------ muM
		self.data['V_9_max'] = 0.3	# ------ nmol min^-1 L^-1
		self.data['K_9_m'] = 17	 # ------ muM
		self.data['V_10_max'] = 1.2	 # ------ nmol min^-1 L^-1
		self.data['K_10_m'] = 38.2	# ------ muM
		self.data['V_F_max'] = 2.2	# ------ nmol min^-1 L^-1
		self.data['K_F_atp'] = 30.4	 # ------ muM
		self.data['K_F_amp'] = 24.7	 # ------ muM
		self.data['V_B_max'] = 2.2	# ------ nmol min^-1 L^-1
		self.data['K_B_adp'] = 61.8	 # ------ muM
		self.data['V_U1_max'] = 0.2	 # ------ nmol min^-1 L^-1
		self.data['K_U1_m'] = 1.2  # ------ muM
		self.data['V_U2_max'] = 0.2	 # ------ nmol min^-1 L^-1
		self.data['K_U2_m'] = 1.2  # ------ muM
		self.data['K_IN_atp'] = 28.4  # ------ muM
		self.data['K_IN_adp'] = 20.4  # ------ muM

		if isinstance(constants, dict):
			for key, val in constants.items():
				if key in self.data:
					self.data[key] = val
				else:
					pass
			
			for cons, val in self.calculated.items():
				if any(k in constants for k in val[0]):
					self.data[cons] = val[1]()

			if any(k in constants for k in self.cotransporter_benj):
				self.data['COT_Z'] = self.calcCOT_Z()
		if init_vals is not None and constants is None:
			self.getSteadyVals(init_vals)

	def __getitem__(self, item):
		return self.data[item]

	def getSteadyVals(self, init_vals):
		"""
		This will calculate the parameters according to the steady state of the system.

		It sets the five equations to 0:
			\\frac{dN_Na^c}{dt} = 0
			\\frac{dN_K^c}{dt} = 0
			\\frac{dN_C^cl}{dt} = 0
			\\frac{dV_a}{dt} = 0
			\\frac{dV_b}{dt} = 0
		
		Uses the Scipy.optimize.least_squares function.
		"""

		# Eq 5.1.1
		par_p_CaCC = lambda p: p[self.indexDict['CACC_PERM_MAX']] / (1 + p[self.indexDict['CACC_ATP_AT_HALF_PERM']] / init_vals["ATP"])

		# Eq 5.1.2
		par_p_CFTR_NL = lambda p: p[self.indexDict['CFTR_PERM_MAX']] / (1 + p[self.indexDict['CFTR_ADO_AT_HALF_PERM']] / init_vals["ADO"] +
											p[self.indexDict['CFTR_ATP_AT_HALF_PERM']] / init_vals["ATP"])

		# Eq 5.2
		par_p_BK = lambda p: p[self.indexDict['BK_PERM_MAX']] / (1 + p[self.indexDict['BK_ATP_AT_HALF_PERM']] / init_vals["ATP"])

		# Eq 5.3 for NL
		par_p_ENaC_NL = lambda p: p[self.indexDict['ENAC_PERM_MAX']] / (1 + init_vals['ATP'] / p[self.indexDict['ENAC_ATP_AT_HALF_PERM']]
											+ init_vals['ADO'] / p[self.indexDict['ENAC_ADO_AT_HALF_PERM']])

		par_p_ENaC_CF = lambda p: p[self.indexDict['ENAC_PERM_MAX']] / (1 +
			init_vals['ATP'] / p[self.indexDict['ENAC_ATP_AT_HALF_PERM']])

		# Eq 5.4
		par_p_CaKC = lambda p: p[self.indexDict['CAKC_PERM_MAX']] / (1 + p[self.indexDict['CAKC_ATP_AT_HALF_PERM']] / init_vals['ATP'])


		if self.isCF:
			par_p_ENaC = par_p_ENaC_CF
			par_p_CFTR = par_p_CFTR_NL
		else:
			par_p_ENaC = par_p_ENaC_NL
			par_p_CFTR = lambda x: 0

		par_aJ_Na = lambda param: par_p_ENaC(param) * init_vals["aV"] * self.data['F_RT'] \
				* (init_vals["aNa"] -
					init_vals["cNa"] * np.exp(init_vals["aV"] * self.data['F_RT'])) \
				/ (np.exp(init_vals["aV"] * self.data['F_RT']) - 1)

		# Eq 4.2
		par_aJ_K = lambda param: par_p_BK(param) * init_vals["aV"] * self.data['F_RT'] / (
					np.exp(init_vals["aV"] * self.data['F_RT']) - 1) \
				* (init_vals["aK"] - init_vals["cK"] * np.exp(
				self.data['F_RT'] * init_vals["aV"]))

		# Eq 4.3
		par_aJ_Cl = lambda param: -(par_p_CaCC(param) + par_p_CFTR(param)) \
				* self.data['F_RT'] * init_vals["aV"] / (np.exp(self.data['F_RT'] * init_vals["aV"]) - 1) \
				* (init_vals["cCl"] - init_vals["aCl"] * np.exp(
				self.data['F_RT'] * init_vals["aV"]))

		# Eq 4.5
		par_bJ_Cl = lambda param: - param[self.indexDict['P_PERM_CL']] * self.data['F_RT'] * init_vals["bV"] \
				/ (np.exp(self.data['F_RT'] * init_vals["bV"]) - 1) \
				* (init_vals['cCl'] - self.data['B_ACT_CL'] * np.exp(
				self.data['F_RT'] * init_vals["bV"]))

		# Eq 4.6
		par_bJ_K = lambda param: param[self.indexDict['P_PERM_K']] * self.data['F_RT'] * init_vals["bV"] \
				/ (np.exp(self.data['F_RT'] * init_vals["bV"]) - 1) \
				* (self.data['B_ACT_K'] - init_vals["cK"] * np.exp(self.data['F_RT'] * init_vals["bV"]))

		# Eq 4.7
		par_J_pump = lambda param: param[self.indexDict['J_Pump_max']] * (init_vals['cNa']
											/ (init_vals['cNa']
											+ self.data['K_Na_In_pump'] * (1 + init_vals['cK']
																			/ self.data['K_K_in_pump']))) ** 3 \
				* (self.data['B_ACT_K'] / (self.data['B_ACT_K'] +
											self.data['K_K_ext_pump'] * (
													1 + self.data['B_ACT_NA'] / self.data['K_Na_ext_pump']))) ** 2

		# Eq 4.10
		par_pJ_Na = lambda param: param[self.indexDict['P_PERM_NA']] * self.data['F_RT'] * init_vals["tV"] \
				/ (np.exp(self.data['F_RT'] * init_vals["tV"]) - 1) \
				* (self.data['B_CONS_NA'] - init_vals["aNa"] * np.exp(
				self.data['F_RT'] * init_vals["tV"]))

		# Eq 4.11
		par_pJ_K = lambda param: param[self.indexDict['P_PERM_K']] * self.data['F_RT'] * init_vals["tV"] \
				/ (np.exp(self.data['F_RT'] * init_vals["tV"]) - 1) \
				* (self.data['B_CONS_K'] - init_vals["aK"] * np.exp(self.data['F_RT'] * init_vals["tV"]))

		# Eq 4.12
		par_pJ_Cl = lambda param: - param[self.indexDict['P_PERM_CL']] * self.data['F_RT'] * init_vals["tV"] \
				/ (np.exp(self.data['F_RT'] * init_vals["tV"]) - 1) \
				* (init_vals["aCl"] - self.data['B_CONS_CL'] * np.exp(
				self.data['F_RT'] * init_vals["tV"]))


		# Eq 3.2.1
		par_dcN_Na = lambda parameters: par_aJ_Na(parameters) + init_vals['J_co'] - 3 * par_J_pump(parameters)

		# Eq 3.2.2
		par_dcN_K = lambda parameters: par_aJ_K(parameters) + par_bJ_K(parameters) + init_vals['J_co'] + 2 * \
				par_J_pump(parameters)

		# Eq 3.2.3
		par_dcN_Cl = lambda parameters: par_aJ_Cl(parameters) + par_bJ_Cl(parameters) + 2 * init_vals['J_co']


		# Eq 7.1
		par_pI = lambda param: - par_pJ_Na(param) + par_pJ_Cl(param) - par_pJ_K(param)

		# Eq 7.2
		par_aI = lambda param: - par_aJ_Na(param) + par_aJ_Cl(param) - par_aJ_K(param)

		# Eq 7.3
		par_bI = lambda param: par_J_pump(param) + par_bJ_Cl(param) - par_aJ_K(param)


		# Eq 6.1
		par_daV = lambda parameters: par_pI(parameters) - par_aI(parameters)

		par_dbV = lambda parameters: par_pI(parameters) + par_bI(parameters)

		# How many dimensions the system has.
		dims = len(self.variable_params)

		# This is the system of 5 euqatons that is supposed to be 0
		def system(p):
			f = np.zeros(5, dtype=np.float64)
			f[0] = par_dcN_Na(p)	# \frac{dN_{Na}^c}{dt}
			f[1] = par_dcN_K(p)		# \frac{dN_K^c}{dt}
			f[2] = par_dcN_Cl(p)	# \frac{dN_{Cl}^c}{dt}
			f[3] = par_daV(p)		# \frac{dV_a}{dt}
			f[4] = par_dbV(p)		# \frac{dV_b}{dt}
			return f
		
		param = least_squares(system, x0=self.variable_params, bounds=(np.zeros(dims), np.full(dims, np.inf)))

		for index, key in enumerate(self.indexDict):
			self.data[key] = param.x[index]

		return param
		





# A class with all the necessary functions to run an ASL Model.
class AirwayModel:
	"""
	Basically a Pandas.DataFrame with all of the necessary functions to
	run the simulation.

	It's recommended to use with a for loop iterating through step in range(init_data['max_steps'])
	where you run AirwayModel.fn_run(step). This is because the model which this class is meant to be used on
	updates the values based on the previous one.

	__init__(self, init_data):
		:arg init_data:
				* a dictionary with:
					:var 'max_steps': How long with the simulation be.
					:var 'time_frame': How many seconds each step is supposed to represent.
					:var 'H': initial water height (meters).
					:var Molar concentration for all the relevant ions (mM).
					:var Concentrations for all nucleotides (ATP, ..., ADO, INO) in muM (micro-moles).
				* path to a .csv:
					It will be used to load the DataFrame from the .csv with the pandas.read_csv() function.
				* DEFAULT: None (it will ask the user for values).
		:exception if there's an issue loading the values, they will be asked manually.
	fn_run(self, step):
		:arg step=1: is the step to load the values to, by default will load the values after the initial.
		:returns nothing
	"""

	# Complete with all variables that will be tracked as column names
	variables = ('Time (min)', 'H', 'H_c', 'dH', 'OSM_a', 'OSM_c', 'aNa', 'aCl', 'aK', 'cNa', 'cCl', 'cK',
				 'aJ_Na', 'pJ_Na', 'aJ_K', 'pJ_K', 'aJ_Cl_CaCC', 'aJ_Cl_CFTR', 'pJ_Cl', 'J_co', 'J_pump',
				 'bJ_Cl', 'bJ_K', 'daV', 'dbV', 'aV', 'bV', 'tV', 'aI', 'bI', 'pI',
				 'p_CaCC', 'p_CFTR', 'p_ENaC', 'p_BK', 'p_CaKC',
				 'ATP', 'ADP', 'AMP', 'ADO', 'INO')

	# Separated variables into sections
	apical_ions = ('aNa', 'aCl', 'aK')
	cel_ions = ('cNa', 'cCl', 'cK')
	apical_flow = ('aJ_Na', 'aJ_K', 'aJ_Cl_CaCC', 'aJ_Cl_CFTR')
	paracellular_flow = ('pJ_Na', 'pJ_K', 'pJ_Cl')
	basolateral_flow = ('J_co', 'J_pump', 'bJ_Cl', 'bJ_K')
	voltage = ('aV', 'bV', 'tV')
	current = ('aI', 'bI', 'pI')
	permeabilities = ('p_CaCC', 'p_CFTR', 'p_ENaC', 'p_BK', 'p_CaKC')
	nucleotide = ('ATP', 'ADP', 'AMP', 'ADO', 'INO')

	sections = ('H', 'H_c', 'dH', apical_ions, cel_ions, apical_flow, paracellular_flow, basolateral_flow,
				voltage, current, permeabilities, nucleotide)

	sections_txt = ('Height', 'Cell Height', 'dH', 'Apical-Ions', 'Cell-Ions', 'Apical-Flow', 'Paracellular-Flow',
					'Basolateral-Flow', 'Voltage', 'Current', 'Permeability', 'Nucleotides')

	# The units correspond to the sections_txt
	sections_units = ('m', 'm', 'm/s', 'mM', 'mM', 'mol/(m^2 sec)', 'mol/(m^2 sec)',
					  'mol/(m^2 sec)', 'V', 'A / m^2', 'm/s', 'micro M (muM)')

	# Required initial values
	initial_vars = ('max_steps', 'time_frame', 'H', 'CF') + apical_ions + cel_ions + nucleotide

	# Text to display when asking for the variables
	ions_user_friendly_text = {'aNa': "apcial sodium", 'aCl': "apical chloride", 'aK': "apical potassium",
						 'cNa': "cellular sodium", 'cCl': "cellular chloride", 'cK': "cellular potassium"}

	# Initial conditions setup
	def __init__(self, init_data=None, special_constants=None):
		"""
		Will initialize the tracked values. If not successfully from init_data, it will prompt the user
		add them manually.

		:arg init_data:
				* a dictionary with:
					:var 'max_steps': How long with the simulation be.
					:var 'time_frame': How many seconds each step is supposed to represent.
					:var 'H': initial water height (meters).
					:var Molar concentration for all the relevant ions (mM).
					:var Concentrations for all nucleotides (ATP, ..., ADO, INO) in muM (micro-moles).
				* path to a .csv:
					It will be used to load the DataFrame from the .csv with the pandas.read_csv() function.
				* DEFAULT: None (it will ask the user for values).
		:arg special_constants (Optional)
				Dictionary with the name of the parameters to overwite. If it is set, it will not calculate
					the parameters that fit the steady state.
		:exception if there's an issue loading the values, they will be asked manually.
		"""

		# Will create the needed dictionary to init the variables
		if isinstance(init_data, str) and '.csv' in init_data:
			self.init_pd = pd.read_csv(init_data)
			init_data = {}

			if all(k in self.init_pd for k in self.initial_vars):
				for key in self.initial_vars:
					if key == 'max_steps':
						init_data[key] = int(self.init_pd[key][0])
					elif key == 'CF':
						init_data[key] = bool(self.init_pd[key][0])
					else:
						init_data[key] = float(self.init_pd[key][0])
			
			for var in self.voltage:
				if var in self.init_pd:
					init_data[var] = float(self.init_pd[var][0])

			del(self.init_pd)

		# If the right argument was passed or it was rightly processed, start the init_vars
		if isinstance(init_data, dict) and all(k in init_data for k in self.initial_vars):
			self.max_steps = int(init_data['max_steps'])
			# Here we are creating the data-frame from the list of 'self.variables'
			self.data = pd.DataFrame(np.zeros((self.max_steps, len(self.variables))), columns=self.variables)

			self.time_frame = float(init_data['time_frame'])
			self.isCF = bool(init_data['CF'])

			self.co = Constants(self.isCF, constants=special_constants)
			
			self.data['H'][0] = float(init_data['H'])

			for ion in self.apical_ions + self.cel_ions:
				self.data[ion][0] = self.to_activity(np.absolute(float(init_data[ion])))

			for nuc in self.nucleotide:
				self.data[nuc][0] = float(init_data[nuc])

			
			if all(v in init_data for v in ('aV', 'bV')):
				self.data['aV'][0] = float(init_data['aV'])
				self.data['bV'][0] = float(init_data['bV'])

				if 'pV' in init_data:
					self.data['pV'][0] = float(init_data['pV'])
				else:
					self.data["tV"][0] = self.data["bV"][0] - self.data["aV"][0]
			else:
				self.data["aV"][0] = self.fn_voltage('a')
				self.data["bV"][0] = self.fn_voltage('b')
				self.data["tV"][0] = self.data["bV"][0] - self.data["aV"][0]
		else:
			self.inputVals(special_constants)
		
		if special_constants is None:
			self.data["J_co"][0] = self.fn_J_co()
			self.co.getSteadyVals(self.data.head(1))

		# Will add all of the accessibility methods
		self.to_csv = self.data.to_csv
		self.drop = self.data.drop

		# Set the fn_run(), constants, and dependent-functions to the corresponding mode, CF or NL
		if self.isCF:
			self.fn_p_CFTR = lambda step=1: 0
			self.fn_p_ENaC = self.fn_p_ENaC_CF
		else:
			self.fn_p_CFTR = self.fn_p_CFTR_NL
			self.fn_p_ENaC = self.fn_p_ENaC_NL

		# If the first line of the data isn't full, then load all of the initial values.
		if not self.isFull():
			self.data['Time (min)'][0] = 0
			self.data["OSM_a"][0] = self.fn_OSM_a()
			self.data["OSM_c"][0] = self.fn_OSM_c()
			self.data["H_c"][0] = self.co['CELL_H']
			self.data["dH"][0] = self.fn_dH()

			# The voltage was calulated before.

			self.data["p_CaCC"][0] = self.fn_p_CaCC()
			self.data["p_CFTR"][0] = self.fn_p_CFTR()
			self.data["p_ENaC"][0] = self.fn_p_ENaC()
			self.data["p_BK"][0] = self.fn_p_BK()
			self.data["p_CaKC"][0] = self.fn_p_CaKC()

			self.data["aJ_Na"][0] = self.fn_aJ_Na()
			self.data["pJ_Na"][0] = self.fn_pJ_Na()
			self.data["aJ_K"][0] = self.fn_aJ_K()
			self.data["pJ_K"][0] = self.fn_pJ_K()
			self.data["bJ_K"][0] = self.fn_bJ_K()
			self.data["aJ_Cl_CaCC"][0] = self.fn_aJ_Cl_CaCC()
			self.data["aJ_Cl_CFTR"][0] = self.fn_aJ_Cl_CFTR()
			self.data["pJ_Cl"][0] = self.fn_pJ_Cl()
			self.data["bJ_Cl"][0] = self.fn_bJ_Cl()
			self.data["J_co"][0] = self.fn_J_co()
			self.data["J_pump"][0] = self.fn_J_pump()

			self.data["aI"][0] = self.fn_aI()
			self.data["bI"][0] = self.fn_bI()
			self.data["pI"][0] = self.fn_pI()
			self.data["daV"][0] = self.fn_daV()
			self.data["dbV"][0] = self.fn_dbV()

	def __getitem__(self, item):
		return self.data[item]

	def __repr__(self):
		"""Since the main object is a the pd.DataFrame, it will be what's printed"""
		return 'ASL Model\n' + str(self.data)

	def __str__(self):
		"""When you're working with the model, the important thing is the data, it will print the
		self.data object which is a pd.DataFrame"""
		return str(self.data)

	def __len__(self):
		return len(self.data)
		
	def to_cons(self, activity):
		"""Return: activity / self.co['GAMMA'] = concentration"""
		return activity / self.co["GAMMA"]

	def to_activity(self, concentration):
		"""Return: activity = concentration * self.co['GAMMA']"""
		return concentration * self.co["GAMMA"]

	def isFull(self):
		"""
		Will check whether the DataFrame was fully loaded into the INITIAL values, there are
		no initial values set to 0.
		:return: True if all the data is complete, false otherwise.
		"""
		for k in self.data:
			if k not in self.variables or k == 'Time (min)' or (k == 'p_CFTR' and self.isCF):
				pass
			elif self.data[k][0] == 0:
				return False
		return True

	def inputVals(self, special_constants=None):
		"""
		This function is used to ask for all of the necessary initial values. It will go one by one
		of the variables in AirwayModel.initial_vars asking the user to manually input then and then
		load it to the corresponding variable.

		:arg special_constants	Dictionary with values with wich to overwrite the default constants.
		"""
		self.max_steps = int(input('Enter how many steps the run will take: '))
		# Here we are creating the data-frame from the list of 'self.variables'
		self.data = pd.DataFrame(np.zeros((self.max_steps, len(self.variables))), columns=self.variables)

		tmp = input('Enter if simulating Cystic Fibrosis (write CF or 1) or Normal Conditions (write NL or 0).\n'
					'\tAny input error will count as Normal Conditions: ')
		if '1' in tmp or ('c' in tmp.lower() and 'f' in tmp.lower()) or 'true' == tmp.lower():
			self.isCF = True
			print('\t\tCystic Fibrosis activated.')
		else:
			self.isCF = False
			print('\t\tNormal Conditions activated.')

		self.co = Constants(self.isCF, constants=special_constants)

		self.time_frame = float(input('Enter how many seconds each step represents: '))

		self.data['H'][0] = float(input('Enter the ASL height in meters: '))

		for ion in self.apical_ions + self.cel_ions:
			self.data[ion][0] = np.absolute(float(input('Enter the value for ' + self.ions_user_friendly_text[ion] + ' (in mM): ')))
		for nucl in self.nucleotide:
			self.data[nucl][0] = np.absolute(float(input('Enter the value for ' + nucl + ' (in micro-M): ')))
		
		aV = input('Enter the value for apical voltage in Volts (leave blank for it to be calculated): ')
		bV = input('Enter the value for basolateral voltage in Volts (leave blank for it to be calculated): ')
		tV = input('Enter the value for transmembrane (paracellular) voltage in Volts (leave blank for it to be calculated): ')
		
		if aV == '':
			self.data["aV"][0] = self.fn_voltage('a')
		else:
			self.data['aV'][0] = float(aV)

		if bV == '':
			self.data["bV"][0] = self.fn_voltage('b')
		else:
			self.data['bV'][0] = float(bV)
		
		if tV == '':
			self.data["tV"][0] = self.data["bV"][0] - self.data["aV"][0]
		else:
			self.data['pV'][0] = float(tV)

	# Get section, will return a dictionary with the data on the step.
	def getSection(self, variables, steps):
		"""
		Will return a dictionary containing the variables as the key
		and the values as the definition of the key
		
		:param variables: Is a list of all the variables to load form self.data. It could also be just one variable.
				But it must be one of the tracked variables, as it needs to be inside self.data
		:param steps: Could be a list or a single step to get the data from. NOT a tuple with the min and max steps.
		"""

		# Dictionary that will contain all of the necessary sections and their steps.
		d = {}
		try:  # to iterate over the sections.
			for variable in variables:
				try:  # to iterate over the steps
					ls = []
					for step in steps:
						ls.append(self.data[variable][step])
					d[variable] = ls
				except TypeError:
					d[variables] = self.data[variables][steps]
		except TypeError:  # if you can't iterate over the sections.
			try:  # to iterate over the steps
				ls = []
				for step in steps:
					ls.append(self.data[variables][step])
				d[variables] = ls
			except TypeError:  # if you can't iterate over the steps
				d[variables] = self.data[variables][steps]
		return d

	# Voltage Calculator, approximate from Nernst Eq.
	def fn_voltage(self, membrane='a', step=1):
		"""
		Uses the Nernst equation to estimate the voltage of 'membrane' at 'step'.
		Uses the main ions in the compartment 'a' or 'b' depending on the membrane and the ones inside the cell.
		Main ions are: Chloride, Sodium, and Potassium.
		:param membrane: The membrane to which to calculate it's voltage, if the transmembrane potential is asked
				('t' or 'p') it will do so by calculating the difference of the apical and basolateral potentials.
		:param step: The step to which to calculate the potential to, it will again use the values from step-1 to
				calculate it's results.
		:return: It will return the estimated potential in Volts.
		"""
		if membrane == 'a':
			v = 0
			for ion in ['Na', 'K']:
				v += np.log(self.data[membrane + ion][step - 1] / self.data['c' + ion][step - 1]) / self.co['F_RT']
			for ion in ['Cl']:
				v -= np.log(self.data[membrane + ion][step - 1] / self.data['c' + ion][step - 1]) / self.co['F_RT']
			return v
		elif membrane == 'b':
			return np.log(self.co['B_CONS_NA'] / self.data['cNa'][step - 1]) / self.co['F_RT'] \
				   + np.log(self.co['B_CONS_K'] / self.data['cK'][step - 1]) / self.co['F_RT'] \
				   - np.log(self.co['B_CONS_CL'] / self.data['cCl'][step - 1]) / self.co['F_RT']
		elif membrane == 'p' or membrane == 't':
			return self.fn_voltage('b', step) - self.fn_voltage('a', step)
		else:
			raise KeyError('Must choose either apical, basolateral or transmembrane (paracellular).')

	# Run Eq
	def fn_run(self, step=1):
		"""
		Will run all the required Eqs and update the values.
		:param step: Is the step to which to calculate the values to. When required to obtain a value, it will
				collect it from step-1.
				Defaults to 1
		"""

		self.data['Time (min)'][step] = step * self.time_frame / 60
		self.data["dH"][step] = self.fn_dH(step)
		self.data["H"][step] = self.data["H"][step - 1] + self.fn_dH(step) * self.time_frame
		self.data["H_c"][step] = self.data["H_c"][step - 1] + self.fn_dH_c(step) * self.time_frame
		self.data["OSM_a"][step] = self.fn_OSM_a(step)
		self.data["OSM_c"][step] = self.fn_OSM_c(step)

		self.data["aNa"][step] = self.data["aNa"][step - 1] + self.fn_daN_Na(step) * self.time_frame \
								 / self.data['H'][step - 1]	 # in mM
		self.data["aCl"][step] = self.data["aCl"][step - 1] + self.fn_daN_Cl(step) * self.time_frame \
								 / self.data['H'][step - 1]	 # in mM
		self.data["aK"][step] = self.data["aK"][step - 1] + self.fn_daN_K(step) * self.time_frame \
								/ self.data['H'][step - 1]	# in mM

		self.data["cNa"][step] = self.data["cNa"][step - 1] + self.fn_dcN_Na(step) * self.time_frame \
								 / self.data['H_c'][step - 1]
		self.data["cCl"][step] = self.data["cCl"][step - 1] + self.fn_dcN_Cl(step) * self.time_frame \
								 / self.data['H_c'][step - 1]
		self.data["cK"][step] = self.data["cK"][step - 1] + self.fn_dcN_K(step) * self.time_frame \
								/ self.data['H_c'][step - 1]

		self.data["aJ_Na"][step] = self.fn_aJ_Na(step)
		self.data["aJ_Cl_CaCC"][step] = self.fn_aJ_Cl_CaCC(step)
		self.data["aJ_Cl_CFTR"][step] = self.fn_aJ_Cl_CFTR(step)
		self.data["aJ_K"][step] = self.fn_aJ_K(step)
		self.data["pJ_Na"][step] = self.fn_pJ_Na(step)
		self.data["pJ_Cl"][step] = self.fn_pJ_Cl(step)
		self.data["pJ_K"][step] = self.fn_pJ_K(step)

		self.data["J_co"][step] = self.fn_J_co(step)
		self.data["J_pump"][step] = self.fn_J_pump(step)
		self.data["bJ_Cl"][step] = self.fn_bJ_Cl(step)
		self.data["bJ_K"][step] = self.fn_bJ_K(step)

		self.data["daV"][step] = self.fn_daV(step)
		self.data["aV"][step] = self.data["aV"][step - 1] + self.data["daV"][step] * self.time_frame
		self.data["dbV"][step] = self.fn_dbV(step)
		self.data["bV"][step] = self.data["bV"][step - 1] + self.data["dbV"][step] * self.time_frame
		self.data["tV"][step] = self.data["bV"][step - 1] - self.data["aV"][step]
		self.data["aI"][step] = self.fn_aI(step)
		self.data["pI"][step] = self.fn_pI(step)
		self.data["bI"][step] = self.fn_bI(step)

		self.data["p_CaCC"][step] = self.fn_p_CaCC(step)
		self.data["p_CFTR"][step] = self.fn_p_CFTR(step)
		self.data["p_ENaC"][step] = self.fn_p_ENaC(step)
		self.data["p_BK"][step] = self.fn_p_BK(step)
		self.data["p_CaKC"][step] = self.fn_p_CaKC(step)

		self.data["ATP"][step] = self.data["ATP"][step - 1] + self.fn_dATP(step) * self.time_frame
		self.data["ADP"][step] = self.data["ADP"][step] + self.fn_dADP(step) * self.time_frame
		self.data["AMP"][step] = self.data["AMP"][step - 1] + self.fn_dAMP(step) * self.time_frame
		self.data["ADO"][step] = self.data["ADO"][step - 1] + self.fn_dADO(step) * self.time_frame
		self.data["INO"][step] = self.data["INO"][step - 1] + self.fn_dINO(step) * self.time_frame

	def fn_runAll(self):
		"""
		Very self explanatory by the name
		it will run the simulation from step 1 till the max ammount of steps
		"""
		
		for step in range(1, self.max_steps):
			self.fn_run(step)

	def setParams(self):
		"""
		This will be used to estimate the correct parameters for the simulation.
		It will be vary basic, go through all the values (as guiven in the paper) between the mean
		and two standard deviations. And will use the best set of parameters that fits the steady state.

		We will need
		J_Na^a + J_Na^b = 0
		J_Cl^a + J_Cl^b = 0
		J_K^a + J_K^b = 0
		-I_Cl^a + I_Na^a + I_K^a = 0
		-I_Cl^b + I_Na^b + I_K^b = 0

		Coeficent of variation = n/mean
		we need to cv to get the value
		So range(mean - mean * 2 * cv, mean + mean * 2 * cv)
			Might need want to create a generator for this speciefic range
		nested for loops
		Will need to create a vector or DataFrame to store all of them beforehand
		"""
		pass

	def fn_aJ_H2O(self, step=1):
		return self.co['A_PERM_H20'] * (self.data["OSM_c"][step - 1] - self.data["OSM_a"][step - 1])

	def fn_pJ_H2O(self, step=1):
		return self.co['P_PERM_H2O'] * (self.data["OSM_a"][step - 1] - self.co['B_OSM'])

	def fn_bJ_H2O(self, step=1):
		return self.co['B_PERM_H20'] * (self.data["OSM_c"][step - 1] - self.co['B_OSM'])

	# Eq 1
	def fn_dH(self, step=1):
		"""Calculates the rate of change of the apical water height in meters per second (m/s)"""
		return self.co['A_V_H2O'] * (self.fn_pJ_H2O(step) - self.fn_aJ_H2O(step))

	def fn_dH_c(self, step=1):
		return self.co['A_V_H2O'] * (self.fn_aJ_H2O(step) + self.fn_bJ_H2O(step))

	# Eq 2.1.2
	def fn_OSM_a(self, step=1):
		"""Will calculate the osmomolarity of the Airway Liquid at step in mili-molar or moles per meter cubed"""
		return self.co['PHI'] / self.co['GAMMA'] * (self.data["aNa"][step - 1] + self.data["aK"][step - 1]
													+ self.data["aCl"][step - 1] + self.co['A_ACT_OI'])

	# Eq 2.2.2
	def fn_OSM_c(self, step=1):
		"""Will calculate the osmomolarity of the cell at step in mili-molar or moles per meter cubed"""
		return self.co['PHI'] / self.co['GAMMA'] * (
				self.data["cNa"][step - 1] + self.data["cK"][step - 1] + self.data["cCl"][step - 1] + self.co[
			'C_ACT_OI'])

	# Eq 3.1.1
	def fn_daN_Na(self, step=1):
		"""Difference of flow of Sodium across the apical membrane"""
		return -self.data['aJ_Na'][step - 1] + self.data['pJ_Na'][step - 1]

	# Eq 3.1.2
	def fn_daN_K(self, step=1):
		"""Difference of flow of Potassium across the apical membrane"""
		return -self.data['aJ_K'][step - 1] + self.data['pJ_K'][step - 1]

	# Eq 3.1.3
	def fn_daN_Cl(self, step=1):
		"""Difference of flow of Chloride across the apical membrane"""
		return -(self.data['aJ_Cl_CaCC'][step - 1] + self.data['aJ_Cl_CFTR'][step - 1]) + self.data['pJ_Cl'][step - 1]

	# Eq 3.2.1
	def fn_dcN_Na(self, step=1):
		"""Difference of flow of Sodium across the basolateral membrane"""
		return self.data['aJ_Na'][step - 1] + self.data['J_co'][step - 1] - 3 * self.data['J_pump'][step - 1]

	# Eq 3.2.2
	def fn_dcN_K(self, step=1):
		"""Difference of flow of Potassium across the basolateral membrane"""
		return self.data['aJ_K'][step - 1] + self.data['bJ_K'][step - 1] + self.data['J_co'][step - 1] + 2 * \
			   self.data['J_pump'][step - 1]

	# Eq 3.2.3
	def fn_dcN_Cl(self, step=1):
		"""Difference of flow of Chloride across the basolateral membrane"""
		return self.data['aJ_Cl_CaCC'][step - 1] + self.data['aJ_Cl_CFTR'][step - 1] \
			   + self.data['bJ_Cl'][step - 1] + 2 * self.data['J_co'][step - 1]

	# Eq 4.1
	def fn_aJ_Na(self, step=1):
		"""
		Will return the flow of Sodium ions across the apical membrane.

		:return: Flow moles per meter squared per second, mol / (m^2 sec)
		"""
		return self.data["p_ENaC"][step - 1] * self.data["aV"][step - 1] * self.co['F_RT'] \
			   * (self.data["aNa"][step - 1] -
				  self.data["cNa"][step - 1] * np.exp(self.data["aV"][step - 1] * self.co['F_RT'])) \
			   / (np.exp(self.data["aV"][step - 1] * self.co['F_RT']) - 1)

	# Eq 4.2
	def fn_aJ_K(self, step=1):
		"""
		:return: Flow of potassium across the apical membrane at a given step (moles / (m^2 sec)).
		"""
		return self.data["p_BK"][step - 1] * self.data["aV"][step - 1] * self.co['F_RT'] / (
				np.exp(self.data["aV"][step - 1] * self.co['F_RT']) - 1) \
			   * (self.data["aK"][step - 1] - self.data["cK"][step - 1] * np.exp(
			self.co['F_RT'] * self.data["aV"][step - 1]))

	# Eq 4.3
	def fn_aJ_Cl(self, step=1):
		"""
		:return: Flow of chloride across the apical membrane at a given step (moles / (m^2 sec)).
		"""
		return -(self.data["p_CaCC"][step - 1] + self.data["p_CFTR"][step - 1]) \
			   * self.co['F_RT'] * self.data["aV"][step - 1] / (np.exp(self.co['F_RT'] * self.data["aV"][step - 1]) - 1) \
			   * (self.data["cCl"][step - 1] - self.data["aCl"][step - 1] * np.exp(
			self.co['F_RT'] * self.data["aV"][step - 1]))

	def fn_aJ_Cl_CaCC(self, step=1):
		"""
		:return: The apical flow of chloride through the CaCC channel.
		"""
		return -self.data["p_CaCC"][step - 1] \
			   * self.co['F_RT'] * self.data["aV"][step - 1] / (np.exp(self.co['F_RT'] * self.data["aV"][step - 1]) - 1) \
			   * (self.data["cCl"][step - 1] - self.data["aCl"][step - 1] * np.exp(
					self.co['F_RT'] * self.data["aV"][step - 1]))

	def fn_aJ_Cl_CFTR(self, step=1):
		"""
		:return: The apical flow of chloride ions through the CFTR channel.
		"""
		return -self.data["p_CFTR"][step - 1] \
			   * self.co['F_RT'] * self.data["aV"][step - 1] / (np.exp(self.co['F_RT'] * self.data["aV"][step - 1]) - 1) \
			   * (self.data["cCl"][step - 1] - self.data["aCl"][step - 1] * np.exp(
			self.co['F_RT'] * self.data["aV"][step - 1]))

	# Eq 4.5
	def fn_bJ_Cl(self, step=1):
		"""moles per second per meter squared - moles / (sec m^2)"""
		return - self.co['P_PERM_CL'] * self.co['F_RT'] * self.data["bV"][step - 1] \
			   / (np.exp(self.co['F_RT'] * self.data["bV"][step - 1]) - 1) \
			   * (self.data['cCl'][step - 1] - self.co['B_ACT_CL'] * np.exp(
			self.co['F_RT'] * self.data["bV"][step - 1]))

	# Eq 4.6
	def fn_bJ_K(self, step=1):
		"""moles per second per meter squared - moles / (sec m^2)"""
		return self.co['P_PERM_K'] * self.co['F_RT'] * self.data["bV"][step - 1] \
			   / (np.exp(self.co['F_RT'] * self.data["bV"][step - 1]) - 1) \
			   * (self.co['B_ACT_K'] - self.data["cK"][step - 1] * np.exp(self.co['F_RT'] * self.data["bV"][step - 1]))

	# Eq 4.7
	def fn_J_pump(self, step=1):
		"""moles per second per meter squared - moles / (sec m^2)"""
		return self.co['J_Pump_max'] * (self.data['cNa'][step - 1]
										/ (self.data['cNa'][step - 1]
										   + self.co['K_Na_In_pump'] * (1 + self.data['cK'][step - 1]
																		/ self.co['K_K_in_pump']))) ** 3 \
			   * (self.co['B_ACT_K'] / (self.co['B_ACT_K'] +
										self.co['K_K_ext_pump'] * (
												1 + self.co['B_ACT_NA'] / self.co['K_Na_ext_pump']))) ** 2

	# NKCC2 Z(1-15)
	def fn_Z_nkcc(self, step=1):
		"""Will return a tuple with all the 'Z's on the Benjamin-Johnson model from 1997. All 16 of them"""
		return (self.co['COT_Z'][0] * self.data['cCl'][step - 1],
				self.co['COT_Z'][1] * self.co['B_ACT_NA'],
				self.co['COT_Z'][2] * self.data['cCl'][step - 1] * self.data['cK'][step - 1],
				self.co['COT_Z'][3] * self.co['B_ACT_CL'] * self.co['B_ACT_K'],
				self.co['COT_Z'][4] * self.data['cCl'][step - 1] ** 2 * self.data['cK'][step - 1],
				self.co['COT_Z'][5] * self.co['B_ACT_CL'] * self.co['B_ACT_K'] * self.co['B_ACT_NA'],
				self.co['COT_Z'][6] * self.data['cCl'][step - 1] ** 2 * self.data['cK'][step - 1] * self.data['cNa'][
					step - 1],
				self.co['COT_Z'][7] * self.co['B_ACT_CL'] ** 2 * self.co['B_ACT_K'] * self.co['B_ACT_NA'],
				self.co['COT_Z'][8] * self.data['cCl'][step - 1] ** 2 * self.data['cK'][step - 1]
				* self.data['cNa'][step - 1] * self.co['B_ACT_NA'],
				self.co['COT_Z'][9] * self.data['cCl'][
					step - 1] * self.co['B_ACT_CL'] ** 2 * self.co['B_ACT_K'] * self.co['B_ACT_NA'],
				self.co['COT_Z'][10] * self.data['cCl'][step - 1] ** 2 * self.data['cK'][step - 1]
				* self.data['cNa'][step - 1] * self.co['B_ACT_CL'] * self.co['B_ACT_NA'],
				self.co['COT_Z'][11] * self.data['cCl'][step - 1] * self.data['cK'][step - 1]
				* self.co['B_ACT_CL'] ** 2 * self.co['B_ACT_K'] * self.co['B_ACT_NA'],
				self.co['COT_Z'][12] * self.data['cCl'][step - 1] ** 2 * self.data['cK'][step - 1]
				* self.co['B_ACT_CL'] ** 2 * self.co['B_ACT_K'] * self.co['B_ACT_NA'],
				self.co['COT_Z'][13] * self.data['cCl'][step - 1] ** 2 * self.data['cK'][step - 1]
				* self.data['cNa'][step - 1] * self.co['B_ACT_CL'] * self.co['B_ACT_K'] * self.co['B_ACT_NA'],
				self.co['COT_Z'][14] * self.data['cCl'][step - 1] ** 2 * self.data['cK'][step - 1]
				* self.data['cNa'][step - 1] * self.co['B_ACT_CL'] ** 2 * self.co['B_ACT_K'] * self.co['B_ACT_NA'],
				self.co['COT_Z'][15])

	# Eq 4.8
	def fn_J_co(self, step=1):
		"""
		Benjamin-Jonson Model of the cotransporter, as found on the thesis (MISSING REFERENCE)
		:return units where originally 10^4 mol / (m sec), that's why we have the 10000 multiplying the expression.
		"""
		return - self.co['COT_D'] * (
				self.co['COT_K_f_full'] * self.co['COT_K_f_empty'] * self.co['B_ACT_CL'] ** 2 * self.co['B_ACT_K'] *
				self.co['B_ACT_NA']
				- self.co['COT_K_b_full'] * self.co['COT_K_b_empty'] * self.data['cCl'][step - 1] ** 2
				* self.data['cNa'][step - 1] * self.data['cK'][step - 1]) \
			   / sum(self.fn_Z_nkcc(step))

	# Eq 4.10
	def fn_pJ_Na(self, step=1):
		"""Paracellular flow of for the Sodium ion, in moles per sec per meters squared - mole / (sec m^2)"""
		return self.co['P_PERM_NA'] * self.co['F_RT'] * self.data["tV"][step - 1] \
			   / (np.exp(self.co['F_RT'] * self.data["tV"][step - 1]) - 1) \
			   * (self.co['B_CONS_NA'] - self.data["aNa"][step - 1] * np.exp(
			self.co['F_RT'] * self.data["tV"][step - 1]))

	# Eq 4.11
	def fn_pJ_K(self, step=1):
		"""Paracellular flow of for the Potassium ion, in moles per sec per meters squared - mole / (sec m^2)"""
		return self.co['P_PERM_K'] * self.co['F_RT'] * self.data["tV"][step - 1] \
			   / (np.exp(self.co['F_RT'] * self.data["tV"][step - 1]) - 1) \
			   * (self.co['B_CONS_K'] - self.data["aK"][step - 1] * np.exp(self.co['F_RT'] * self.data["tV"][step - 1]))

	# Eq 4.12
	def fn_pJ_Cl(self, step=1):
		"""Paracellular flow of for the Chloride ion, in moles per sec per meters squared - mole / (sec m^2)"""
		return - self.co['P_PERM_CL'] * self.co['F_RT'] * self.data["tV"][step - 1] \
			   / (np.exp(self.co['F_RT'] * self.data["tV"][step - 1]) - 1) \
			   * (self.data["aCl"][step - 1] - self.co['B_CONS_CL'] * np.exp(
			self.co['F_RT'] * self.data["tV"][step - 1]))

	# Eq 5.1.1
	def fn_p_CaCC(self, step=1):
		"""
		Returns the permeability of the CaCC Channel.
		:return: in meters per second (m / sec)
		"""
		return self.co['CACC_PERM_MAX'] / (1 + self.co['CACC_ATP_AT_HALF_PERM'] / self.data["ATP"][step - 1])

	# Eq 5.1.2
	def fn_p_CFTR_NL(self, step=1):
		"""
		Returns the permeability of the CFTR channel
		:return: meters per second (m / sec)
		"""
		return self.co['CFTR_PERM_MAX'] / (1 + self.co['CFTR_ADO_AT_HALF_PERM'] / self.data["ADO"][step - 1] +
										   self.co['CFTR_ATP_AT_HALF_PERM'] / self.data["ATP"][step - 1])

	# Eq 5.2
	def fn_p_BK(self, step=1):
		"""
		Returns the permeability of the BK channel.
		:returns: in meters per second (m / sec).
		"""
		return self.co['BK_PERM_MAX'] / (1 + self.co['BK_ATP_AT_HALF_PERM'] / self.data["ATP"][step - 1])

	# Eq 5.3 for NL
	def fn_p_ENaC_NL(self, step=1):
		"""
		Is the permeability of the ENaC channel (Sodium)
		Now it's different for C.F.
		:return: The permeability for step in m/sec.
		"""
		return self.co['ENAC_PERM_MAX'] / (1 + self.data['ATP'][step - 1] / self.co['ENAC_ATP_AT_HALF_PERM']
										   + self.data['ADO'][step - 1] / self.co['ENAC_ADO_AT_HALF_PERM'])

	# Eq 5.3 for CF
	def fn_p_ENaC_CF(self, step=1):
		"""
		Is the permeability of the ENaC channel (Sodium)
		Now it's different for C.F.
		:return: The permeability for step in m/sec.
		"""

		return self.co['ENAC_PERM_MAX'] / (1 + self.data['ATP'][step - 1] / self.co['ENAC_ATP_AT_HALF_PERM'])

	# Eq 5.4
	def fn_p_CaKC(self, step=1):
		"""Permeability of the CaKC channel in meters per second (m/sec)"""

		return self.co['CAKC_PERM_MAX'] / (1 + self.co['CAKC_ATP_AT_HALF_PERM'] / self.data['ATP'][step - 1])

	# Eq 6.1
	def fn_daV(self, step=1):
		"""Will return the rate of change of the apical membrane's voltage. In Volts per second"""
		return (self.data["pI"][step - 1] - self.data["aI"][step - 1]) / (self.co['A_CAPACI'] * 100)

	# Eq 6.2
	def fn_dbV(self, step=1):
		"""Will return the rate of change of the basolateral membrane's voltage. In Volts per second"""
		return - (self.data["pI"][step - 1] + self.data["bI"][step - 1]) / (self.co['B_CAPACI'] * 100)

	# Eq 7.1
	def fn_pI(self, step=1):
		"""Will return paracellular current in A / m^2 or C / (sec m^2)"""
		return self.co['FARADAY'] * (- self.data["pJ_Na"][step - 1]
									 + self.data["pJ_Cl"][step - 1] - self.data["pJ_K"][step - 1])

	# Eq 7.2
	def fn_aI(self, step=1):
		"""Apical current in A / m^2 or C / (sec m^2)"""
		return self.co['FARADAY'] * (- self.data["aJ_Na"][step - 1]
									 + self.data["aJ_Cl_CaCC"][step - 1] + self.data["aJ_Cl_CFTR"][step - 1]
									 - self.data["aJ_K"][step - 1])

	# Eq 7.3
	def fn_bI(self, step=1):
		"""Basolateral Current in A / m^2 or C / (sec m^2)"""
		return self.co['FARADAY'] * (self.data["J_pump"][step - 1]
									 + self.data["bJ_Cl"][step - 1] - self.data["aJ_K"][step - 1])

	# Eq 8
	def fn_dATP(self, step=1):
		"""
		ATP rate of change in micro-Molar per second (muM / sec)
		"""

		return (self.co['J_atp'] - self.co['V_1_max'] * self.data["ATP"][step - 1] / (
				self.co['K_1_m'] + self.data["ATP"][step - 1])
				- self.co['V_2_max'] * self.data["ATP"][step - 1] / (self.co['K_2_m'] + self.data["ATP"][step - 1])
				- self.co['V_3_max'] * self.data["ATP"][step - 1] / (self.co['K_3_m'] + self.data["ATP"][step - 1])
				- self.co['V_10_max'] * self.data["ATP"][step - 1] / (self.co['K_10_m'] + self.data["ATP"][step - 1])
				- self.co['V_F_max'] / (1 +
										self.co['K_F_atp'] * self.co['K_F_amp'] /
										(self.data["ATP"][step - 1] * self.data["AMP"][step - 1])
										+ self.co['K_F_atp'] / self.data["ATP"][step - 1]
										+ self.co['K_F_amp'] / self.data["AMP"][step - 1])
				+ self.co['V_B_max'] / (1 + (self.co['K_B_adp'] / self.data["ADP"][step - 1]) ** 2
										+ 2 * self.co['K_B_adp'] / self.data["ADP"][step - 1])
				- self.data["ATP"][step - 1] * self.data["dH"][step - 1] / self.data["H"][step - 1]) / 60

	# Eq 9
	def fn_dADP(self, step=1):
		"""
		ADP rate of change in micro-Molar per second (muM / sec)
		"""

		return (self.co['J_adp'] / self.data['H'][step - 1]
				+ self.co['V_1_max'] * self.data["ATP"][step - 1] / (self.co['K_1_m'] + self.data["ATP"][step - 1])
				+ self.co['V_2_max'] * self.data["ATP"][step - 1] / (self.co['K_2_m'] + self.data["ATP"][step - 1])
				+ self.co['V_3_max'] * self.data["ATP"][step - 1] / (self.co['K_3_m'] + self.data["ATP"][step - 1])
				- self.co['V_4_max'] * self.data["ADP"][step - 1] / (self.co['K_4_m'] + self.data["ADP"][step - 1])
				- self.co['V_5_max'] * self.data["ATP"][step - 1] / (self.co['K_5_m'] + self.data["ATP"][step - 1])
				- 2 * self.co['V_F_max'] / (
						1 + self.co['K_F_atp'] * self.co['K_F_amp'] / (
						self.data["ATP"][step - 1] * self.data["AMP"][step - 1])
						+ self.co['K_F_atp'] / self.data["ATP"][step - 1] + self.co['K_F_amp'] / self.data["AMP"][
							step - 1])
				+ self.co['V_B_max'] / (1 + (self.co['K_B_adp'] / self.data["ADP"][step - 1]) ** 2)
				- 2 * self.co['V_B_max'] / (1 + self.co['K_B_adp'] / self.data["ADP"][step - 1])
				- self.data["ADP"][step - 1] * self.data['dH'][step - 1]) / 60

	# Eq 10
	def fn_dAMP(self, step=1):
		"""
		AMP's apical rate of change in micro-Molar per second (muM / sec)
		"""

		return (self.data['dH'][step - 1] * (self.data['AMP'][step - 1] / self.data['H'][step - 1])
				- self.co['V_8_max'] * self.data['AMP'][step - 1] / (
						self.co['K_8_m'] * (
						1 + self.data['ATP'][step - 1] / self.co['K_IN_atp'] + self.data['ADP'][
					step - 1] / self.co['K_IN_adp'])
						+ self.data['AMP'][step - 1])
				- self.co['V_7_max'] * self.data["AMP"][step - 1] / (
						self.co['K_7_m'] * (
						1 + self.data["ATP"][step - 1] / self.co['K_IN_atp']
						+ self.data["ADP"][step - 1] / self.co['K_IN_adp'] + self.data["AMP"][step - 1])
				)
				- self.co['V_6_max'] * self.data["AMP"][step - 1] / (
						self.co['K_6_m'] * (
						1 + self.data["ATP"][step - 1] / self.co['K_IN_atp']
						+ self.data["ADP"][step - 1] / self.co['K_IN_adp']
				) + self.data["AMP"][step - 1])
				+ self.co['V_B_max'] / (1 + self.co['K_B_adp'] / self.data["ADP"][step - 1])
				- self.co['V_F_max'] / (
						1 + self.co['K_F_amp'] * self.co['K_B_adp'] / (
						self.data["ADP"][step - 1] * self.data["AMP"][step - 1]) +
						self.co['K_F_atp'] / self.data["ATP"][step - 1] + self.co['K_F_amp'] / self.data["AMP"][
							step - 1])
				+ self.co['V_10_max'] * self.data["ATP"][step - 1] / (self.co['K_10_m'] + self.data["ATP"][step - 1])
				+ self.co['V_5_max'] * self.data["ATP"][step - 1] / (self.co['K_5_m'] + self.data["ATP"][step - 1])
				+ self.co['V_4_max'] / (self.co['K_4_m'] + self.data["ADP"][step - 1])
				+ self.co['J_amp']) / 60

	# Eq 11
	def fn_dADO(self, step=1):
		"""ADO's apical rate of change in micro-Molar per second (muM / sec)"""

		return (- self.data["ADO"][step - 1] * self.data["dH"][step - 1] / self.data["H"][step - 1]
				+ self.co['V_8_max'] * self.data["AMP"][step - 1] / (
						self.co['K_8_m'] * (
						1 + self.data["ATP"][step - 1] / self.co['K_IN_atp']
						+ self.data["ADP"][step - 1] / self.co['K_IN_adp']
				) + self.data["AMP"][step - 1])
				+ self.co['V_7_max'] * self.data["AMP"][step - 1] / (
						self.co['K_7_m'] * (
						1 + self.data["ATP"][step - 1] * self.data["ADP"][step - 1] / (
						self.co['K_IN_adp'] * self.co['K_IN_atp'])
						+ self.data["AMP"][step - 1])
				)
				+ self.co['V_6_max'] * self.data["AMP"][step - 1] / (
						self.co['K_6_m'] * (
						1 + self.data["ATP"][step - 1] / self.co['K_IN_atp']
						+ self.data["ADP"][step - 1] / self.co['K_IN_adp']
				) + self.data["AMP"][step - 1])
				- self.co['V_U1_max'] * self.data["ADO"][step - 1] / (
						self.co['K_U1_m'] + self.data["ADO"][step - 1]
						+ self.data["INO"][step - 1] * self.co['V_U1_max'] / self.co['K_U2_m']
				)
				- self.co['V_9_max'] * self.data["ADO"][step - 1] / (
						self.co['K_9_m'] + self.data["ADO"][step - 1])) / 60

	# Eq 12
	def fn_dINO(self, step=1):
		"""ADO's apical rate of change in micro-Molar per second (muM / sec)"""

		return ((self.co['V_9_max'] * self.data["ADO"][step - 1] / (self.co['K_9_m'] + self.data["ADO"][step - 1]))
				- self.co['V_U2_max'] * self.data["INO"][step - 1] / (
						self.co['K_U2_m'] + self.data["INO"][step - 1]
						+ self.data["ADO"][step - 1] * self.co['K_U2_m'] / self.co['K_U1_m'])
				- self.data["INO"][step - 1] * self.data["dH"][step - 1]) / 60
