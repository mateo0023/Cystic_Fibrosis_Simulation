import pandas as pd
import numpy as np


class Constants:
    """
    This is a class that will be used to store the values of the constants.
    """

    def __init__(self, cf=False, constants=None):
        self.CF = cf

        self.data = {}

        # General Naming: where_what_ofWhat

        self.data['A_V_H2O'] = 1.8e-6  # Molar volume of H2O. ------ m^3 mol^-1 (meters cube / moles)
        self.data['TEMP'] = 310.15  # Temperature is assumed to be constant. ------ K (degrees kelvin)

        self.data[
            'CELL_VOL'] = 1450  # Cell volume of water, used for init cell height. ------ (mu m)^3 (micro-meters cubed)
        self.data['CELL_H'] = (self.data['CELL_VOL'] / 2) ** (
                1 / 3) * 2e-6  # Cell shape assumed to be around 2NxNxN. ------ m (meters)

        # Permeability Constants
        self.data['P_PERM_H2O'] = 3.1e-5  # Permeability of water through the paracellular. ------ m/s
        self.data['A_PERM_H20'] = 2.4e-4  # Permeability of H2O through the apical. ------ m/s
        self.data['B_PERM_H20'] = 2.4e-5  # Permeability of H2O through the basolateral. ------ m/s

        # Osmomolarity & Concentration Constants
        self.data['B_OSM'] = 279.1  # Osmomolarity on the basolateral membrane. ------ mM
        self.data['PHI'] = 0.93
        self.data['GAMMA'] = 0.76
        self.data[
            'B_CONS_NA'] = 106.4  # Na+ concentration in the basolateral compartment. ------ mM (milli-moles / litter)
        self.data[
            'B_CONS_K'] = 4.0  # Potassium concentration in the basolateral compartment. ------ mM (milli-moles / litter)
        self.data[
            'B_CONS_CL'] = 91.2  # Chloride-ion concentration in the basolateral compartment. ------ mM (milli-moles / litter)
        self.data['A_CONS_NCL'] = 48.7  # Apical concentration on non-Chloride anions ------ mM
        self.data['A_CONS_IO'] = 2.7  # Concentration of impermeable osmolytes in the apical compartment. ------ mM (milli-moles / litter)
        self.data['C_CONS_NCL'] = 69.1  # Cellular concentration on non-Chloride anions ------ mM
        self.data['C_CONS_IO'] = 30.4  # Concentration of impermeable osmolytes in the cellular compartment. ------ mM (milli-moles / litter)
        self.data['B_ACT_NA'] = self.data['B_CONS_NA'] * self.data['GAMMA']  # Activity of the ion, see above.
        self.data['B_ACT_K'] = self.data['B_CONS_K'] * self.data['GAMMA']
        self.data['B_ACT_CL'] = self.data['B_CONS_CL'] * self.data['GAMMA']
        self.data['A_ACT_OI'] = (self.data['A_CONS_IO'] + self.data['A_CONS_NCL']
                                    * self.data['GAMMA'])  # Activity of the non-Chloride ions and the impermeable osmolytes.
        self.data['C_ACT_OI'] = (self.data['C_CONS_IO '] + self.data['C_CONS_NCL']
                                    * self.data['GAMMA'])  # Activity of the non-Chloride ions and the impermeable osmolytes.

        # Flow
        self.data['FARADAY'] = 96485  # Faraday's Constant. ------ C / mol
        self.data['R'] = 8.31447  # Gas constant. ------ J / (K mol)
        self.data['F_RT'] = self.data['FARADAY'] / (
                self.data['R'] * self.data['TEMP'])  # F/(R*T), To save computations.
        self.data['K_Na_In_pump'] = 0.99  # Look Ref 44 for more info. ------ mM
        self.data['K_K_in_pump'] = 9.1  # Look Ref 44 for more info.  -   mM
        self.data['K_K_ext_pump'] = 0.11  # Look Ref 44 for more info. ------ mM
        self.data['K_Na_ext_pump'] = 24.3  # Look Ref 44 for more info. ------ mM

        # REF 22
        # self.data['COT_K_1'] = 1000  # Cotransporter (Cot) off binding constant Na+. ------ mol^-1
        # self.data['COT_K_2'] = 0.7428  # Cot off binding constant Cl. ------ mol^-1
        # self.data['COT_K_03'] = 892.3  # Luminal off binding constant K. ------ mol^-1
        # self.data['COT_K_I3'] = 458.8  # Cytosolic off binding constant K. ------ mol^-1
        # self.data['COT_K_4'] = 46.71  # Cot of binding rate constant for Cl_2. ------ mol^-1
        # self.data['COT_K_FF'] = 1824  # Cot translocation rate constant. ------ s^-1
        # self.data['COT_K_BF'] = 2724  # Cot translocation rate constant. ------ s^-1
        # self.data['COT_K_FE'] = 18781  # Cot translocation rate constant. ------ s^-1
        # self.data['COT_K_BE'] = 13259  # ------ s^-1
        # self.data['COT_K_ON'] = 1e8  # ------ mol^-1 s^-1

        # NKCC2 Cotransporter - from Benjamin-Jonson (1997)
        self.data['COT_D'] = 0.4e-6  # Density of the NKCC Cotransporter - mol / m^2
        self.data['COT_K_Cl'] = 2.42  # mM
        self.data['COT_K_Na'] = 22.38  # mM
        self.data['COT_K_K'] = 234.74  # mM
        self.data['COT_K_f_empty'] = 37767  # s^-1
        self.data['COT_K_f_full'] = 1406  # s^-1
        self.data['COT_K_b_empty'] = 13196  # s^-1
        self.data['COT_K_b_full'] = 4025  # s^-1
        self.data['COT_Z'] = (self.data[''] * self.data['COT_K_K'] * self.data['COT_K_Na'] * self.data['COT_K_b_empty'],
                              self.data['COT_K_Cl'] ** 2 * self.data['COT_K_K'], self.data['COT_K_f_empty'],
                              self.data['COT_K_Cl'] * self.data['COT_K_Na'] * self.data['COT_K_b_empty'],
                              self.data['COT_K_Cl'] * self.data['COT_K_K'] * self.data['COT_K_f_empty'],
                              self.data['COT_K_Na'] * self.data['COT_K_b_empty'],
                              self.data['COT_K_Cl'] * self.data['COT_K_b_empty'],
                              self.data['COT_K_b_empty'] + self.data['COT_K_b_full'],
                              self.data['COT_K_f_empty'] + self.data['COT_K_f_full'],
                              self.data['COT_K_f_full'] / self.data['COT_K_Na'],
                              self.data['COT_K_f_full'] / self.data['COT_K_Cl'],
                              self.data['COT_K_b_full'] / (self.data['COT_K_Cl'] * self.data['COT_K_Na)']),
                              self.data['COT_K_f_full'] / (self.data['COT_K_Cl'] * self.data['COT_K_K)']),
                              self.data['COT_K_b_full'] / (self.data['COT_K_Cl'] ** 2 * self.data['COT_K_K)']),
                              self.data['COT_K_b_full'] / (self.data['COT_K_Cl'] * self.data['COT_K_K']
                                                           * self.data['COT_K_Na']),
                              (self.data['COT_K_b_full'] + self.data['COT_K_f_full']) /
                              (self.data['COT_K_Cl'] ** 2 * self.data['COT_K_K']
                               * self.data['COT_K_Na']),
                              self.data['COT_K_Cl'] ** 2 * self.data['COT_K_K'] * self.data['COT_K_Na']
                              * (self.data['COT_K_b_empty'] + self.data['COT_K_f_empty']))

        # Voltage (Capacitance)
        self.data['A_CAPACI'] = 3.23  # Capacitance of the apical membrane. ------ muF m^-2 (micro-F / meter squared)
        self.data[
            'B_CAPACI'] = 33.4  # Capacitance of the basolateral membrane. ------ muF m^-2 (micro-F / meter squared)

        # Current uses only Faraday's constant, is shown at FLOW

        # ATP, ADP, AMP, etc. NUCLEOTIDE REGULATION
        self.data['J_atp'] = 0.0011  # Nucleotide release. ------ nmol min^-1 L^-1 (nano-moles, minutes, liters)
        self.data['J_adp'] = 0.0131  # Nucleotide release. ------ nmol min^-1 L^-1
        self.data['J_amp'] = 0.0125  # Nucleotide release. ------ nmol min^-1 L^-1
        self.data['V_1_max'] = 6.3  # ------ nmol min^-1 L^-1
        self.data['K_1_m'] = 16.3  # ------ muM (micro-moles)
        self.data['V_2_max'] = 15.5  # ------ nmol min^-1 L^-1
        self.data['K_2_m'] = 114.9  # ------ muM
        self.data['V_3_max'] = 20  # ------ nmol min^-1 L^-1
        self.data['K_3_m'] = 418  # ------ muM
        self.data['V_4_max'] = 0.5  # ------ nmol min^-1 L^-1
        self.data['K_4_m'] = 2.8  # ------ muM
        self.data['V_5_max'] = 10.7  # ------ nmol min^-1 L^-1
        self.data['K_5_m'] = 83.9  # ------ muM
        self.data['V_6_max'] = 1.7  # ------ nmol min^-1 L^-1
        self.data['K_6_m'] = 13  # ------ muM
        self.data['V_7_max'] = 6.2  # ------ nmol min^-1 L^-1
        self.data['K_7_m'] = 27.2  # ------ muM
        self.data['V_8_max'] = 11.9  # ------ nmol min^-1 L^-1
        self.data['K_8_m'] = 694.9  # ------ muM
        self.data['V_9_max'] = 0.3  # ------ nmol min^-1 L^-1
        self.data['K_9_m'] = 17  # ------ muM
        self.data['V_10_max'] = 1.2  # ------ nmol min^-1 L^-1
        self.data['K_10_m'] = 38.2  # ------ muM
        self.data['V_F_max'] = 2.2  # ------ nmol min^-1 L^-1
        self.data['K_F_atp'] = 30.4  # ------ muM
        self.data['K_F_amp'] = 24.7  # ------ muM
        self.data['V_B_max'] = 2.2  # ------ nmol min^-1 L^-1
        self.data['K_B_adp'] = 61.8  # ------ muM
        self.data['V_U1_max'] = 0.2  # ------ nmol min^-1 L^-1
        self.data['K_U1_m'] = 1.2  # ------ muM
        self.data['V_U2_max'] = 0.2  # ------ nmol min^-1 L^-1
        self.data['K_U2_m'] = 1.2  # ------ muM
        self.data['K_IN_atp'] = 28.4  # ------ muM
        self.data['K_IN_adp'] = 20.4  # ------ muM

        if cf:
            # Permeability Constants
            self.data['P_PERM_NA'] = 5e-8  # Permeability of sodium ions (Na+) through the paracellular. ------ m/s
            self.data['P_PERM_K'] = 7.2e-10  # Permeability of potassium ions (K+) through the P. ------ m/s
            self.data['P_PERM_CL'] = 1.2e-8  # Permeability chloride ions (Cl-) through the paracellular. ------ m/s

            # Osmomolarity Constants
            self.data['A_CONS_NCL'] = 48.7  # Apical concentration on non-Chloride anions ------ mM
            self.data[
                'A_CONS_IO'] = 2.7  # Concentration of impermeable osmolytes in the apical compartment. ------ mM (milli-moles / litter)
            self.data['C_CONS_NCL'] = 69.1  # Cellular concentration on non-Chloride anions ------ mM
            self.data[
                'C_CONS_IO'] = 30.4  # Concentration of impermeable osmolytes in the cellular compartment. ------ mM (milli-moles / litter)

            # Flow
            self.data['J_Pump_max'] = 4.4e-6  # Max flow of the Na-K-ATPase pump. ------ milli-mol m^-2 s^-1

            # Permeability EQUATIONS
            self.data['CACC_PERM_MAX'] = 6.9e-9  # The estimated MAX permeability of the CaCC Channel. ------ m / s
            self.data['CACC_ATP_AT_HALF_PERM'] = 3.6e-5  # [ATP] when half receptors occupied (half perm). ------ muM

            self.data['BK_PERM_MAX'] = 4.4e-9  # Similar concept, with the same units.
            self.data['BK_ATP_AT_HALF_PERM'] = 1.5e-3  # Similar concept, with the same units.

            self.data['ENAC_PERM_MAX'] = 2.3e-8  # Similar concept, with the same units.
            self.data['ENAC_ATP_AT_HALF_PERM'] = 5.1e-2  # Similar concept, with the same units.

            self.data['CAKC_PERM_MAX'] = 2.3e-7  # Similar concept, with the same units.
            self.data['CAKC_ATP_AT_HALF_PERM'] = 3.3e-3  # Similar concept, with the same units.
        else:
            # Permeability Constants
            self.data['P_PERM_NA'] = 5e-8  # Permeability of sodium ions (Na+) through the paracellular. ------ m/s
            self.data['P_PERM_K'] = 7.2e-10  # Permeability of potassium ions (K+) through the P. ------ m/s
            self.data['P_PERM_CL'] = 1.2e-8  # Permeability chloride ions (Cl-) through the paracellular. ------ m/s

            # Osmomolarity Constants
            self.data['A_CONS_NCL'] = 48.7  # Apical concentration on non-Chloride anions ------ mM
            self.data[
                'A_CONS_IO'] = 2.7  # Concentration of impermeable osmolytes in the apical compartment. ------ mM (milli-moles / litter)
            self.data['C_CONS_NCL'] = 69.1  # Cellular concentration on non-Chloride anions ------ mM
            self.data[
                'C_CONS_IO'] = 30.4  # Concentration of impermeable osmolytes in the cellular compartment. ------ mM (milli-moles / litter)

            # Flow
            self.data['J_Pump_max'] = 4.4e-6  # Max flow of the Na-K-ATPase pump. ------ milli-mol m^-2 s^-1

            # Permeability EQUATIONS
            self.data['CACC_PERM_MAX'] = 6.9e-9  # The estimated MAX permeability of the CaCC Channel. ------ m / s
            self.data['CACC_ATP_AT_HALF_PERM'] = 3.6e-5  # [ATP] when half receptors occupied (half perm). ------ muM

            self.data['BK_PERM_MAX'] = 4.4e-9  # Similar concept, with the same units.
            self.data['BK_ATP_AT_HALF_PERM'] = 1.5e-3  # Similar concept, with the same units.

            self.data['ENAC_PERM_MAX'] = 2.3e-8  # Similar concept, with the same units.
            self.data['ENAC_ATP_AT_HALF_PERM'] = 5.1e-2  # Similar concept, with the same units.

            self.data['CAKC_PERM_MAX'] = 2.3e-7  # Similar concept, with the same units.
            self.data['CAKC_ATP_AT_HALF_PERM'] = 3.3e-3  # Similar concept, with the same units.

        if constants is not None:
            for val, key in enumerate(constants):
                if key in self.data:
                    self.data[key] = val
                else:
                    pass

    def __getitem__(self, item):
        return self.data[item]


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
                * a pandas.DataFrame:
                    It will be used to load an already full DataFrame.
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
                * a pandas.DataFrame:
                    It will be used to load an already full DataFrame.
                * path to a .csv:
                    It will be used to load the DataFrame from the .csv with the pandas.read_csv() function.
                * DEFAULT: None (it will ask the user for values).
        :exception if there's an issue loading the values, they will be asked manually.
        """

        # Will create the needed dictionary to init the variables
        if isinstance(init_data, str) and '.csv' in init_data:
            self.init_pd = pd.read_csv(init_data)
            init_data = {}

            if all(k in self.init_pd for k in self.initial_vars):
                for var in self.initial_vars:
                    if var == 'max_steps':
                        init_data[var] = int(self.init_pd[var][0])
                    elif var == 'CF':
                        init_data[var] = bool(self.init_pd[var][0])
                    else:
                        init_data[var] = float(self.init_pd[var][0])
            
            for var in self.voltage:
                if var in self.init_pd:
                    init_data[var] = float(self.init_pd[var])

            del(self.init_pd)

        # If the right argument was passed or it was rightly processed, start the init_vars
        if isinstance(init_data, dict) and all(k in init_data for k in self.initial_vars):
            self.max_steps = int(init_data['max_steps'])
            # Here we are creating the data-frame from the list of 'self.variables'
            self.data = pd.DataFrame(np.zeros((self.max_steps, len(self.variables))), columns=self.variables)

            self.time_frame = float(init_data['time_frame'])
            self.isCF = bool(init_data['CF'])

            self.co = Constants(self.isCF, special_constants)

            for ion in self.apical_ions + self.cel_ions:
                self.data[ion][0] = self.to_activity(np.absolute(float(init_data[ion])))

            if all(v in init_data for v in ('aV', 'bV')):
                self.data['aV'][0] = float(init_data['aV'])
                self.data['bV'][0] = float(init_data['bV'])

                if 'pV' in init_data:
                    self.data['pV'][0] = float(init_data['pV'])
                else:
                    self.data["tV"][0] = self.data["bV"][0] - self.data["aV"][0]
        else:
            self.inputVals()
            self.co = Constants(self.isCF, special_constants)

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

            if not all(v in init_data for v in ('aV', 'bV')):
                self.data["aV"][0] = self.fn_voltage('a')
                self.data["bV"][0] = self.fn_voltage('b')
                self.data["tV"][0] = self.data["bV"][0] - self.data["aV"][0]

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

    def inputVals(self):
        """
        This function is used to ask for all of the necessary initial values. It will go one by one
        of the variables in AirwayModel.initial_vars asking the user to manually input then and then
        load it to the corresponding variable. 
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

        self.time_frame = float(input('Enter how many seconds each step represents: '))

        self.data['H'][0] = float(input('Enter the ASL height in meters: '))

        for ion in self.initial_vars[3:9]:
            self.data[ion][0] = np.absolute(float(input('Enter the value for ' + ion + ' (in mM): ')))
        for nucl in self.initial_vars[9:-1]:
            self.data[nucl][0] = np.absolute(float(input('Enter the value for ' + nucl + ' (in micro-M): ')))

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
                                 / self.data['H'][step - 1]  # in mM
        self.data["aCl"][step] = self.data["aCl"][step - 1] + self.fn_daN_Cl(step) * self.time_frame \
                                 / self.data['H'][step - 1]  # in mM
        self.data["aK"][step] = self.data["aK"][step - 1] + self.fn_daN_K(step) * self.time_frame \
                                / self.data['H'][step - 1]  # in mM

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
               * (self.data['cCl'][step - 1] - self.co['B_ACT_CL '] * np.exp(
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
                + self['V_8_max'] * self.data["AMP"][step - 1] / (
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
