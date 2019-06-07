from multiprocessing import Pool
import constants_NL as co_NL
import constants_CF as co_CF
import pandas as pd
import numpy as np
import datetime
import zipfile
import sys
import os

# We need to stop the simulations after any error, else data will be harder to analyze.
prev_np_err = np.seterr(all='raise')
# Will be the folder where all of the temporary files will be saved.
SUB_DIR = './data'
# Will be the MAX files that can be on the SUB_DIR before they are sent to a ZIP file.
MAX_FILES = 5
# Max number of parallel processes.
MAX_PROCESSES = 4


# A class with all the necessary functions to run an ASL Model.
class AirwayModel:
    """
    Basically a Pandas.DataFrame with all of the necessary functions to
    run the simulation.

    It's recommended to use with a for loop iterating through step in range(init_data['max_steps'])
    where you run AirwayModel.fn_run(step). This is because the model which this class is meant to be used on
    updates the values based on the previous one.

    __init__(self, init_data):
        :arg init_data, a dictionary with:
            :var 'max_steps': How long with the simulation be.
            :var 'H': initial water height (meters)
            :var Molar concentration for all the relevant ions (mM)
            :var Concentrations for all nucleotides (ATP, ..., ADO, INO) in muM (micro-moles)
        :exception if there's an issue loading the CSV file, the program will exit with code 3.
    fn_run(self, step):
        :arg step=1: is the step to load the values to, by default will load the values after the initial.
        :returns nothing
    """
    # Complete with all variables that will be tracked as column names
    variables = ['Time (min)', 'H', 'dH', 'OSM_a', 'OSM_c', 'aNa', 'aCl', 'aK', 'cNa', 'cCl', 'cK',
                 'aJ_Na', 'pJ_Na', 'aJ_K', 'pJ_K', 'aJ_Cl', 'pJ_Cl', 'J_co', 'J_pump',
                 'bJ_Cl', 'bJ_K', 'daV', 'dbV', 'aV', 'bV', 'tV', 'aI', 'bI', 'pI',
                 'p_CaCC', 'p_CFTR', 'p_ENaC', 'p_BK', 'p_CaKC',
                 'p_CaCC_rate', 'p_CFTR_rate', 'p_ENaC_rate', 'p_BK_rate', 'p_CaKC_rate',
                 'ATP', 'dATP', 'ADP', 'dADP', 'AMP', 'dAMP', 'ADO', 'dADO', 'INO', 'dINO']

    # Separated variables into sections
    apical_ions = ['aNa', 'aCl', 'aK']
    cel_ions = ['cNa', 'cCl', 'cK']
    apical_flow = ['aJ_Na', 'aJ_K', 'aJ_Cl']
    paracellular_flow = ['pJ_Na', 'pJ_K', 'pJ_Cl']
    basolateral_flow = ['J_co', 'J_pump', 'bJ_Cl', 'bJ_K']
    voltage = ['aV', 'bV', 'tV']
    current = ['aI', 'bI', 'pI']
    permeabilities = ['p_CaCC', 'p_CFTR', 'p_ENaC', 'p_BK', 'p_CaKC']
    permeabilities_rate = ['p_CaCC_rate', 'p_CFTR_rate', 'p_ENaC_rate', 'p_BK_rate', 'p_CaKC_rate']
    nucleotide = ['ATP', 'ADP', 'AMP', 'ADO', 'INO']
    nucleotide_dt = ['dATP', 'dADP', 'dAMP', 'dADO', 'dINO']

    sections = ['H', 'dH', apical_ions, cel_ions, apical_flow, paracellular_flow, basolateral_flow,
                voltage, current, permeabilities, permeabilities_rate, nucleotide, nucleotide_dt]

    # Required initial values
    initial_vars = ['max_steps', 'time_frame', 'H', 'aNa', 'aCl', 'aK', 'cNa', 'cCl', 'cK',
                    'ATP', 'ADP', 'AMP', 'ADO', 'INO', 'CF']

    # Initial conditions setup
    def __init__(self, init_data=None):
        """
        Will initialize the tracked values. If not successfully from init_data, it will prompt the user
        add them manually.

        :param init_data: is a dictionary with the initial values.
        """
        if init_data is None:
            self.inputVals()
        else:
            try:
                self.max_steps = int(init_data['max_steps'])
                # Here we are creating the data-frame from the list of 'self.variables'
                self.data = pd.DataFrame(np.zeros((self.max_steps, len(self.variables))), columns=self.variables)

                self.time_frame = float(init_data['time_frame'])
                self.isCF = bool(init_data['CF'])

                for ion in self.initial_vars[2:-1]:
                    self.data[ion][0] = np.absolute(float(init_data[ion]))
            except KeyError or TypeError or ValueError:
                # If there's any error ask the user to input the values
                self.inputVals()

        # Set the fn_run() and constants to the corresponding mode, CF or NL
        if self.isCF:
            self.fn_run = self.fn_run_CF
            self.co = co_CF
        else:
            self.fn_run = self.fn_run_NL
            self.co = co_NL

        self.data['Time (min)'][0] = 0
        self.data["OSM_a"][0] = self.co.PHI * (
                self.data["aNa"][0] + self.data["aK"][0] + self.data["aCl"][0] + self.co.A_ACT_OI) / self.co.GAMMA
        self.data["OSM_c"][0] = self.co.PHI * (
                self.data["cNa"][0] + self.data["cK"][0] + self.data["cCl"][0] + self.co.C_ACT_OI) / self.co.GAMMA
        self.data["dH"][0] = self.fn_dH(1)

        self.data["aV"][0] = self.fn_voltage('a')
        self.data["bV"][0] = self.fn_voltage('b')
        self.data["tV"][0] = self.data["bV"][0] - self.data["aV"][0]

        self.data["p_CaCC"][0] = self.fn_p_CaCC(1)
        if not self.isCF:
            self.data["p_CFTR"][0] = self.fn_p_CFTR(1)
            self.data["p_ENaC"][0] = self.fn_p_ENaC(1)
        else:
            self.data["p_ENaC"][0] = self.fn_p_ENaC_CF(1)
        self.data["p_BK"][0] = self.fn_p_BK(1)
        self.data["p_CaKC"][0] = self.fn_p_CaKC(1)

        self.data["p_CaCC_rate"][0] = self.co.CACC_PERM_MAX / self.data["p_CaCC"][0]
        if not self.isCF:
            self.data["p_CFTR_rate"][0] = self.co.CFTR_PERM_MAX / self.data["p_CFTR"][0]
        self.data["p_ENaC_rate"][0] = self.co.ENAC_PERM_MAX / self.data["p_ENaC"][0]
        self.data["p_BK_rate"][0] = self.co.BK_PERM_MAX / self.data["p_BK"][0]
        self.data["p_CaKC_rate"][0] = self.co.CAKC_PERM_MAX / self.data["p_CaKC"][0]

        self.data["aJ_Na"][0] = self.fn_aJ_Na(1)
        self.data["pJ_Na"][0] = self.fn_pJ_Na(1)
        self.data["aJ_K"][0] = self.fn_aJ_K(1)
        self.data["pJ_K"][0] = self.fn_pJ_K(1)
        self.data["bJ_K"][0] = self.fn_bJ_K(1)
        self.data["aJ_Cl"][0] = self.fn_aJ_Cl(1)
        self.data["pJ_Cl"][0] = self.fn_pJ_Cl(1)
        self.data["bJ_Cl"][0] = self.fn_bJ_Cl(1)
        self.data["J_co"][0] = self.fn_J_co(1)
        self.data["J_pump"][0] = self.fn_J_pump(1)

        self.data["aI"][0] = self.co.F_RT * (- self.data["aJ_Na"][0] + self.data["aJ_Cl"][0] - self.data["aJ_K"][0])
        self.data["bI"][0] = self.co.F_RT * (- self.data["J_pump"][0] + self.data["bJ_Cl"][0] - self.data["aJ_K"][0])
        self.data["pI"][0] = self.co.F_RT * (- self.data["pJ_Na"][0] + self.data["pJ_Cl"][0] - self.data["pJ_K"][0])
        self.data["daV"][0] = (self.data["pI"][0] - self.data["aI"][0]) / (self.co.A_CAPACI * 100)
        self.data["dbV"][0] = - (self.data["pI"][0] + self.data["bI"][0]) / (self.co.B_CAPACI * 100)

        self.data["dATP"][0] = self.fn_dATP(1) / 60  # Remember that it's in min
        self.data["dADP"][0] = self.fn_dADP(1) / 60  # Remember that it's in min
        self.data["dAMP"][0] = self.fn_dAMP(1) / 60  # Remember that it's in min
        self.data["dADO"][0] = self.fn_dADO(1) / 60  # Remember that it's in min
        self.data["dINO"][0] = self.fn_dINO(1) / 60  # Remember that it's in min

    def __getitem__(self, item):
        return self.data[item]

    def __repr__(self):
        """Since the main object is a the pd.DataFrame, it will be what's printed"""
        return 'ASL Model\n' + str(self.data)

    def __str__(self):
        """When you're working with the model, the important thing is the data, it will print the
        self.data object which is a pd.DataFrame"""
        return str(self.data)

    # Again, to make the pd.DataFrame more accessible.
    def drop(self, labels=None, axis=0, index=None,
             columns=None, level=None, inplace=False, errors='raise'):
        return self.data.drop(labels=labels, axis=axis, index=index, columns=columns, level=level, inplace=inplace,
                              errors=errors)

    def to_cons(self, activity):
        """Return: activity / self.co.GAMMA = concentration"""
        return activity / self.co.GAMMA

    def to_activity(self, concentration):
        """Return: activity = concentration * self.co.GAMMA"""
        return concentration * self.co.GAMMA

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
                     'Any input error will count as Normal Conditions: ')
        if '1' in tmp or ('C' in tmp and 'F' in tmp):
            self.isCF = True
            print('\tCystic Fibrosis activated.')
        else:
            self.isCF = False
            print('\tNormal Conditions activated.')

        self.time_frame = float(input('Enter how many seconds each step represents: '))

        self.data['H'][0] = int(input('Enter the ASL height in meters: '))

        for ion in self.initial_vars[2:8]:
            self.data[ion][0] = np.absolute(float(input('Enter the value for ' + ion + ' (in mM): ')))
        for nucl in self.initial_vars[8:-1]:
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

    # It's here to keep the old method of working going on. data.to_csv('fileName.csv')
    def to_csv(self,  name, sub_dir='./'):
        """
        Will save self.data into a .CSV file. Defaults to the current directory.
        If the directory is not found, it will use os.mkdir() to create the desired directory.

        :param name: Is the name of the file, if the '.csv' is not in the name, it will be added.
                NOTE: It should not contain the path of the file, only its name.
        :param sub_dir: Is the directory to which to save the file.
        :return:
        """
        if sub_dir[-1] != '/' and sub_dir[-1] != '\\':
            sub_dir += '/'
        if '.csv' not in name:
            name += '.csv'
        try:
            self.data.to_csv(sub_dir + name)
        except FileNotFoundError:
            os.mkdir(sub_dir)
            self.data.to_csv(sub_dir + name)

    # Voltage Calculator, approximate!
    def fn_voltage(self, membrane='a', step=1):
        """
        Uses the Ernst equation to estimate the voltage of 'membrane' at 'step'.
        Uses the main ions in the compartment 'a' or 'b' depending on the membrane and the ones inside the cell.
        Main ions are: Chloride, Sodium, and Potassium.
        :param membrane: The membrane to which to calculate it's voltage, if the transmembrane potential is asked
                ('t' or 'p') it will do so by calculating the sum of the apical and basolateral potentials.
        :param step: The step to which to calculate the potential to, it will again use the values from step-1 to
                calculate it's results.
        :return: It will return the estimated potential in Volts.
        """
        if membrane == 'a':
            v = 0
            for ion in ['Na', 'K']:
                v += np.log(self.data[membrane + ion][step-1] / self.data['c' + ion][step-1]) / self.co.F_RT
            for ion in ['Cl']:
                v -= np.log(self.data[membrane + ion][step-1] / self.data['c' + ion][step-1]) / self.co.F_RT
            return v
        elif membrane == 'b':
            return np.log(self.co.B_CONS_NA / self.data['cNa'][step-1]) / self.co.F_RT \
                   + np.log(self.co.B_CONS_K / self.data['cK'][step-1]) / self.co.F_RT \
                   - np.log(self.co.B_CONS_CL / self.data['cCl'][step-1]) / self.co.F_RT
        elif membrane == 'p' or membrane == 't':
            return self.fn_voltage('b', step) - self.fn_voltage('a', step)
        else:
            raise KeyError('Must choose either apical, basolateral or transmembrane (paracellular).')

    # Run Eq in Normal Conditions
    def fn_run_NL(self, step=1):
        """
        Will run all the required Eqs and update the values.
        :param step: Is the step to which to calculate the values to. When required to obtain a value, it will
                collect it from step-1.
                Defaults to 1
        """

        self.data['Time (min)'][step] = step * self.time_frame / 60
        self.data["dH"][step] = self.fn_dH(step)
        self.data["H"][step] = self.data["H"][step - 1] + self.data["dH"][step] * self.time_frame
        self.data["OSM_a"][step] = self.fn_OSM_a(step)
        self.data["OSM_c"][step] = self.fn_OSM_c(step)

        self.data["aNa"][step] = self.data["aNa"][step - 1] + self.fn_daN_Na(step) * self.time_frame / \
                                 self.data['H'][
                                     step - 1]  # in mM
        self.data["aCl"][step] = self.data["aCl"][step - 1] + self.fn_daN_Cl(step) * self.time_frame / \
                                 self.data['H'][
                                     step - 1]  # in mM
        self.data["aK"][step] = self.data["aK"][step - 1] + self.fn_daN_K(step) * self.time_frame / self.data['H'][
            step - 1]  # in mM

        # It's now in mol / m^2
        self.data["cNa"][step] = self.data["cNa"][step - 1] + self.fn_dcN_Na(step) * self.time_frame
        self.data["cCl"][step] = self.data["cCl"][step - 1] + self.fn_dcN_Cl(step) * self.time_frame
        self.data["cK"][step] = self.data["cK"][step - 1] + self.fn_dcN_K(step) * self.time_frame

        self.data["aJ_Na"][step] = self.fn_aJ_Na(step)
        self.data["aJ_Cl"][step] = self.fn_aJ_Cl(step)
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

        # Will be to know how far they are from max open (closer to '1')
        self.data["p_CaCC_rate"][step] = self.co.CACC_PERM_MAX / self.data["p_CaCC"][step]
        self.data["p_CFTR_rate"][step] = self.co.CFTR_PERM_MAX / self.data["p_CFTR"][step]
        self.data["p_ENaC_rate"][step] = self.co.ENAC_PERM_MAX / self.data["p_ENaC"][step]
        self.data["p_BK_rate"][step] = self.co.BK_PERM_MAX / self.data["p_BK"][step]
        self.data["p_CaKC_rate"][step] = self.co.CAKC_PERM_MAX / self.data["p_CaKC"][step]

        self.data["dATP"][step] = self.fn_dATP(step) / 60  # To get it to muM/sec instead of muM/min
        self.data["ATP"][step] = self.data["ATP"][step - 1] + self.data["dATP"][step] * self.time_frame
        self.data["dADP"][step] = self.fn_dADP(step) / 60  # To get it to muM/sec instead of muM/min
        self.data["ADP"][step] = self.data["ADP"][step] + self.data["dADP"][step] * self.time_frame
        self.data["dAMP"][step] = self.fn_dAMP(step) / 60  # To get it to muM/sec instead of muM/min
        self.data["AMP"][step] = self.data["AMP"][step - 1] + self.data["dAMP"][step] * self.time_frame
        self.data["dADO"][step] = self.fn_dADO(step) / 60  # To get it to muM/sec instead of muM/min
        self.data["ADO"][step] = self.data["ADO"][step - 1] + self.data["dADO"][step] * self.time_frame
        self.data["dINO"][step] = self.fn_dINO(step) / 60  # To get it to muM/sec instead of muM/min
        self.data["INO"][step] = self.data["INO"][step - 1] + self.data["dINO"][step] * self.time_frame

    # Run Eq in Cystic Fibrosis
    def fn_run_CF(self, step=1):
        """
        Cystic Fibrosis.
        Will run all the required Eqs and update the values.
        :param step: Is the step to which to calculate the values to. When required to obtain a value, it will
                collect it from step-1.
                Defaults to 1
        """

        self.data['Time (min)'][step] = step * self.time_frame / 60
        self.data["dH"][step] = self.fn_dH(step)
        self.data["H"][step] = self.data["H"][step - 1] + self.data["dH"][step] * self.time_frame
        self.data["OSM_a"][step] = self.fn_OSM_a(step)
        self.data["OSM_c"][step] = self.fn_OSM_c(step)

        self.data["aNa"][step] = self.data["aNa"][step - 1] + self.fn_daN_Na(step) * self.time_frame / \
                                 self.data['H'][
                                     step - 1]  # in mM
        self.data["aCl"][step] = self.data["aCl"][step - 1] + self.fn_daN_Cl(step) * self.time_frame / \
                                 self.data['H'][
                                     step - 1]  # in mM
        self.data["aK"][step] = self.data["aK"][step - 1] + self.fn_daN_K(step) * self.time_frame / self.data['H'][
            step - 1]  # in mM

        # It's now in mol / m^2
        self.data["cNa"][step] = self.data["cNa"][step - 1] + self.fn_dcN_Na(step) * self.time_frame
        self.data["cCl"][step] = self.data["cCl"][step - 1] + self.fn_dcN_Cl(step) * self.time_frame
        self.data["cK"][step] = self.data["cK"][step - 1] + self.fn_dcN_K(step) * self.time_frame

        self.data["aJ_Na"][step] = self.fn_aJ_Na(step)
        self.data["aJ_Cl"][step] = self.fn_aJ_Cl(step)
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
        self.data["p_ENaC"][step] = self.fn_p_ENaC_CF(step)
        self.data["p_BK"][step] = self.fn_p_BK(step)
        self.data["p_CaKC"][step] = self.fn_p_CaKC(step)

        # Will be to know how far they are from max open (closer to '1')
        self.data["p_CaCC_rate"][step] = self.co.CACC_PERM_MAX / self.data["p_CaCC"][step]
        self.data["p_ENaC_rate"][step] = self.co.ENAC_PERM_MAX / self.data["p_ENaC"][step]
        self.data["p_BK_rate"][step] = self.co.BK_PERM_MAX / self.data["p_BK"][step]
        self.data["p_CaKC_rate"][step] = self.co.CAKC_PERM_MAX / self.data["p_CaKC"][step]

        self.data["dATP"][step] = self.fn_dATP(step) / 60  # To get it to muM/sec instead of muM/min
        self.data["ATP"][step] = self.data["ATP"][step - 1] + self.data["dATP"][step] * self.time_frame
        self.data["dADP"][step] = self.fn_dADP(step) / 60  # To get it to muM/sec instead of muM/min
        self.data["ADP"][step] = self.data["ADP"][step] + self.data["dADP"][step] * self.time_frame
        self.data["dAMP"][step] = self.fn_dAMP(step) / 60  # To get it to muM/sec instead of muM/min
        self.data["AMP"][step] = self.data["AMP"][step - 1] + self.data["dAMP"][step] * self.time_frame
        self.data["dADO"][step] = self.fn_dADO(step) / 60  # To get it to muM/sec instead of muM/min
        self.data["ADO"][step] = self.data["ADO"][step - 1] + self.data["dADO"][step] * self.time_frame
        self.data["dINO"][step] = self.fn_dINO(step) / 60  # To get it to muM/sec instead of muM/min
        self.data["INO"][step] = self.data["INO"][step - 1] + self.data["dINO"][step] * self.time_frame

    # Eq 1
    def fn_dH(self, step=1):
        """Calculates the rate of change of the apical water height in meters per second (m/s)"""
        return self.co.A_V_H2O * (self.co.P_PERM_H2O * (self.data["OSM_a"][step-1] - self.co.B_OSM) +
                                  self.co.A_PERM_H20 * (self.data["OSM_a"][step-1] - self.data["OSM_c"][step-1]))

    # Eq 2.1.2
    def fn_OSM_a(self, step=1):
        """Will calculate the osmomolarity of the Airway Liquid at step in mili-molar or moles per meter cubed"""
        return self.co.PHI / self.co.GAMMA * (self.data["aNa"][step-1] + self.data["aK"][step-1]
                                              + self.data["aCl"][step-1] + self.co.A_ACT_OI)

    # Eq 2.2.2
    def fn_OSM_c(self, step=1):
        """Will calculate the osmomolarity of the cell at step in mili-molar or moles per meter cubed"""
        return self.co.PHI / self.co.GAMMA * (
                self.data["cNa"][step-1] + self.data["cK"][step-1] + self.data["cCl"][step-1] + self.co.C_ACT_OI)

    # Eq 3.1.1
    def fn_daN_Na(self, step=1):
        """Difference of flow of Sodium across the apical membrane"""
        return -self.data['aJ_Na'][step-1] + self.data['pJ_Na'][step-1]

    # Eq 3.1.2
    def fn_daN_K(self, step=1):
        """Difference of flow of Potassium across the apical membrane"""
        return -self.data['aJ_K'][step-1] + self.data['pJ_K'][step-1]

    # Eq 3.1.3
    def fn_daN_Cl(self, step=1):
        """Difference of flow of Chloride across the apical membrane"""
        return -self.data['aJ_Cl'][step-1] + self.data['pJ_Cl'][step-1]

    # Eq 3.2.1
    def fn_dcN_Na(self, step=1):
        """Difference of flow of Sodium across the basolateral membrane"""
        return self.data['aJ_Na'][step-1] + self.data['J_co'][step-1] - 3 * self.data['J_pump'][step-1]

    # Eq 3.2.2
    def fn_dcN_K(self, step=1):
        """Difference of flow of Potassium across the basolateral membrane"""
        return self.data['aJ_K'][step-1] + self.data['bJ_K'][step-1] + self.data['J_co'][step-1] + 2 * \
               self.data['J_pump'][step-1]

    # Eq 3.2.3
    def fn_dcN_Cl(self, step=1):
        """Difference of flow of Chloride across the basolateral membrane"""
        return self.data['aJ_Cl'][step-1] + self.data['bJ_Cl'][step-1] + 2 * self.data['J_co'][step-1]

    # Eq 4.1
    def fn_aJ_Na(self, step=1):
        """
        Will return the flow of Sodium ions across the apical membrane.

        :return: Flow moles per meter squared per second, mol / (m^2 sec)
        """
        return self.data["p_ENaC"][step-1] * self.data["aV"][step-1] * self.co.F_RT \
               * (self.data["aNa"][step-1] -
                  self.data["cNa"][step-1] * np.exp(self.data["aV"][step-1] * self.co.F_RT)) \
               / (np.exp(self.data["aV"][step-1] * self.co.F_RT) - 1)

    # Eq 4.2
    def fn_aJ_K(self, step=1):
        """
        :return: Flow of potassium across the apical membrane at a given step (moles / (m^2 sec)).
        """
        return self.data["p_BK"][step-1] * self.data["aV"][step-1] * self.co.F_RT / (
                np.exp(self.data["aV"][step-1] * self.co.F_RT) - 1) \
               * (self.data["aK"][step-1] - self.data["cK"][step-1] * np.exp(self.co.F_RT * self.data["aV"][step-1]))

    # Eq 4.3
    def fn_aJ_Cl(self, step=1):
        """
        :return: Flow of chloride across the apical membrane at a given step (moles / (m^2 sec)).
        """
        return -(self.data["p_CaCC"][step-1] + self.data["p_CFTR"][step-1]) \
               * self.co.F_RT * self.data["aV"][step-1] / (np.exp(self.co.F_RT * self.data["aV"][step-1]) - 1) \
               * (self.data["cCl"][step-1] - self.data["aCl"][step-1] * np.exp(self.co.F_RT * self.data["aV"][step-1]))

    # Eq 4.4
    def fn_aJ_H2O(self, step=1):
        """moles per second per meter squared - moles / (sec m^2)"""
        return self.co.A_PERM_H20 * (self.data["OSM_c"][step-1] - self.data["OSM-a"][step-1])

    # Eq 4.5
    def fn_bJ_Cl(self, step=1):
        """moles per second per meter squared - moles / (sec m^2)"""
        return - self.co.P_PERM_CL * self.co.F_RT * self.data["bV"][step-1] * (np.exp(self.co.F_RT *
                                                                                      self.data["bV"][step-1]) - 1) \
               * (self.data['cCl'][step-1] - self.co.B_ACT_CL * np.exp(self.co.F_RT * self.data["bV"][step-1]))

    # Eq 4.6
    def fn_bJ_K(self, step=1):
        """moles per second per meter squared - moles / (sec m^2)"""
        return - self.co.P_PERM_K * self.co.F_RT * self.data["bV"][step-1] / (np.exp(self.co.F_RT
                                                                                     * self.data["bV"][step-1]) - 1) \
               * (self.co.B_ACT_K - self.data["cK"][step-1] * np.exp(self.co.F_RT * self.data["bV"][step-1]))

    # Eq 4.7
    def fn_J_pump(self, step=1):
        """moles per second per meter squared - moles / (sec m^2)"""
        return self.co.J_Pump_max * (self.data['cNa'][step-1] / (self.data['cNa'][step-1]
                                                                 + self.co.K_Na_In_pump
                                                                 * (1 + self.data['cK'][step-1]
                                                                        / self.co.K_K_in_pump))) ** 3 \
               * self.co.B_ACT_K / (self.co.B_ACT_K +
                                    self.co.K_K_ext_pump * (1 + self.co.B_ACT_NA / self.co.K_Na_ext_pump)) ** 2

    # Eq 4.8: System
    def fn_J_co_sys(self, step=1):
        """
        Will return a list with all the values of E.
        The first array is the 10x10 matrix for the cotransporter equation.
        It uses Numpy to solve the homogeneous linear system.
        """

        return np.linalg.solve(np.array(
            [[-self.co.COT_K_ON * self.co.B_CONS_NA - self.co.COT_K_BE, self.co.COT_K_1, 0,
              0, 0, 0, 0, 0, 0, self.co.COT_K_FE],
             [self.co.COT_K_ON * self.co.B_CONS_NA, -self.co.COT_K_ON * self.co.B_CONS_CL - self.co.COT_K_1,
              self.co.COT_K_2, 0, 0, 0, 0, 0, 0, 0],
             [0, self.co.COT_K_ON * self.co.B_CONS_CL, -self.co.COT_K_ON * self.co.B_CONS_K - self.co.COT_K_2,
              self.co.COT_K_03, 0, 0, 0, 0, 0, 0],
             [0, 0, self.co.COT_K_ON * self.co.B_CONS_K, -self.co.COT_K_ON * self.co.B_CONS_CL - self.co.COT_K_03,
              self.co.COT_K_4, 0, 0, 0, 0, 0],
             [0, 0, 0, self.co.COT_K_ON * self.co.B_CONS_CL, -self.co.COT_K_4 - self.co.COT_K_FF,
              self.co.COT_K_BF, 0, 0, 0, 0],
             [0, 0, 0, 0, self.co.COT_K_FF, -self.co.COT_K_BF - self.co.COT_K_1,
              self.co.COT_K_ON * self.to_cons(self.data["cNa"][step-1]), 0, 0, 0],
             [0, 0, 0, 0, 0, self.co.COT_K_1,
              -self.co.COT_K_ON * self.to_cons(self.data["cNa"][step-1]) - self.co.COT_K_2,
              self.co.COT_K_ON * self.to_cons(self.data["cCl"][step-1]), 0, 0],
             [0, 0, 0, 0, 0, 0, self.co.COT_K_2,
              -self.co.COT_K_ON * self.to_cons(self.data["cCl"][step-1]) - self.co.COT_K_I3,
              self.co.COT_K_ON * self.to_cons(self.data["cK"][step-1]), 0],
             [0, 0, 0, 0, 0, 0, 0, self.co.COT_K_I3,
              -self.co.COT_K_ON * self.to_cons(self.data["cK"][step-1]) - self.co.COT_K_4,
              self.co.COT_K_ON * self.to_cons(self.data["cCl"][step-1])],
             [self.co.COT_K_BE, 0, 0, 0, 0, 0, 0, 0,
              self.co.COT_K_4, -self.co.COT_K_FE - self.co.COT_K_ON * self.to_cons(self.data["cCl"][step-1])]
             ]), np.array([0] * 10))

    # Eq 4.8
    def fn_J_co(self, step=1):
        """
        Will return the flow of the cotransporter.
        The function fn_J_co_sys() will be called to solve the system.
        """

        e_var = self.fn_J_co_sys(step)
        return self.co.COT_K_FF * e_var[4] - self.co.COT_K_BF * e_var[5]

    # Eq 4.9
    def fn_pJ_H2O(self, step=1):
        """
        Water flow through the paracellular membrane. This equation is no needed because its effects
        are accounted for in AirwayModel.fn_dH()
        """
        return self.co.B_PERM_H20 * (self.data['OSM_c'][step-1] - self.co.B_OSM)

    # Eq 4.10
    def fn_pJ_Na(self, step=1):
        """Paracellular flow of for the Sodium ion, in moles per sec per meters squared - mole / (sec m^2)"""
        return self.co.P_PERM_NA * (self.co.F_RT * self.data["tV"][step-1] /
                                    ((self.co.R * self.co.TEMP) * (np.exp(self.co.F_RT * self.data["tV"][step-1] /
                                                                          (self.co.R * self.co.TEMP)) - 1))) * (
                       self.co.B_CONS_NA - self.data["aNa"][step-1] *
                       np.exp(self.co.F_RT * self.data["tV"][step-1] / (self.co.R * self.co.TEMP)))

    # Eq 4.11
    def fn_pJ_K(self, step=1):
        """Paracellular flow of for the Potassium ion, in moles per sec per meters squared - mole / (sec m^2)"""
        return self.co.P_PERM_K * (self.co.F_RT * self.data["tV"][step-1] /
                                   ((self.co.R * self.co.TEMP) * (np.exp(self.co.F_RT * self.data["tV"][step-1] /
                                                                         (self.co.R * self.co.TEMP)) - 1))) * (
                       self.co.B_CONS_K - self.data["aK"][step-1] *
                       np.exp(self.co.F_RT * self.data["tV"][step-1] / (self.co.R * self.co.TEMP)))

    # Eq 4.12
    def fn_pJ_Cl(self, step=1):
        """Paracellular flow of for the Chloride ion, in moles per sec per meters squared - mole / (sec m^2)"""
        return -self.co.P_PERM_CL * (self.co.F_RT * self.data["tV"][step-1] /
                                     ((self.co.R * self.co.TEMP) * (np.exp(self.co.F_RT * self.data["tV"][step-1] /
                                                                           (self.co.R * self.co.TEMP)) - 1))) * (
                       self.data["aCl"][step-1] - self.co.B_CONS_CL *
                       np.exp(self.co.F_RT * self.data["tV"][step-1] / (self.co.R * self.co.TEMP)))

    # Eq 5.1.1
    def fn_p_CaCC(self, step=1):
        """
        Returns the permeability of the CaCC Channel.
        :return: in meters per second (m / sec)
        """
        return self.co.CACC_PERM_MAX / (1 + self.co.CACC_ATP_AT_HALF_PERM / self.data["ATP"][step - 1])

    # Eq 5.1.2
    def fn_p_CFTR(self, step=1):
        """
        Returns the permeability of the CFTR channel
        :return: meters per second (m / sec)
        """
        return self.co.CFTR_PERM_MAX / (1 + self.co.CFTR_ADO_AT_HALF_PERM / self.data["ADO"][step-1] +
                                        self.co.CFTR_ATP_AT_HALF_PERM / self.data["ATP"][step-1])

    # Eq 5.2
    def fn_p_BK(self, step=1):
        """
        Returns the permeability of the BK channel.
        :returns: in meters per second (m / sec).
        """
        return self.co.BK_PERM_MAX / (1 + self.co.BK_ATP_AT_HALF_PERM / self.data["ATP"][step-1])

    # Eq 5.3 for NL
    def fn_p_ENaC(self, step=1):
        """
        Is the permeability of the ENaC channel (Sodium)
        Now it's different for C.F.
        :return: The permeability for step in m/sec.
        """
        return self.co.ENAC_PERM_MAX / (1 + self.co.ENAC_ATP_AT_HALF_PERM / self.data['ATP'][step-1]
                                        + self.co.ENAC_ADO_AT_HALF_PERM / self.data['ADO'][step-1])

    # Eq 5.3 for CF
    def fn_p_ENaC_CF(self, step=1):
        """
        Is the permeability of the ENaC channel (Sodium)
        Now it's different for C.F.
        :return: The permeability for step in m/sec.
        """

        return self.co.ENAC_PERM_MAX / (1 + self.co.ENAC_ATP_AT_HALF_PERM / self.data['ATP'][step - 1])

    # Eq 5.4
    def fn_p_CaKC(self, step=1):
        """Permeability of the CaKC channel in meters per second (m/sec)"""

        return self.co.CAKC_PERM_MAX / (1 + self.co.CAKC_ATP_AT_HALF_PERM / self.data['ATP'][step-1])

    # Eq 6.1
    def fn_daV(self, step=1):
        """Will return the rate of change of the apical membrane's voltage. In Volts per second"""
        return (self.data["pI"][step-1] - self.data["aI"][step-1]) / (self.co.A_CAPACI * 100)

    # Eq 6.2
    def fn_dbV(self, step=1):
        """Will return the rate of change of the basolateral membrane's voltage. In Volts per second"""
        return - (self.data["pI"][step-1] + self.data["bI"][step-1]) / (self.co.B_CAPACI * 100)

    # Eq 7.1
    def fn_pI(self, step=1):
        """Will return paracellular current in A / m^2 or C / (sec m^2)"""
        return self.co.F_RT * (- self.data["pJ_Na"][step-1] + self.data["pJ_Cl"][step-1] - self.data["pJ_K"][step-1])

    # Eq 7.2
    def fn_aI(self, step=1):
        """Apical current in A / m^2 or C / (sec m^2)"""
        return self.co.F_RT * (- self.data["aJ_Na"][step-1] + self.data["aJ_Cl"][step-1] - self.data["aJ_K"][step-1])

    # Eq 7.3
    def fn_bI(self, step=1):
        """Basolateral Current in A / m^2 or C / (sec m^2)"""
        return self.co.F_RT * (- self.data["J_pump"][step-1] + self.data["bJ_Cl"][step-1] - self.data["aJ_K"][step-1])

    # Eq 8
    def fn_dATP(self, step=1):
        """ATP rate of change in micro-Molar per minute (muM / min)"""

        return self.co.J_atp / self.data["H"][step-1] - self.co.V_1_max * self.data["ATP"][step-1] / (
                self.co.K_1_m + self.data["ATP"][step-1]) \
               - self.co.V_2_max * self.data["ATP"][step-1] / (self.co.K_2_m + self.data["ATP"][step-1]) \
               - self.co.V_3_max * self.data["ATP"][step-1] / (self.co.K_3_m + self.data["ATP"][step-1]) \
               - self.co.V_10_max * self.data["ATP"][step-1] / (self.co.K_10_m + self.data["ATP"][step-1]) \
               - self.co.V_F_max / (1 +
                                    self.co.K_F_atp * self.co.K_F_amp /
                                        (self.data["ATP"][step-1] * self.data["AMP"][step-1])
                                    + self.co.K_F_atp / self.data["ATP"][step-1]
                                    + self.co.K_F_amp / self.data["AMP"][step-1]) \
               + self.co.V_B_max / (1 + (self.co.K_B_adp / self.data["ADP"][step-1]) ** 2
                                    + 2 * self.co.K_B_adp / self.data["ADP"][step-1]) \
               - self.data["ATP"][step-1] * self.data["dH"][step-1] / self.data["H"][step-1]

    # Eq 9
    def fn_dADP(self, step=1):
        """ADP rate of change in micro-Molar per minute (muM / min)"""

        return self.co.J_adp / self.data['H'][step-1] \
               + self.co.V_1_max * self.data["ATP"][step-1] / (self.co.K_1_m + self.data["ATP"][step-1]) \
               + self.co.V_2_max * self.data["ATP"][step-1] / (self.co.K_2_m + self.data["ATP"][step-1]) \
               + self.co.V_3_max * self.data["ATP"][step-1] / (self.co.K_3_m + self.data["ATP"][step-1]) \
               - self.co.V_4_max * self.data["ADP"][step-1] / (self.co.K_4_m + self.data["ADP"][step-1]) \
               - self.co.V_5_max * self.data["ATP"][step-1] / (self.co.K_5_m + self.data["ATP"][step-1]) \
               - 2 * self.co.V_F_max / (
                       1 + self.co.K_F_atp * self.co.K_F_amp / (self.data["ATP"][step-1] * self.data["AMP"][step-1])
                       + self.co.K_F_atp / self.data["ATP"][step-1] + self.co.K_F_amp / self.data["AMP"][step-1]) \
               + self.co.V_B_max / (1 + (self.co.K_B_adp / self.data["ADP"][step-1]) ** 2) \
               - 2 * self.co.V_B_max / (1 + self.co.K_B_adp / self.data["ADP"][step-1]) \
               - self.data["ADP"][step-1] * self.data['dH'][step-1] / self.data["H"][step-1]

    # Eq 10
    def fn_dAMP(self, step=1):
        """AMP's apical rate of change in micro-Molar per minute (muM / min)"""

        return self.data['dH'][step-1] * (self.data['AMP'][step-1] / self.data['H'][step-1]) \
               - self.co.V_8_max * self.data['AMP'][step-1] / (
                       self.co.K_8_m * (
                       1 + self.data['ATP'][step-1] / self.co.K_IN_atp + self.data['ADP'][step-1] / self.co.K_IN_adp)
                       + self.data['AMP'][step-1]) \
               - self.co.V_7_max * self.data["AMP"][step-1] / (
                       self.co.K_7_m * (
                       1 + self.data["ATP"][step-1] / self.co.K_IN_atp
                       + self.data["ADP"][step-1] / self.co.K_IN_adp + self.data["AMP"][step-1])
               ) \
               - self.co.V_6_max * self.data["AMP"][step-1] / (
                       self.co.K_6_m * (
                       1 + self.data["ATP"][step-1] / self.co.K_IN_atp + self.data["ADP"][step-1] / self.co.K_IN_adp
               ) + self.data["AMP"][step-1]) \
               + self.co.V_B_max / (1 + self.co.K_B_adp / self.data["ADP"][step-1]) \
               - self.co.V_F_max / (
                       1 + self.co.K_F_amp * self.co.K_B_adp / (self.data["ADP"][step-1] * self.data["AMP"][step-1]) +
                       self.co.K_F_atp / self.data["ATP"][step-1] + self.co.K_F_amp / self.data["AMP"][step-1]) \
               + self.co.V_10_max * self.data["ATP"][step-1] / (self.co.K_10_m + self.data["ATP"][step-1]) \
               + self.co.V_5_max * self.data["ATP"][step-1] / (self.co.K_5_m + self.data["ATP"][step-1]) \
               + self.co.V_4_max / (self.co.K_4_m + self.data["ADP"][step-1]) \
               + self.co.J_amp / self.data["H"][step-1]

    # Eq 11
    def fn_dADO(self, step=1):
        """ADO's apical rate of change in micro-Molar per minute (muM / min)"""

        return - self.data["ADO"][step-1] * self.data["dH"][step-1] / self.data["H"][step-1] \
               + self.co.V_8_max * self.data["AMP"][step-1] / (
                       self.co.K_8_m * (
                       1 + self.data["ATP"][step-1] / self.co.K_IN_atp + self.data["ADP"][step-1] / self.co.K_IN_adp
               ) + self.data["AMP"][step-1]) \
               + self.co.V_7_max * self.data["AMP"][step-1] / (
                       self.co.K_7_m * (
                       1 + self.data["ATP"][step-1] * self.data["ADP"][step-1] / (self.co.K_IN_adp * self.co.K_IN_atp)
                       + self.data["AMP"][step-1])
               ) \
               + self.co.V_6_max * self.data["AMP"][step-1] / (
                       self.co.K_6_m * (
                       1 + self.data["ATP"][step-1] / self.co.K_IN_atp + self.data["ADP"][step-1] / self.co.K_IN_adp
               ) + self.data["AMP"][step-1]) \
               - self.co.V_U1_max * self.data["ADO"][step-1] / (
                       self.co.K_U1_m + self.data["ADO"][step-1]
                       + self.data["INO"][step-1] * self.co.V_U1_max / self.co.K_U2_m
               ) \
               - self.co.V_9_max * self.data["ADO"][step-1] / (self.co.K_9_m + self.data["ADO"][step-1])

    # Eq 12
    def fn_dINO(self, step=1):
        """ADO's apical rate of change in micro-Molar per minute (muM / min)"""

        return (self.co.V_9_max * self.data["ADO"][step-1] / (self.co.K_9_m + self.data["ADO"][step-1])) \
               - self.co.V_U2_max * self.data["INO"][step-1] / (
                       self.co.K_U2_m + self.data["INO"][step-1]
                       + self.data["ADO"][step-1] * self.co.K_U2_m / self.co.K_U1_m) \
               - self.data["INO"][step-1] * self.data["dH"][step-1]


# Run the simulation with the parameters in init_d and save its data. Will fill the complete pd.DataFrame
def runAll(init_d=None, save_extra=True, sub_dir=SUB_DIR):
    """
    Will create a AirwayModel object and fill up all of it. Once the model is complete it
    will save the results into a csv file. If an error is encountered the run will exit at that point.
    Optionally it will save extra info into a secondary CSV (if there are no errors):
        * Initial values of the variables.
        * Average values of the variables.
        * Runtime.
        * How many seconds each step represented.
        * Total steps the simulation took.
    If there is any error (independently of save_extra), the secondary CSV will contain the following:
        * Initial values of the variables.
        * Values for the step previous to the error.
        * The values of the variables when the error was raised.
        * Runtime.
        * How many seconds each step represented.
        * Total steps the simulation took.
        * And an error code with the step at which it was received.

    :param sub_dir: The directory to which to save the CSV files, defaults to SUB_DIR (global in the simulation).
    :param init_d: The initial values of the DataFrame. If None, it will still create the model, which will ask for the
            values individually.
    :param save_extra: Boolean, whether to save the CSV with the extra info. If there are any issues encountered,
            this will be overwritten to True.
    :return: The name of the file to which the data was saved (NOT the path, just the file name).
    """

    date_time = datetime.datetime.now()
    time = datetime.time(hour=date_time.hour, minute=date_time.minute, second=date_time.second)
    print('Starting run.\n\tTime:', time)
    # Create d_extr to save all the averages and runtime
    d_extr = {}
    gen_data = AirwayModel(init_d)

    init_time = datetime.datetime.now()

    # Run all, but if an error is encountered, exit the run. It'll be easier to see the error in that case.
    # And save time if the error was caused only in one of the values due to initial conditions.
    for step in range(1, init_d['max_steps']):
        # noinspection PyBroadException
        try:
            gen_data.fn_run(step)
        except:
            print('Unexpected error', sys.exc_info()[0])
            d_extr['err'] = [sys.exc_info()[0], step, None]
            # Drop the unfilled rows, they consume space and just complicate things.
            gen_data.drop(labels=range(i+1, gen_data.max_steps))
            break

    end_time = datetime.datetime.now()
    runtime = end_time - init_time

    file_name = 'Data_' + str(datetime.datetime.now())[:-5].replace(':', '-').replace(' ', '_')

    if gen_data.isCF:
        file_name += '_CF'
    else:
        file_name += '_NL'

    # Save the file
    gen_data.to_csv(file_name + '.csv', sub_dir + '/')

    if 'err' in d_extr:
        # The step at which the error occurred
        step_err = d_extr['err'][1]
        # If exited from an error, save the relevant data: initial, final, and final-1
        for variable in (AirwayModel.variables[1::])[::-1]:
            d_extr[variable] = [gen_data[variable][0], gen_data[variable][step_err - 1],
                                gen_data[variable][step_err]]
        d_extr['runtime'] = [str(runtime), None, None]
        d_extr['time_frame'] = [gen_data.time_frame, None, None]
        d_extr['max_step'] = [gen_data.max_steps, None, None]
        # Save the extra data about the error
        pd.DataFrame(d_extr, index=['Initial', 'Previous', 'Final (err)']).to_csv(
            sub_dir + '/' + file_name + '__error-info.csv')
    elif save_extra:
        # To ignore 'Time (min)' we start from the index 1
        for variable in (AirwayModel.variables[1::])[::-1]:
            av = np.sum(gen_data[variable]) / gen_data.max_steps
            d_extr[variable] = [gen_data[variable][0], av]
        d_extr['runtime'] = str(runtime)
        d_extr['time_frame - max_step'] = [gen_data.time_frame, gen_data.max_steps]
        # Save the extra data about the run
        pd.DataFrame(d_extr, index=['initial', 'average']).to_csv(
            sub_dir + '/' + file_name + '__extra-info.csv')

    date_time = datetime.datetime.now()
    time = datetime.time(hour=date_time.hour, minute=date_time.minute, second=date_time.second)
    print('Done with the run.', '\n\tRuntime:', runtime, '\n\tAt time:', time)
    return file_name


# Will return all the files inside a directory.
def flsInDir(directory=SUB_DIR):
    """
    Will look at all the files in a directory and collect them into a list. Will walk using the OS Library
    :param directory: The path to the directory
    :return: A list with the paths to the files inside the directory. None if there are no files inside.
    """
    files_in = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            files_in.append(os.path.join(root, file))
    if len(files_in) == 0:
        return None
    else:
        return files_in


# Will delete the given folder, all of its contents and sub-contents.
def deleteFolder(directory=SUB_DIR):
    """
    Will delete the folder and all of it's contents, including sub-directories.
    :param directory: The directory to which to delete.
    :return: True if successful.
    """
    for root, dirs, files in os.walk(directory):
        for d in dirs:
            deleteFolder(d)
        for f in files:
            os.remove(os.path.join(root, f))
    os.rmdir(directory)
    return True


# Will save the file_name.csv and file_name_extra-info.csv files into one ZIP file.
def zipFiles(file_names=flsInDir(), zip_name=None, delete_after=False):
    """
    Will create a zip file in the current directory which contains the files in the file_names
    :param file_names: Should be an iterable with the path and file names of the files to compress.
            Each item must be the path and name of the file.
            The files should have an extension, if a dot character ('.') is not detected in the 4th or 3rd position
            from the back, it will attach the file extension '.csv'
    :param zip_name: (string) The name of the ZIP file, it will add a '.zip' extension if it's not already in the name.
            It defaults to None, in which case it will generate one with 'Data_' and the current date.
    :param delete_after: boolean whether or not to delete the files after compressing. Defaults to False
    :return: True if completed successfully
    """
    if zip_name is None:
        zip_name = 'Data_' + str(datetime.datetime.now())[0:-10].replace(':', '-').replace(' ', '_') + '.zip'
    elif '.zip' not in zip_name:
        zip_name += '.zip'
    zip_file = zipfile.ZipFile(zip_name, 'w')

    for file in file_names:
        if '.' not in file[-4:-2]:
            try:
                zip_file.write(file + '.csv')
            except FileNotFoundError:
                print('Could not find the file', file + '.csv')
        else:
            try:
                zip_file.write(file)
            except FileNotFoundError:
                print('Could not find the file', file)

    print(zip_name, 'was created.')

    if delete_after:
        for file in file_names:
            try:
                # Check to see if there's a file format attached
                if '.' not in file[-4::-2]:
                    os.remove(file + '.csv')
                else:
                    os.remove(file)
            except FileNotFoundError:
                pass
        return True
    else:
        return True


if __name__ == '__main__':
    # Empty the SUB_DIR folder before starting. If there are files inside, try to zip and delete them.
    if os.path.isdir(SUB_DIR):
        fls = flsInDir(SUB_DIR)
        if fls is not None:
            zipFiles(fls, 'OldData' + str(datetime.date.today()), delete_after=True)
        os.rmdir(SUB_DIR)

    # Try to load the 'initial-values.csv' file and run all of the simulations in parallel in sets of
    # MAX_PROCESSES using the multiprocess.Pool module.
    # If it can't be done (FileNotFound error), then run only one simulation, asking for all of the input values.
    try:
        # Prepare the initial values into a list of dictionaries.
        initial_values = []
        initial_pd = pd.read_csv('initial-values.csv')
        for i in range(len(initial_pd)):
            initial_values.append({})
            for var in initial_pd:
                if var == 'max_steps':
                    initial_values[-1][var] = int(initial_pd[var][i])
                elif var == 'CF':
                    initial_values[-1][var] = bool(initial_pd[var][i])
                else:
                    initial_values[-1][var] = float(initial_pd[var][i])

        print('Loaded the initial values')
        print('There will be', len(initial_values), 'simulations run. Will run parallel in sets of', MAX_PROCESSES)

        init = datetime.datetime.now()

        with Pool(MAX_PROCESSES) as p:
            p.map(runAll, initial_values)

        end_t = datetime.datetime.now()

        zipFiles(flsInDir(SUB_DIR), delete_after=False)

        av_steps = 0
        for r in initial_values:
            av_steps += r['max_steps']
        av_steps /= len(initial_values)

        print('Did', len(initial_values), 'simulations, with an average of', av_steps,
              'steps in', end_t - init)
    except FileNotFoundError:
        runAll()

    print('\nDone!')
