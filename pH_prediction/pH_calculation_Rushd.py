# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 16:04:35 2019

@author: jeromlu2

based on Rushd Fortran script
"""


# Python imports
import re
import subprocess
import os.path


# Global variables
FOLDER = "./testing/"
# Input filenames
FNAME_IN_TEMP = "out_temp.dat"
FNAME_IN_COMP = "out_A&B_sngl.dat"
# Output filename
FNAME_OUT = "param_all.dat"
FORTRAN_EXE = "calc_pH_values.exe"


class FortranCommunication(object):
    def __init__(self, parent=None, folder=FOLDER):
        super(FortranCommunication, self).__init__()
        self.folder = folder

    def write_parameters(self, write_lst):
        """Writes to "param_all.dat" """
        with open(self.folder + FNAME_OUT, "w") as file:
            for element in write_lst:
                write_string = element + "   " + element + "\n"
                file.write(write_string)

    def run_subprocess(self):
        """Runs Fortran subprocess. Current executable: calc_pH_values.exe"""
        subprocess.call(self.folder + FORTRAN_EXE, cwd=FOLDER)

    def read_data(self, parent):

        parent.clear_data()

        with open(self.folder + FNAME_IN_COMP, "r") as file:
            for line in file:
                reg_exp = "\s+"
                value = re.split(reg_exp, line.strip())[0]
                if value == "NaN":
                    return False
                frtd_val = "{0:.2f} g".format(float(value))
                parent.composition.append(frtd_val)

        with open(self.folder + FNAME_IN_TEMP, "r") as file:
            for line in file:
                reg_exp = "\s+"
                values = re.split(reg_exp, line.strip())
                parent.temp_dep["Temp"].append(float(values[0]) - 273.15)
                parent.temp_dep["IS"].append(float(values[1]))
                parent.temp_dep["pH"].append(float(values[2]))
                parent.temp_dep["Cond"].append(float(values[3]))

        return True

    def check_if_exe_is_accesible(self):
        """
        Check if fortran exe file is present in the required self.folder
        """
        if not os.path.isfile(self.folder + FORTRAN_EXE):
            print("No exe file")
            return
