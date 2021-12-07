#
# Created on Sat Mar 06 2021
#
# Copyright (c) 2021 Your Company
# Name: Luka Jeromel
#
# ******************************Python imports***********************************
from typing import Union

# ******************************PyQt5 imports*************************************
# ******************************Other third party imports************************
import numpy as np

# ******************************My modules***************************************

# Required data.
pKa_dict = {
    "Acetic acid": np.array([4.756], dtype=np.float32),
    "Citric acid": np.array([3.128, 4.761, 6.396], dtype=np.float32),
    "Hydrocloric acid": np.array([-4], dtype=np.float32),
    "Hydorcyanic acid": np.array([9.31], dtype=np.float32),
    "Nitric acid": np.array([-1], dtype=np.float32),
    "Phosphate": np.array([2.147, 7.207, 12.346], dtype=np.float32),
    "Succinic acid": np.array([4.2, 5.6], dtype=np.float32),
    "Malate": np.array([3.4, 5.13], dtype=np.float32),
    "TRIS": np.array([8.07], dtype=np.float32),
    "MES": np.array([6.21], dtype=np.int32),
    "EDTA": np.array([2.0, 2.7, 6.2, 10.31], dtype=np.float32),
    "Cysteine": np.array([1.91, 10.28, 8.14], dtype=np.float32),
    "Histidine": np.array([1.7, 9.09, 6.04], dtype=np.float32),
    "Arginine": np.array([2.09, 9, 12.1], dtype=np.float32),
    "Glycine": np.array([2.35], dtype=np.float32),
}

# charges of conjugate (bases acid -> conjuget_base + H^+ ---- HA -> A^- + H^+)
charges_dict = {
    "OH": np.array([-1], dtype=np.int32),
    "H": np.array([1], dtype=np.int32),
    "Na": np.array([1], dtype=np.int32),
    "Ca": np.array([-1], dtype=np.int32),
    "Acetic acid": np.array([-1], dtype=np.int32),
    "Citric acid": np.array([-1, -2, -3], dtype=np.int32),
    "Hydrocloric acid": np.array([-1], dtype=np.int32),
    "Hydorcyanic acid": np.array([-1], dtype=np.int32),
    "Nitric acid": np.array([-1], dtype=np.int32),
    "Phosphate": np.array([-1, -2, -3], dtype=np.int32),
    "Succinic acid": np.array([-1, -2], dtype=np.int32),
    "Malate": np.array([0, -1], dtype=np.int32),  # not ok - have subtract 1, check
    "TRIS": np.array([1], dtype=np.int32),
    "MES": np.array([-1], dtype=np.int32),
    "EDTA": np.array([0, -1, -2, -3], dtype=np.int32),  # not ok
    "Cysteine": np.array([0, 1, 0], dtype=np.int32),  # not ok
    "Histidine": np.array([0, 1, 1], dtype=np.int32),  # not ok
    "Arginine": np.array([0, 1, 1], dtype=np.int32),  # not ok
    "Glycine": np.array([-1], dtype=np.int32),
}


class Buffer(object):
    def __init__(
        self, added_species_molarities_M: dict, pH: np.float32, volume_L: np.float32
    ) -> None:
        # Initial we create water solution with no Na or Cl ions. We will later change
        # this so we reach desired pH. In reality other ions may be used.
        self.__water_solution = WaterSolution(added_species_molarities_M, 0.0, 0.0)
        self.__volume_L = volume_L
        self.__pH = pH

        # TODO Old, to be removed.
        # self.__added_species_conc_M = added_species_molarities_M
        # self.__Na_conc_M = 0.0
        # self.__Cl_conc_M = 0.0

    def __calc_needed_Na_Cl(self):
        """Based on current water solution, it calculates how much additional Na+
        or Cl- ions we need to reach pH of buffer.
        """
        added_mol_M = self.__water_solution.get_added_chemical_species_M()
        added_NaCl_M = self.__water_solution.get_added_NaCl_M()
        phosphate_type = 1
        # This is dependant on specific chemicals that would be used.
        Na_min = (
            added_NaCl_M + phosphate_type * added_mol_M["Phosphate"] + 2 * added_mol_M["EDTA"]
        )
        Cl_min = added_NaCl_M + added_mol_M["Cysteine"] + added_mol_M["Arginine"]
        ion_str_M = 0.5 * (Na_min + Cl_min)
        pass

    def calculate_recipe(self, chemicals: dict) -> None:
        pass


class WaterSolution(object):
    """Class representing water solution. When water soultion is created we pour
    in different chemical species. Based
    Changing concentration (molar) of sodium and chlorine ions mimics up and down
    titration of the solution.
    """

    def __init__(
        self, molarities_dict: dict, conc_Na_M: np.float32, conc_Cl_M: np.float32
    ) -> None:
        self.__added_species_conc_M = molarities_dict
        self.__added_NaCl_M = 0.0
        self.__pKa_data = pKa_dict
        self.__charges_data = charges_dict
        self.__conc_Na_M = conc_Na_M
        self.__conc_Cl_M = conc_Cl_M
        self.__ion_molarities = None
        self.__pH = None

        self.__config = {
            "min_pH": 0.0,
            "max_pH": 14,
            "min_ion_strength_M": 0.00001,
            "max_ion_strength_M": 5.5,
        }

        # We calculate initial solution's ion molarities.
        # self.__calculate_molarities()

    def get_Na_conc_M(self) -> np.float32:
        return self.__conc_Na_M

    def get_Cl_conc_M(self) -> np.float32:
        return self.__conc_Cl_M

    def get_added_chemical_species_M(self) -> dict:
        return self.__added_species_conc_M

    def set_Na_conc_M(self, conc_Na_M: np.float32) -> None:
        self.__conc_Na_M = conc_Na_M

    def set_Cl_conc_M(self, conc_Cl_M: np.float32) -> None:
        self.__conc_Cl_M = conc_Cl_M

    def add_NaCl_M(self, conc_NaCl_M: np.float32) -> None:
        self.__added_NaCl_M = conc_NaCl_M

    def __ion_strength(self) -> np.float32:
        """Calculates ion strength of current water solution.
        """
        ion_strength = 0.0
        for species in self.__ion_molarities:
            ion_strength += (
                self.__charges_data[species] ** 2 * self.__ion_molarities[species]
            ).sum()

        return (ion_strength + self.__conc_Na_M + self.__conc_Cl_M) / 2.0

    def calc_electro_neutrality_error(self) -> np.float32:
        return self.__electro_neutrality_error()

    def __electro_neutrality_error(self) -> np.float32:
        error = 0.0
        for species in self.__ion_molarities:
            error += (self.__charges_data[species] * self.__ion_molarities[species]).sum()

        return error + self.__conc_Na_M - self.__conc_Cl_M

    def __ion_strength_error(self, test_ion_strength: np.float32) -> np.float32:
        return test_ion_strength - self.__ion_strength()

    @staticmethod
    def __effective_pKa(pKa: np.float32, IS_M: np.float32, charge: np.float32) -> np.float32:
        """Corrects dissociation constant defined at 298 K and infinite dilution
        for  solution ionic strength using extended Debye-Huckel relationship proposed
        by Davies (see Pabst Carta 2007 articel page 21).

        @param pKa Origina equilibrium constant.
        @param IS_M Molar ionic strength of solution.

        
        """
        # Constant incorporating temperature effects
        A = 0.5114  # at 298 K
        # Buffer-dependent constant, e. g. 0.16 for acetate and 0.16 for phosphate
        b = 0.1
        # b was set to 0.1 if no measurements were found.
        return pKa + (2 * (charge) * (A * np.sqrt(IS_M) / (1 + np.sqrt(IS_M)) - b * IS_M))

    def __calculate_molarities(
        self, pH: np.float32, ion_strength: np.float32, correct_for_IS: bool = True,
    ) -> None:
        """Calculates molarites of all charged species (ions) in water solution 
        at given pH and ion_stregth.
        @param species_molarities"""
        # concentration of H+ ions
        conc_H_M = 10 ** (-pH)
        # Water dissociation constant (could move it be a class constant).
        K_w = 10 ** (-14)
        # Concentration of hydroxyl ions.
        conc_OH_M = K_w / conc_H_M

        # Output dictionary.
        ion_molarities = {
            "H": np.array([conc_H_M]),
            "OH": np.array([conc_OH_M]),
        }

        for species in self.__added_species_conc_M:
            ion_molarities[species] = np.zeros(len(pKa_dict[species]))
            denom = 1
            term = 1

            for charge, pKa in zip(self.__charges_data[species], self.__pKa_data[species]):
                pKa_eff = pKa
                if correct_for_IS:
                    pKa_eff = self.__effective_pKa(pKa, ion_strength, charge)
                term *= 10 ** (-pKa_eff) / conc_H_M
                denom += term
            molarity = self.__added_species_conc_M[species] / denom

            for i, pKa in enumerate(self.__pKa_data[species]):
                pKa_eff = pKa
                if correct_for_IS:
                    pKa_eff = self.__effective_pKa(pKa, ion_strength, charge)
                molarity *= 10 ** (-pKa_eff) / conc_H_M
                ion_molarities[species][i] = molarity
        # print(f"{ion_molarities=}")
        self.__ion_molarities = ion_molarities

    def __calc_pH_bisect(
        self,
        pH_low: np.float32,
        pH_high: np.float32,
        ion_strength_M: np.float32,
        correct_for_IS=False,
    ) -> np.float32:
        """Calculates pH from input molarities.
        """
        eps = 1e-8
        step_counter = 0
        max_steps = 100
        error = 1
        pH_low = pH_low
        pH_high = pH_high
        while abs(error) > eps and step_counter < max_steps:
            pH = (pH_low + pH_high) / 2
            self.__calculate_molarities(
                pH, ion_strength_M, correct_for_IS=correct_for_IS,
            )
            # print(f"{ion_strength=}\n{c_Na=}\n{c_Cl=}")
            # print(f"{self.__ion_molarities=}")
            error = self.__electro_neutrality_error()
            if error > 0:
                pH_low = pH
            else:
                pH_high = pH
            # print("i: {} pH {:.4f} error: {:.4f}".format(step_counter, pH, error))
            step_counter += 1
        return pH

    def __calc_ion_strength_bisect(
        self,
        ion_strength_low_M: np.float32,
        ion_strength_high_M: np.float32,
        pH: np.float32,
        correct_for_IS=True,
    ) -> np.float32:
        # TODO Change algorithm to Newton.
        eps = 1e-8
        step_counter = 0
        max_steps = 100
        error = 1
        while abs(error) > eps and step_counter < max_steps:
            ion_strength_M = (ion_strength_low_M + ion_strength_high_M) / 2
            self.__calculate_molarities(
                pH, ion_strength_M, correct_for_IS=correct_for_IS,
            )
            error = self.__ion_strength_error(ion_strength_M)
            if error < 0:
                ion_strength_low_M = ion_strength_M
            else:
                ion_strength_high_M = ion_strength_M
            step_counter += 1
        return ion_strength_M

    def __calculate_pH_IS(
        self, initial_pH: np.float32, verbose=False,
    ):
        IS_min = 0.00001
        IS_max = 5.5
        pH_min = 0.0
        pH_max = 14.0
        error = 1e-5

        # TODO some smart guess of initial_pH
        IS = self.__calc_ion_strength_bisect(IS_min, IS_max, initial_pH, True,)
        pH = self.__calc_pH_bisect(pH_min, pH_max, IS, True,)
        if verbose:
            print(f"\nInitial values, pH: {pH:8.5f} IS: {IS:8.5f}")

        iteration = 1
        previous_pH = initial_pH
        while abs(pH - previous_pH) > error and iteration < 1000:
            IS = self.__calc_ion_strength_bisect(IS_min, IS_max, pH, True)
            previous_pH = pH
            pH = self.__calc_pH_bisect(pH_min, pH_max, IS, True)

            if verbose:
                print(
                    f"step: {iteration} pH: {pH:8.5f} IS: {IS:8.5f}, error: "
                    f"{abs(pH-previous_pH):9.7f}"
                )
            iteration += 1
        return pH, IS

    def calculate(self, init_pH, verbose=False):
        self.__calculate_pH_IS(init_pH, verbose=verbose)
        if verbose:
            print(f"{self.__ion_molarities=}")


def test():
    c_Na = 0.0
    c_Cl = 0.02

    # molarities are in mol/L
    molarities_dict = {
        "Acetic acid": 0.1,
        # "Citric acid": 0.1,
        # "Phosphate": 0.3,
    }

    water_solution = WaterSolution(molarities_dict, 0.0, 0.0)
    water_solution.calculate(init_pH=7.0, verbose=True)
    water_solution.set_Na_conc_M(c_Na)
    water_solution.set_Cl_conc_M(c_Cl)
    water_solution.calculate(init_pH=7.0, verbose=True)
    #%timeit pH, IS = calculate_pH_IS(molarities_dict, c_Na, c_Cl, 1, pKa_dict, charges_dict, verbose =False)
    # pH, IS = calculate_pH_IS(
    #    molarities_dict, c_Na, c_Cl, 1, pKa_dict, charges_dict, verbose=False
    # )
    # print("result: ", pH, IS)


if __name__ == "__main__":
    test()
