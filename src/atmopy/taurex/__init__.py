"""Module containing TauREx Plugin."""
from astropy import units as u
from taurex.cache import GlobalCache
from taurex.chemistry import AutoChemistry
from taurex.core import fitparam

from atmopy.io import baseline_molecules, baseline_ratio
from atmopy.util import compute_element_factor
from atmopy.runner import ATMORunner
import typing as t
from atmopy.util import compute_element_factor


class ATMOChemistry(AutoChemistry):
    """Taurex Atmospheric Chemistry ATMO class."""

    def __init__(
        self,
        metallicity=1.0,
        ratio_elements=None,
        ratio_factors=None,
        baseline_element="O",
        cond_h2o=False,
        cond_nh3=False,
    ):
        """Initiialize."""
        super().__init__("ATMO")
        self.atmorunner = ATMORunner(
            atmo_path=GlobalCache()["atmo_path"],
            atmo_executable=GlobalCache()["atmo_executable"],
        )

        ratio_elements = ratio_elements or []
        ratio_factors = ratio_factors or []
        self.baseline_element = baseline_element

        self.current_ratios = baseline_ratio(self.baseline_element)

        if len(ratio_factors) != len(ratio_elements):
            raise ValueError("Ratio elements and factors must be same size.")

        element_factors = list(zip(ratio_elements, ratio_factors))

        self.atmorunner.chemistry.metallicity = 1.0

        self._metallicity = metallicity

        element_factors = element_factors or []

        for elem, val in element_factors:
            self.current_ratios[elem] = val
        # set element ratio factors to give to TR
        self.atmorunner.chemistry.condensation_h2o = cond_h2o
        self.atmorunner.chemistry.condensation_nh3 = cond_nh3
        # set condensation factors to give to TR
        self.mix = None

        self.recompute_elements()
        # compute element ratios again post setting a metallicity
        self.determine_active_inactive()
        self.build_ratio_params()

    def recompute_elements(self):
        """Remcompute metallicity and ratios"""
        element_factors = compute_element_factor(
            self._metallicity, self.current_ratios.items()
        )
        self.atmorunner.chemistry.element_factor = element_factors

    def build_ratio_params(self):
        """Create fitting parameters for the elements"""
        from taurex.util.util import molecule_texlabel

        for element in self.current_ratios.keys():
            if element in ("H", "He", self.baseline_element):
                continue
            mol_name = f"{element}_{self.baseline_element}_ratio"
            param_name = mol_name
            param_tex = "{}/{}".format(
                molecule_texlabel(element), molecule_texlabel(self.baseline_element)
            )
            # for-loop to create fitting parameter for each element ratio that you have selected
            def read_mol(self, element=element):
                return self.current_ratios[element]

            def write_mol(self, value, element=element):
                self.current_ratios[element] = value

            # functions to read, write so that they're passed back to TR
            read_mol.__doc__ = f"{element}/{self.baseline_element} ratio"

            fget = read_mol
            fset = write_mol

            bounds = [1.0e-12, 0.1]
            # create default bounds
            default_fit = False
            self.add_fittable_param(
                param_name, param_tex, fget, fset, "linear", default_fit, bounds
            )

    @property
    def gases(self):
        """Return gases in atmosphere."""
        return baseline_molecules()

    def initialize_chemistry(
        self,
        nlayers=100,
        temperature_profile=None,
        pressure_profile=None,
        altitude_profile=None,
    ):
        """Build chemistry."""
        self.recompute_elements()
        temperature = temperature_profile << u.K
        pressure = pressure_profile << u.Pa
        result = self.atmorunner.run(temperature, pressure)

        self.mix = result["abundances_vmr"]
        super().compute_mu_profile(nlayers=nlayers)

    @property
    def mixProfile(self):  # noqa: N802
        """Get VMR profiles for each molecule."""
        return self.mix

    @fitparam(param_name="metallicity", param_latex=r"$Z_{o}$")
    def metallicity(self):
        """Metallicity of atmosphere in solar metallicity units."""
        return self._metallicity

    @metallicity.setter
    def metallicity(self, value):
        """Setter."""
        self._metallicity = value

    @classmethod
    def input_keywords(cls):
        """Taurex detection keywords."""
        return ["atmopy", "ATMO", "atmo", "frogr"]


class ATMONeqChemistry(ATMOChemistry):
    pass
