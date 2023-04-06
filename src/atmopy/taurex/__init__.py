"""Module containing TauREx Plugin."""
from astropy import units as u
from taurex.cache import GlobalCache
from taurex.chemistry import AutoChemistry
from taurex.core import fitparam

from atmopy.io import baseline_molecules
from atmopy.runner import ATMORunner


class ATMOChemistry(AutoChemistry):
    """Taurex Atmospheric Chemistry ATMO class."""

    def __init__(self, metallicity=0.0, cond_h2o=False):
        """Initiialize."""
        super().__init__("ATMO")
        self.atmorunner = ATMORunner(
            atmo_path=GlobalCache()["atmo_path"],
            atmo_executable=GlobalCache()["atmo_executable"],
        )
        self.atmorunner.chemistry.metallicity = metallicity
        self.atmorunner.chemistry.condensation_h2o = cond_h2o

        self.mix = None

        self.element_ratios = []

        self.determine_active_inactive()

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
        return self.atmorunner.chemistry.metallicity

    @metallicity.setter
    def metallicity(self, value):
        """Setter."""
        self.atmorunner.chemistry.metallicity = value

    @classmethod
    def input_keywords(cls):
        """Taurex detection keywords."""
        return ["ATMO", "atmo", "frogr"]
