"""Handles building input files and namelists."""
import typing as t
from dataclasses import asdict
from dataclasses import dataclass
from enum import Enum
from functools import lru_cache

import numpy as np
import numpy.typing as npt
import pkg_resources
from astropy import units as u


def to_fortran_bool(value: bool) -> str:
    """Convert Python bool to Fortran bool."""
    if value:
        return ".true."
    else:
        return ".false."


def create_namelist(
    title: str, params: t.Optional[t.Dict[str, t.Union[float, str, int, bool]]] = None
) -> str:
    """Creates namelist text string.

    Args:
        title: header of namelist
        params: values to put in

    Returns:
        Namelist format of parameters.

    """
    base = [f"&{title.upper()}"]
    namelist_params = []

    params = params or {}

    for k, v in params.items():
        value_text = None
        if isinstance(v, str):
            value_text = f"'{v}'"
        elif isinstance(v, bool):
            value_text = to_fortran_bool(v)
        else:
            value_text = f"{v}"
        namelist_params.append(f"{k} = {value_text}")

    return "\n".join(base + namelist_params + ["/"])


def create_element_factors(
    element_factors: t.List[t.Tuple[str, float]],
    elements: t.Optional[t.List[str]] = None,
) -> str:
    """Build ele_xfactor in namelist."""
    elements = elements or baseline_elements()

    return "\n".join(
        [
            f"ele_xfactor({elements.index(elem)+1}) = {value}"
            for elem, value in element_factors
        ]
    )


class ChemistryType(str, Enum):
    """Define chemsitry type."""

    EQUIL = ("eq",)
    NONEQUIL = ("neq",)


@dataclass
class ChemistryInputSection:
    """Builds chemistry input.

    Takes in arguments for chem namelist.
    Uses create_namelist to convert from user input to Fortran.
    """

    chem: t.Optional[ChemistryType] = ChemistryType.EQUIL
    """The switch for whether you want NEquil or not"""

    fcoeff: t.Optional[str] = "coeff_NASA_sc.dat"
    """NASA polynomial file: users should already have this from their atmo install"""
    fcoeffnine: t.Optional[str] = "coeff_NASA_sc.dat"
    """NASA 9 polynomial file: users should already have this from their atmo install"""

    element_factor: t.Optional[t.List[t.Tuple[str, float]]] = None

    metallicity: t.Optional[float] = 0.0
    """In log 10 units mate"""

    tfreeze_eq: t.Optional[float] = 0.0

    condensation_nh3: t.Optional[bool] = False
    condensation_h2o: t.Optional[bool] = False

    gibbs_step_size: t.Optional[float] = 2

    rainout: t.Optional[bool] = False

    def move_dummy(
        self, output_directory: str, extra_molecules: t.Optional[t.List[str]] = None
    ):
        """Moves chem_dummy to working folder."""
        import os
        import shutil

        extra_molecules = extra_molecules or []

        chem_dummy_path = os.path.join(output_directory, "chem_dummy.ncdf")
        base_chem_dummy_path = pkg_resources.resource_filename(
            "atmopy", os.path.join("data", "chem", "chem_dummy.ncdf")
        )

        if not extra_molecules:
            shutil.copyfile(base_chem_dummy_path, chem_dummy_path)
            shutil.copymode(base_chem_dummy_path, chem_dummy_path)
        else:
            raise NotImplementedError

        return None

    def copy_nasas(self, output_directory: str) -> None:
        """Copy NASA files to working folder."""
        copy_pkg_file(self.fcoeff, output_directory)
        copy_pkg_file(self.fcoeffnine, output_directory)

    def build_section(
        self,
        output_directory: str,
        output_name: t.Optional[str] = "chem_out.ncdf",
        chem_dummy_file: t.Optional[str] = None,
        extra_molecules: t.Optional[t.List[str]] = None,
        flat: t.Optional[bool] = True,
    ) -> str:
        """Convert plugin chemistry input into FORTRAN namelist."""
        import os

        my_output = asdict(self)

        output_filename = (
            os.path.join(output_directory, output_name) if not flat else output_name
        )

        self.move_dummy(
            output_directory=output_directory, extra_molecules=extra_molecules
        )

        self.copy_nasas(output_directory)

        my_output["fAin"] = "chem_dummy.ncdf"
        my_output = {k: v for k, v in my_output.items() if v is not None}
        # Remove nones

        if self.chem == ChemistryType.EQUIL:
            my_output["fAeqout"] = output_filename
        else:
            my_output["fAneqout"] = output_filename

        mappings = {
            "metallicity": "MdH",
            "condensation_nh3": "cond_NH3",
            "condensation_h2o": "cond_H2O",
            "gibbs_step_size": "chem_conv",
        }
        # dictionary from what we called it to what the namelist param is in atmo

        my_output = {mappings.get(k, k): v for k, v in my_output.items()}
        # check atmo inputs against the namelist params from atmo
        if self.element_factor:
            my_output["elem_xfactor"] = create_element_factors(self.element_factor)

        return create_namelist("chemistry", my_output)


@dataclass
class ParamInputSection:
    """Generate PARAM section for ATMO."""

    pressure: u.Quantity
    temperature: u.Quantity
    debug: t.Optional[bool] = False

    def build_section(self, directory: str, flat: t.Optional[bool] = True) -> str:
        """Build input temp/pressure ncdf for PARAM input to ATMO."""
        import os

        from scipy.io import netcdf

        if self.pressure.shape != self.temperature.shape:
            raise ValueError("Temperature and pressure shapes do not match")

        filename_path = os.path.join(directory, "pt_example.ncdf")
        with netcdf.netcdf_file(filename_path, "w") as f:
            f.createDimension("nlevel", self.temperature.size)
            temperature = f.createVariable("temperature", np.float64, ("nlevel",))
            temperature[:] = self.temperature.to(u.K).value

            pressure = f.createVariable("pressure", np.float64, ("nlevel",))
            pressure[:] = self.pressure.to(u.dyn / u.cm**2).value
        return create_namelist(
            "param",
            {
                "debug": int(self.debug),
                "fin": filename_path if not flat else "pt_example.ncdf",
            },
        )


def generate_input_file(
    directory: str,
    nlevels: t.Optional[int] = 100,
    output_name: str = "input.in",
    sections: t.List[str] = None,
) -> t.Tuple[str, str]:
    """Generate ATMO input file.

    Args:
        directory: Directory to write input file
        nlevels: Number of atmospheric levels in atmosphere.
        output_name: Output filename
        sections: List of sections to include

    Returns:
        Returns full path and just filename.

    """
    import os

    sections = sections or []

    full_string = "\n\n".join(
        [baseline_input_file().replace("__NLEVELS__", f"{nlevels}")] + sections
    )

    filename = os.path.join(directory, output_name)
    with open(filename, "w") as f:
        f.write(full_string + "\n")

    return filename, output_name


def convert_molecule_name_to_string(string_names) -> t.List[str]:
    """Converts S1 2D array to list of strings."""
    return [
        "".join(list(map(lambda x: x.decode("utf-8"), molecule))).strip()
        for molecule in string_names
    ]


@lru_cache(maxsize=1)
def baseline_molecules() -> t.List[str]:
    """Return base molecules in ATMO."""
    import os

    from scipy.io import netcdf_file

    with netcdf_file(
        pkg_resources.resource_filename(
            "atmopy", os.path.join("data", "chem", "chem_dummy.ncdf")
        )
    ) as f:
        molecules = f.variables["molname"][:].copy()

    return convert_molecule_name_to_string(molecules)


@lru_cache(maxsize=1)
def baseline_masses() -> npt.NDArray[np.float64]:
    """Return base molecule masses in AMU."""
    from taurex.util.util import calculate_weight

    return np.array([calculate_weight(m) for m in baseline_molecules()])


@lru_cache(maxsize=1)
def baseline_input_file() -> str:
    """Return base input file for build new ones."""
    import os

    base_file = pkg_resources.resource_filename(
        "atmopy", os.path.join("data", "chem", "general_input.in")
    )

    with open(base_file) as f:
        return f.read()


def baseline_elements():
    """Return baseline elements in ATMO."""
    return [
        "H",
        "He",
        "C",
        "N",
        "O",
        "Na",
        "K",
        "Si",
        "Ar",
        "Ti",
        "V",
        "S",
        "Cl",
        "Mg",
        "Al",
        "Ca",
        "Fe",
        "Cr",
        "Li",
        "Cs",
        "Rb",
        "F",
        "P",
    ]


def copy_pkg_file(filename: str, output_directory: str) -> None:
    """Copy a file from the package into a directory."""
    import os
    import shutil

    import pkg_resources

    file_path = pkg_resources.resource_filename(
        "atmopy", os.path.join("data", "chem", filename)
    )

    shutil.copy(file_path, os.path.join(output_directory, filename))
    shutil.copymode(file_path, os.path.join(output_directory, filename))
    return None
