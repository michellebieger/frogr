"""Handles running the pre-installed atmo for TauREx."""
import os
import subprocess as sp
import typing as t
from subprocess import PIPE

import numpy as np
import numpy.typing as npt
from astropy import units as u

from .io import baseline_masses


class ATMOOutput(t.TypedDict):
    """Output format for ATMO runs."""

    abundances: npt.NDArray[np.float64]
    abundances_vmr: npt.NDArray[np.float64]
    pressure: u.Quantity
    temperature: u.Quantity
    mean_mol_mass: npt.NDArray[np.float64]
    molname: t.List[str]
    element_abundances: npt.NDArray[np.float64]


class ATMOFailureError(Exception):
    """Exception when ATMO fails to generate abundances."""

    pass


def run_atmo(
    input_file: str,
    working_directory: t.Optional[str] = None,
    atmo_path: t.Optional[str] = None,
    atmo_executable: t.Optional[str] = "atmo.x",
    copy_atmo: t.Optional[bool] = True,
):
    """Run ATMO for a given input."""
    import os
    import shutil

    atmo_path = atmo_path or os.environ.get("ATMO_PATH", None)
    atmo_executable = atmo_executable or os.environ.get("ATMO_EXE", "atmo.x")

    working_directory = working_directory or atmo_path

    total_atmo_path = os.path.join(atmo_path, atmo_executable)
    # print(total_atmo_path)
    if copy_atmo:
        new_path = os.path.join(working_directory, atmo_executable)
        shutil.copyfile(total_atmo_path, new_path)
        shutil.copymode(total_atmo_path, new_path)
        total_atmo_path = new_path

    process = sp.Popen(  # noqa: S404, S603
        [total_atmo_path, input_file],
        cwd=working_directory,
        stdout=PIPE,
        stderr=PIPE,
    )  # noqa: S404

    stdout, stderr = process.communicate()
    if process.returncode != 0:
        # if a program exits on an error or fail condition its exitcode is usually not 0
        raise RuntimeError(f"ATMO did not run right!\n\n{process.stderr}")

    if f"{input_file} does not exist" in stdout.decode("utf-8"):
        raise ValueError("Incorrect input file specified")

    return stdout.decode("utf-8"), stderr.decode("utf-8")


class ATMORunner:
    """Runs ATMO."""

    def __init__(
        self,
        atmo_path: t.Optional[str] = None,
        atmo_executable: t.Optional[str] = "atmo.x",
    ):
        """Initialize."""
        from .io import ChemistryInputSection

        self.chemistry = ChemistryInputSection()
        self.atmo_path = atmo_path
        self.atmo_executable = atmo_executable

    def _run(
        self,
        temperature: u.Quantity,
        pressure: u.Quantity,
        run_directory: str,
        output_filename: t.Optional[str] = "chem_out.ncdf",
    ) -> ATMOOutput:
        """Run atmo in a certain directory."""
        from scipy.io import netcdf_file

        from .io import ParamInputSection
        from .io import convert_molecule_name_to_string
        from .io import generate_input_file

        param_input = ParamInputSection(pressure=pressure, temperature=temperature)

        param_section = param_input.build_section(run_directory)

        nlevels = temperature.size - 1
        chem_section = self.chemistry.build_section(
            run_directory, output_name=output_filename
        )

        full_input_file, input_file = generate_input_file(
            run_directory, nlevels=nlevels, sections=[param_section, chem_section]
        )

        stdout, stderr = run_atmo(
            input_file,
            working_directory=run_directory,
            atmo_path=self.atmo_path,
            atmo_executable=self.atmo_executable,
        )

        expected_output_path = os.path.join(run_directory, output_filename)
        if not os.path.isfile(expected_output_path):
            #
            print(stdout)
            print(stderr)
            raise ATMOFailureError("ATMO did not finish or run properly")

        with netcdf_file(expected_output_path, "r") as f:
            abundances = f.variables["abundances"][...].copy()
            pressure = f.variables["pressure"][...].copy() << u.dyn / u.cm**2
            temperature = f.variables["temperature"][...].copy() << u.K
            mean_mol_mass = f.variables["mean_mol_mass"][...].copy()
            molname = convert_molecule_name_to_string(
                f.variables["molname"][...].copy()
            )
            element_abundances = f.variables["element_abundances"][...].copy()

        mol_masses = baseline_masses()

        abund_vmr = abundances * mean_mol_mass[None, :] / mol_masses[:, None]

        abund_vmr = abund_vmr / abund_vmr.sum(axis=0)[None, :]

        return {
            "abundances": abundances,
            "abundances_vmr": abund_vmr,
            "pressure": pressure,
            "temperature": temperature,
            "mean_mol_mass": mean_mol_mass,
            "molname": molname,
            "element_abundances": element_abundances,
        }

    def run(
        self,
        temperature: u.Quantity,
        pressure: u.Quantity,
        run_directory: t.Optional[str] = None,
        output_filename: t.Optional[str] = "chem_out.ncdf",
    ) -> ATMOOutput:
        """Run ATMO."""
        import tempfile

        if run_directory:
            return self._run(
                temperature,
                pressure,
                run_directory,
                output_filename,
            )
        else:
            with tempfile.TemporaryDirectory() as f:
                return self._run(
                    temperature,
                    pressure,
                    str(f),
                    output_filename,
                )
