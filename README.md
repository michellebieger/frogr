# frogr: from gas to retrieval

[![PyPI](https://img.shields.io/pypi/v/atmopy.svg)][pypi_]
[![Status](https://img.shields.io/pypi/status/atmopy.svg)][status]
[![Python Version](https://img.shields.io/pypi/pyversions/atmopy)][python version]
[![License](https://img.shields.io/pypi/l/atmopy)][license]

[![Read the documentation at https://frogr-from-gas-to-retrieval.readthedocs.io/](https://img.shields.io/readthedocs/atmopy/latest.svg?label=Read%20the%20Docs)][read the docs]
[![Tests](https://github.com/mfbieger/atmopy/workflows/Tests/badge.svg)][tests]
[![Codecov](https://codecov.io/gh/mfbieger/atmopy/branch/main/graph/badge.svg)][codecov]

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)][pre-commit]
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)][black]

[pypi_]: https://pypi.org/project/frogr/
[status]: https://pypi.org/project/frogr/
[python version]: https://pypi.org/project/frogr
[read the docs]: https://frogr-from-gas-to-retrieval.readthedocs.io/
[tests]: https://github.com/michellebieger/frogr/actions?workflow=Tests
[codecov]: https://app.codecov.io/gh/michellebieger/frogr
[pre-commit]: https://github.com/pre-commit/pre-commit
[black]: https://github.com/psf/black

## Features

- Use the ATMO equilibrium chemistry plugin in your TauREx 3.1 forward models and retrievals.

## Requirements

- Ensure that you have a working installation of ATMO on your machine (https://www.gitlab.erc-atmo.eu/1-2d_tools/atmo). If you need access to ATMO, please contact Pascal Tremblin: pascal.tremblin@cea.fr
- You will need a working installation of TauREx 3.1 on your machine, with associated required Python packages: https://taurex3-public.readthedocs.io/

## Installation

<!-- You can install _frogr_ via [pip] from [PyPI]:

```console
$ pip install frogr
```

Alternatively, clone this Github and install via terminal with. This is done by: -->

You can install _frogr_ by cloning this Github and installing via the terminal. This is done by:

Cloning the directory using:

```console
$ git clone https://github.com/michellebieger/frogr.git
```

Move into frogr's folder:

```console
$ cd frogr
```

Install by then typing in:

```console
$ pip install .
```

If you are unable to install, a common error for HPC systems can be the Poetry log install requirements, which dictate a high-level version of Python. Try creating a new conda environment and then try pip installing once more:

```console
$ conda create -n [insertpreferredenvname] python=3.8
```

You can check the installation by importing frogr into Python:

```console
$ python -c "import atmopy"
```

To check that TauREx 3.1 has correctly registered your plugin:

```console
$ taurex --plugins
```

If there are no errors, you have been successful!

## Usage

Please see the [Command-line Reference] for details.

## License

Distributed under the terms of the [GPL 3.0 license][license],
_frogr_ is free and open source software.

## Issues

If you encounter any problems, please email michellebieger@live.com with a detailed description of the issue.

<!-- please [file an issue] along with a detailed description. -->

## Credits

This project was generated from [@cjolowicz]'s [Hypermodern Python Cookiecutter] template.

[@cjolowicz]: https://github.com/cjolowicz
[pypi]: https://pypi.org/
[hypermodern python cookiecutter]: https://github.com/cjolowicz/cookiecutter-hypermodern-python

<!-- [file an issue]: https://github.com/mfbieger/frogr/issues -->

[pip]: https://pip.pypa.io/

<!-- github-only -->

[license]: https://github.com/michellebieger/frogr/blob/main/LICENSE
[contributor guide]: https://github.com/michellebieger/frogr/blob/main/CONTRIBUTING.md
[command-line reference]: https://frogr-from-gas-to-retrieval.readthedocs.io/en/latest/usage.html
