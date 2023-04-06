# Usage

![frogr Logo](img/frogr_logo_reverse_nobckgd.png "frogr Logo")
<br>

_Above: The logo for *frogr*, designed with the PhyloPic CC0 1.0 Universal Public Domain image of **bufo bufo (Linnaeus 1758)**. The name *frogr* comes as a nod to the 1997 Hasbro Interactive game **Frogger: He's Back!**._

The following documentation assumes the user is familiar with TauREx 3.1, and has marginal familiarity with basic Python programming, the command line, and basic Anaconda Python editors. Documentation for TauREx 3.1 can be found [here](https://taurex3-public.readthedocs.io/). For any issues with the use of _frogr_, please contact michellebieger@live.com.

There are two ways in which you can use _frogr_: either through passing a .par file through to TauREx 3.1, or with Python with a preferred code editor such as Jupyter Notebooks or the command line.

## Par file use

To use _frogr_ when running a retrieval or forward model on command line, in a HPC server or similar, you can use the example .par file included in the Git repository of _frogr_, under the folder called _input_examples_. Note that when requesting to use _frogr_ as your chemistry model, the following keywords are equally acceptable: "ATMO", "atmo","frogr".

## Using with Python in a code editor

If you have Git cloned the _frogr_ directory, you will find a folder called _input_examples_, which provides a Jupyter Notebook tutorial in which _frogr_ is first used as a forward model for an example hot Jupiter exoplanet, WASP-79b. This Notebook then subsequently demonstrates how to incorporate _frogr_ in a trial retrieval, which can be run on a laptop. This syntax shown in the Notebook demonstrates how _frogr_ can be used in Python.

## What parameters can _frogr_ contribute to TauREx 3.1?

### Chemical equilibrium

The following parameters can be retrieved using ATMO:

| Retrievable parameter | Units                   |
| --------------------- | ----------------------- |
| Condensation of water | Title                   |
| Metallicity           | Solar metallicity units |

ATMO has the following molecules for which volume-mixing ratios can be calculated in the atmosphere: H, He, C, No, O, Na, K, Si, Ar, Ti, V, S, Cl, Mg, Al, Ca, Fe, Cr, Li, Cs, Rb, F, P.

### Chemical disequilibrium

The chemical disequilibrium code of ATMO cannot as yet be used with TauREx; however construction is ongoing and the next version of _frogr_ will likely implement this functionality.

## Reminder: please cite properly when using _frogr_

Please cite the following if using _frogr_ in your publication:

- Bieger et al. 2023 (in prep)
- The citations associated with [TauRex](https://taurex3-public.readthedocs.io/en/latest/citation.html)
- The citations associated with [ATMO](https://www.wiki.erc-atmo.eu/index.php?title=Credit) (link requires access to ERC ATMO wiki)

```{eval-rst}
.. click:: frogr.__main__:main
    :prog: frogr
    :nested: full
```
