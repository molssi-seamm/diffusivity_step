# -*- coding: utf-8 -*-

"""This file contains metadata describing the results from Diffusivity
"""

metadata = {}

"""Description of the computational models for Diffusivity.

Hamiltonians, approximations, and basis set or parameterizations,
only if appropriate for this code. For example::

    metadata["computational models"] = {
        "Hartree-Fock": {
            "models": {
                "PM7": {
                    "parameterizations": {
                        "PM7": {
                            "elements": "1-60,62-83",
                            "periodic": True,
                            "reactions": True,
                            "optimization": True,
                            "code": "mopac",
                        },
                        "PM7-TS": {
                            "elements": "1-60,62-83",
                            "periodic": True,
                            "reactions": True,
                            "optimization": False,
                            "code": "mopac",
                        },
                    },
                },
            },
        },
    }
"""
# metadata["computational models"] = {
# }

"""Description of the Diffusivity keywords.

(Only needed if this code uses keywords)

Fields
------
description : str
    A human readable description of the keyword.
takes values : int (optional)
    Number of values the keyword takes. If missing the keyword takes no values.
default : str (optional)
    The default value(s) if the keyword takes values.
format : str (optional)
    How the keyword is formatted in the MOPAC input.

For example::
    metadata["keywords"] = {
        "0SCF": {
            "description": "Read in data, then stop",
        },
        "ALT_A": {
            "description": "In PDB files with alternative atoms, select atoms A",
            "takes values": 1,
            "default": "A",
            "format": "{}={}",
        },
    }
"""
# metadata["keywords"] = {
# }

"""Properties that Diffusivity produces.
`metadata["results"]` describes the results that this step can produce. It is a
dictionary where the keys are the internal names of the results within this step, and
the values are a dictionary describing the result. For example::

    metadata["results"] = {
        "total_energy": {
            "calculation": [
                "energy",
                "optimization",
            ],
            "description": "The total energy",
            "dimensionality": "scalar",
            "methods": [
                "ccsd",
                "ccsd(t)",
                "dft",
                "hf",
            ],
            "property": "total energy#Psi4#{model}",
            "type": "float",
            "units": "E_h",
        },
    }

Fields
______

calculation : [str]
    Optional metadata describing what subtype of the step produces this result.
    The subtypes are completely arbitrary, but often they are types of calculations
    which is why this is name `calculation`. To use this, the step or a substep
    define `self._calculation` as a value. That value is used to select only the
    results with that value in this field.

description : str
    A human-readable description of the result.

dimensionality : str
    The dimensions of the data. The value can be "scalar" or an array definition
    of the form "[dim1, dim2,...]". Symmetric tringular matrices are denoted
    "triangular[n,n]". The dimensions can be integers, other scalar
    results, or standard parameters such as `n_atoms`. For example, '[3]',
    [3, n_atoms], or "triangular[n_aos, n_aos]".

methods : str
    Optional metadata like the `calculation` data. `methods` provides a second
    level of filtering, often used for the Hamiltionian for *ab initio* calculations
    where some properties may or may not be calculated depending on the type of
    theory.

property : str
    An optional definition of the property for storing this result. Must be one of
    the standard properties defined either in SEAMM or in this steps property
    metadata in `data/properties.csv`.

type : str
    The type of the data: string, integer, or float.

units : str
    Optional units for the result. If present, the value should be in these units.
"""
metadata["results"] = {
    # Control parameters
    "approach": {
        "description": "Algorithm for diffusion calculation",
        "dimensionality": "scalar",
        "type": "string",
    },
    "nruns": {
        "description": "Number of runs to average",
        "dimensionality": "scalar",
        "type": "integer",
    },
    "msd_fit_start": {
        "description": "Where to start the fit of the MSD curve",
        "dimensionality": "scalar",
        "type": "float",
    },
    "msd_fit_end": {
        "description": "Where to end the fit of the MSD curve",
        "dimensionality": "scalar",
        "type": "float",
    },
    "msd samples": {
        "description": "Number of samples to take for the MSD calculation",
        "dimensionality": "scalar",
        "type": "float",
    },
    "helfand_fit_start": {
        "description": "Where to start the fit of the Helfand moments curve",
        "dimensionality": "scalar",
        "type": "float",
    },
    "helfand_fit_end": {
        "description": "Where to end the fit of the Helfand moments curve",
        "dimensionality": "scalar",
        "type": "float",
    },
    "Helfand integral length": {
        "description": "The length of the Helfand numerical integration",
        "dimensionality": "scalar",
        "type": "integer",
    },
    # Computed results
    "T": {
        "description": "The temperature of the diffusion calculation",
        "dimensionality": "scalar",
        "type": "float",
        "units": "K",
    },
    "T,stderr": {
        "description": (
            "The standard error of the temperature of the diffusion calculation"
        ),
        "dimensionality": "scalar",
        "type": "float",
        "units": "K",
    },
    "P": {
        "description": "The pressure of the diffusion calculation",
        "dimensionality": "scalar",
        "type": "float",
        "units": "atm",
    },
    "P,stderr": {
        "description": (
            "The standard error of the pressure of the diffusion calculation"
        ),
        "dimensionality": "scalar",
        "type": "float",
        "units": "atm",
    },
    "density": {
        "description": "The density of the diffusion calculation",
        "dimensionality": "scalar",
        "type": "float",
        "units": "g/mL",
    },
    "density,stderr": {
        "description": (
            "The standard error of the density of the diffusion calculation"
        ),
        "dimensionality": "scalar",
        "type": "float",
        "units": "g/mL",
    },
    "1/L": {
        "description": "1/length of cell in the diffusion calculation",
        "dimensionality": "scalar",
        "type": "float",
        "units": "1/Å",
    },
    "1/L,stderr": {
        "description": (
            "Standard error if 1/length of cell in the diffusion calculation"
        ),
        "dimensionality": "scalar",
        "type": "float",
        "units": "1/Å",
    },
    "D {key} (MSD)": {
        "description": "The total diffusion coefficient from MSD",
        "dimensionality": "{species: value}",
        "property": "D#MSD#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "D {key} (MSD),stderr": {
        "description": "The standard error of the total diffusion coefficient from MSD",
        "dimensionality": "{species: value}",
        "property": "D,stderr#MSD#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dx {key} (MSD)": {
        "description": "The diffusion coefficient in x from MSD",
        "dimensionality": "{species: value}",
        "property": "Dx#MSD#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dx {key} (MSD),stderr": {
        "description": "The standard error if thediffusion coefficient in x from MSD",
        "dimensionality": "{species: value}",
        "property": "Dx,stderr#MSD#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dy {key} (MSD)": {
        "description": "The diffusion coefficient in y from MSD",
        "dimensionality": "{species: value}",
        "property": "Dy#MSD#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dy {key} (MSD),stderr": {
        "description": "The standard error if thediffusion coefficient in y from MSD",
        "dimensionality": "{species: value}",
        "property": "Dy,stderr#MSD#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dz {key} (MSD)": {
        "description": "The diffusion coefficient in z from MSD",
        "dimensionality": "{species: value}",
        "property": "Dz#MSD#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dz {key} (MSD),stderr": {
        "description": "The standard error if thediffusion coefficient in z from MSD",
        "dimensionality": "{species: value}",
        "property": "Dz,stderr#MSD#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "D {key} (HM)": {
        "description": "The total diffusion coefficient from Helfand Moments",
        "dimensionality": "{species: value}",
        "property": "D#Helfand Moments#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "D {key} (HM),stderr": {
        "description": (
            "The standard error of the total diffusion coefficient from Helfand Moments"
        ),
        "dimensionality": "{species: value}",
        "property": "D,stderr#Helfand Moments#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dx {key} (HM)": {
        "description": "The diffusion coefficient in x from Helfand Moments",
        "dimensionality": "{species: value}",
        "property": "Dx#Helfand Moments#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dx {key} (HM),stderr": {
        "description": (
            "The standard error of the diffusion coefficient in x from Helfand Moments"
        ),
        "property": "Dx,stderr#Helfand Moments#{model}",
        "dimensionality": "{species: value}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dy {key} (HM)": {
        "description": "The diffusion coefficient in y from Helfand Moments",
        "dimensionality": "{species: value}",
        "property": "Dy#Helfand Moments#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dy {key} (HM),stderr": {
        "description": (
            "The standard error of the diffusion coefficient in y from Helfand Moments"
        ),
        "property": "Dy,stderr#Helfand Moments#{model}",
        "dimensionality": "{species: value}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dz {key} (HM)": {
        "description": "The diffusion coefficient in z from Helfand Moments",
        "dimensionality": "{species: value}",
        "property": "Dz#Helfand Moments#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dz {key} (HM),stderr": {
        "description": (
            "The standard error of the diffusion coefficient in z from Helfand Moments"
        ),
        "property": "Dz,stderr#Helfand Moments#{model}",
        "dimensionality": "{species: value}",
        "type": "float",
        "units": "m^2/s",
    },
    "D {key}": {
        "description": "The total diffusion coefficient",
        "dimensionality": "{species: value}",
        "property": "D#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "D {key},stderr": {
        "description": ("The standard error of the total diffusion coefficient"),
        "dimensionality": "{species: value}",
        "property": "D,stderr#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dx {key}": {
        "description": "The diffusion coefficient in x",
        "dimensionality": "{species: value}",
        "property": "Dx#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dx {key},stderr": {
        "description": ("The standard error of the diffusion coefficient in x"),
        "property": "Dx,stderr#{model}",
        "dimensionality": "{species: value}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dy {key}": {
        "description": "The diffusion coefficient in y",
        "dimensionality": "{species: value}",
        "property": "Dy#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dy {key},stderr": {
        "description": ("The standard error of the diffusion coefficient in y"),
        "property": "Dy,stderr#{model}",
        "dimensionality": "{species: value}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dz {key}": {
        "description": "The diffusion coefficient in z",
        "dimensionality": "{species: value}",
        "property": "Dz#{model}",
        "type": "float",
        "units": "m^2/s",
    },
    "Dz {key},stderr": {
        "description": ("The standard error of the diffusion coefficient in z"),
        "property": "Dz,stderr#{model}",
        "dimensionality": "{species: value}",
        "type": "float",
        "units": "m^2/s",
    },
    # Timings
    "t_msd": {
        "description": "The time taken in the MSD analysis",
        "dimensionality": "scalar",
        "type": "float",
        "units": "s",
    },
    "t_hm": {
        "description": "The time taken in the Helfand moments analysis",
        "dimensionality": "scalar",
        "type": "float",
        "units": "s",
    },
}
