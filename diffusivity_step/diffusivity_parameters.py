# -*- coding: utf-8 -*-
"""
Control parameters for the Diffusivity step in a SEAMM flowchart
"""

import logging
import seamm

logger = logging.getLogger(__name__)


class DiffusivityParameters(seamm.Parameters):
    """
    The control parameters for Diffusivity.

    You need to replace the "time" entry in dictionary below these comments with the
    definitions of parameters to control this step. The keys are parameters for the
    current plugin,the values are dictionaries as outlined below.

    Examples
    --------
    ::

        parameters = {
            "time": {
                "default": 100.0,
                "kind": "float",
                "default_units": "ps",
                "enumeration": tuple(),
                "format_string": ".1f",
                "description": "Simulation time:",
                "help_text": ("The time to simulate in the dynamics run.")
            },
        }

    parameters : {str: {str: str}}
        A dictionary containing the parameters for the current step.
        Each key of the dictionary is a dictionary that contains the
        the following keys:

    parameters["default"] :
        The default value of the parameter, used to reset it.

    parameters["kind"] : enum()
        Specifies the kind of a variable. One of  "integer", "float", "string",
        "boolean", or "enum"

        While the "kind" of a variable might be a numeric value, it may still have
        enumerated custom values meaningful to the user. For instance, if the parameter
        is a convergence criterion for an optimizer, custom values like "normal",
        "precise", etc, might be adequate. In addition, any parameter can be set to a
        variable of expression, indicated by having "$" as the first character in the
        field. For example, $OPTIMIZER_CONV.

    parameters["default_units"] : str
        The default units, used for resetting the value.

    parameters["enumeration"] : tuple
        A tuple of enumerated values.

    parameters["format_string"] : str
        A format string for "pretty" output.

    parameters["description"] : str
        A short string used as a prompt in the GUI.

    parameters["help_text"] : str
        A longer string to display as help for the user.

    See Also
    --------
    Diffusivity, TkDiffusivity, Diffusivity DiffusivityParameters, DiffusivityStep
    """

    parameters = {
        "approach": {
            "default": "both",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "both",
                "Helfand moments",
                "Mean Square Displacement (MSD)",
            ),
            "format_string": "",
            "description": "Approach:",
            "help_text": "The approach or method for determining the diffusivity.",
        },
        "nruns": {
            "default": "5",
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of runs to average:",
            "help_text": "The number for separate runs to average.",
        },
        "hydrodynamic correction": {
            "default": "no",
            "kind": "boolean",
            "format_string": "s",
            "enumeration": ("no", "yes"),
            "description": "Correction for cell size:",
            "help_text": (
                "Apply the Yeh-Hummer hydrodynamic correction for the cell size."
            ),
        },
        "viscosity": {
            "default": "0.0",
            "kind": "float",
            "default_units": "cP",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Viscosity:",
            "help_text": "The viscosity for the Yeh-Hummer cell size correction.",
        },
        "msd_fit_start": {
            "default": "0.1",
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Start point of MSD fit (0.0..1.0):",
            "help_text": "Starting point of fit the MSD this far into the data.",
        },
        "msd_fit_end": {
            "default": "0.6",
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "End point of MSD fit (0.0..1.0):",
            "help_text": "End point of fit the MSD this far into the data.",
        },
        "helfand_fit_start": {
            "default": "0.35",
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Start point of Helfand fit (0.0..1.0):",
            "help_text": "Starting point of Helfand fit this far into the data.",
        },
        "helfand_fit_end": {
            "default": "1.00",
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "End point of Helfand fit (0.0..1.0):",
            "help_text": "End point of Helfand fit this far into the data.",
        },
        "maximum Helfand Integral length": {
            "default": 1000,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Maximum Helfand Integral length:",
            "help_text": "Maximum length of the Integral used in the Helfand moments.",
        },
        "errors": {
            "default": "continue to next run",
            "kind": "string",
            "default_units": "",
            "enumeration": (
                "continue to next run",
                "exit the thermal conductivity step",
                "stop the job",
            ),
            "format_string": "s",
            "description": "On errors",
            "help_text": "How to handle errors in the runs",
        },
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": None,
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": "The results to save to variables or in tables.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("DiffusivityParameters.__init__")

        super().__init__(
            defaults={**DiffusivityParameters.parameters, **defaults}, data=data
        )
