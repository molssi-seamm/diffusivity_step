# -*- coding: utf-8 -*-

"""Non-graphical part of the Diffusivity step in a SEAMM flowchart
"""

import logging
import math
from math import log10, ceil, floor, sqrt
from pathlib import Path
import pkg_resources
import sys
import time
import traceback

import numpy as np
from tabulate import tabulate

from .analysis import (
    read_vector_trajectory,
    compute_msd,
    add_msd_trace,
    fit_slope,
    add_helfand_trace,
    create_helfand_moments,
)
import diffusivity_step
import molsystem
import seamm
import seamm_util
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Diffusivity")

# Add this module's properties to the standard properties
path = Path(pkg_resources.resource_filename(__name__, "data/"))
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


def fmt_err(value, err, precision=2, as_float=False):
    try:
        decimals = -ceil(log10(err)) + precision
    except Exception:
        e = "--"
        try:
            v = f"{value:.3f}"
            if as_float:
                v = round(float(v), 3)
        except Exception:
            v = value
    else:
        if decimals < 0:
            decimals = 0
        if as_float:
            e = round(err, decimals)
        else:
            fmt = f".{decimals}f"
            e = f"{err:{fmt}}"
        try:
            if as_float:
                if abs(value / e) > 1.0e6:
                    decimals = 6
                v = round(float(value), decimals)
                e = round(float(e), decimals)
            else:
                v = f"{value:{fmt}}"
        except Exception:
            v = value
    return v, e


class Diffusivity(seamm.Node):
    """
    The non-graphical part of a Diffusivity step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : DiffusivityParameters
        The control parameters for Diffusivity.

    See Also
    --------
    TkDiffusivity,
    Diffusivity, DiffusivityParameters
    """

    def __init__(
        self,
        flowchart=None,
        title="Diffusivity",
        namespace="org.molssi.seamm",
        extension=None,
        logger=logger,
    ):
        """A step for Diffusivity in a SEAMM flowchart.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        namespace : str
            The namespace for the plug-ins of the subflowchart
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Diffusivity {self}")

        self.subflowchart = seamm.Flowchart(
            parent=self, name="Diffusivity", namespace=namespace
        )

        super().__init__(
            flowchart=flowchart,
            title="Diffusivity",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = diffusivity_step.metadata
        self.parameters = diffusivity_step.DiffusivityParameters()
        self._file_handler = None

        # Set our citation level to 1!
        self.citation_level = 1

        self._tensor_labels = [
            ("x", "red", "rgba(255,0,0,0.1)"),
            ("y", "green", "rgba(0,255,0,0.1)"),
            ("z", "blue", "rgba(0,0,255,0.1)"),
            ("", "black", "rgba(0,0,0,0.1)"),
        ]

        self._configuration = None  # The initial configuration
        self._n_molecules = None
        self._species = None  # The molecule indices for each species

        self._state_vars = {}  # The state of the system in the MD
        self._results = {}  # The results data.

        self._use_msd = False
        self._msd_dt = None
        self._msds = None
        self._msd_errs = None
        self._msd = None
        self._msd_err = None
        self._scale = None  # Factor to make the results nice numbers
        self._msd_data = {}

        self._use_velocity = False
        self._velocity_dt = None
        self._Ms = None
        self._M_errs = None
        self._M = None
        self._M_err = None
        self._M_coeff = []
        self._M_stderr = []
        self._M_data = {}

        self._state = {}

        # Internal statistics
        self._msd_samples = None
        self._msd_times = []
        self._helfand_integral_length = None
        self._helfand_integral_times = []

    @property
    def version(self):
        """The semantic version of this module."""
        return diffusivity_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return diffusivity_step.__git_revision__

    @property
    def configuration(self):
        """The initial configuration for this run."""
        if self._configuration is None:
            _, self._configuration = self.get_system_configuration
        return self._configuration

    @property
    def n_molecules(self):
        """The number of molecules in the system."""
        return self._n_molecules

    @property
    def results(self):
        """The results data."""
        return self._results

    @property
    def species(self):
        """The species and list of molecules for each."""
        if self._species is None:
            self._species = self.parse_molecules()
        return self._species

    def analyze_run(
        self,
        indent="",
        P=None,
        run=None,
        weighted=False,
        **kwargs,
    ):
        """Do any analysis of the output from this step.

        Also print important results to the local step.out file using
        "printer".

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        weighted : bool
            Whether to use the stderr to weight the fit, default False
        """
        # Print useful information about the analysis cost before the first run results
        if run == 1:
            if self._msd_samples is not None:
                self.results["msd samples"] = self._msd_samples
                text = (
                    f"The MSD analysis over {self._msd_samples} steps took "
                    f"{self._msd_times[0]:.3f} seconds."
                )
                printer.normal(__(text, indent=8 * " "))
            if self._helfand_integral_length is not None:
                self.results["Helfand integral length"] = self._helfand_integral_length
                text = (
                    f"The Helfand moment analysis over {self._helfand_integral_length} "
                    f"steps took {self._helfand_integral_times[0]:.3f} seconds."
                )
                printer.normal(__(text, indent=8 * " "))
            printer.normal("")

        # If requested, calculate the Yeh-Hummer correction for cell size
        correction = None
        if P["hydrodynamic correction"]:
            if "T" in self._state_vars:
                T = self._state["T"][-1]
                if "a" in self._state_vars:
                    L = self._state["a"][-1]
                else:
                    L = self._configuration.cell.a
                L = Q_(L, "Å")

                correction = (
                    2.837297
                    * Q_(1.0, "k_B")
                    * Q_(T, "K")
                    / (6 * math.pi * P["viscosity"] * L)
                )
                correction = float(correction.m_as("m^2/s"))
            else:
                if run == 1:
                    text = (
                        "\nThe temperature is not available, so cannot calculate the "
                        "Yeh-Hummer correction. You need to save it to the variable "
                        "'T' in the molecular dynamics."
                    )
                printer.normal(__(text, indent=8 * " "))

        table = {
            "Run": [],
            "Method": [],
            "Species": [],
            "Dx": [],
            "ex": [],
            "Dy": [],
            "ey": [],
            "Dz": [],
            "ez": [],
            "D": [],
            "e": [],
        }
        if correction is not None:
            table["Correction"] = []

        # Add state variables at end
        for var in self._state_vars:
            var = var.replace(",stderr", " ±")
            table[var] = []

        # On to the results
        if self._use_msd:
            # Fit the MSD to a line
            nframes, nalpha = self._msd[0].shape
            ts = np.arange(nframes) * self._msd_dt.m_as("ps")  # Scale to ps
            ts = ts.tolist()

            # Set up storage for each run the first time
            msd_data = self._msd_data
            if run == 1:
                # Space for the results from each run. Arrays indexed by smiles & alpha
                per_run = msd_data["per_run"] = {}
                for smiles in self.species.keys():
                    per_run[smiles] = {}
                    for i in range(nalpha):
                        per_run[smiles][i] = []

                # Space for the mean and stderr so far
                mean = msd_data["mean"] = {}
                stderr = msd_data["stderr"] = {}
                for smiles in self.species.keys():
                    mean[smiles] = {}
                    stderr[smiles] = {}
                    for i in range(nalpha):
                        mean[smiles][i] = []
                        stderr[smiles][i] = []
                    mean["D"] = []
                    stderr["D"] = []

            # Create the plot for the MSD
            figure = self.create_figure(
                module_path=("seamm",),
                template="line.graph_template",
                title="Mean Square Displacement",
            )
            plot = figure.add_plot("MSD")

            x_axis = plot.add_axis("x", label="Time (ps)")
            y_axis = plot.add_axis("y", label="MSD (Å^2)", anchor=x_axis)
            x_axis.anchor = y_axis

            for spec, smiles in enumerate(self.species.keys()):
                coeffs = msd_data["per_run"][smiles]
                # Fit the slopes
                fit = []
                # Convert units and remember the factor of 1/6 in the Einstein equation
                factor = Q_("Å^2/ps").m_as("m^2/s") / 6
                for i in range(nalpha):
                    if weighted:
                        slope, err, xs, ys = fit_slope(
                            self._msd[spec][:, i],
                            ts,
                            sigma=self._msd_err[spec][:, i],
                            start=P["msd_fit_start"],
                            end=P["msd_fit_end"],
                        )
                    else:
                        slope, err, xs, ys = fit_slope(
                            self._msd[spec][:, i],
                            ts,
                            start=P["msd_fit_start"],
                            end=P["msd_fit_end"],
                        )
                    d_coeff = slope * factor
                    err = err * factor

                    # Set a scale factor to make the numbers managable
                    if self._scale is None:
                        self._scale = 10 ** floor(log10(d_coeff))

                    if correction is not None:
                        if i < 3:
                            d_coeff += correction / 3
                        else:
                            d_coeff += correction

                    msd_data["per_run"][smiles][i].append(d_coeff)

                    if run == 1:
                        v = round(d_coeff / self._scale, 3)
                        msd_data["mean"][smiles][i].append(v)
                        msd_data["stderr"][smiles][i].append(None)
                    else:
                        coeffs = msd_data["per_run"][smiles][i]
                        tmp = np.array(coeffs)
                        mean = float(tmp.mean())
                        std = float(tmp.std())
                        stderr = std / sqrt(run - 1)
                        v_mean, e_mean = fmt_err(
                            mean / self._scale, 2 * stderr / self._scale, as_float=True
                        )
                        v, e = fmt_err(
                            d_coeff / self._scale,
                            2 * stderr / self._scale,
                            as_float=True,
                        )
                        msd_data["mean"][smiles][i].append(v_mean)
                        msd_data["stderr"][smiles][i].append(e_mean)

                    fit.append(
                        {
                            "Species": smiles,
                            "D": d_coeff,
                            "Dmean": mean,
                            "stderr": stderr,
                            "xs": xs,
                            "ys": ys,
                            "scale": self._scale,
                            "D_s": v,
                        }
                    )

                    if i == 0:
                        if spec == 0:
                            table["Run"].append(run)
                            if correction is not None:
                                table["Correction"].append(
                                    f"{correction / self._scale:.3f}"
                                )
                            for var in self._state_vars:
                                value = self._state[var][-1]
                                if value != "":
                                    if "density" in var:
                                        value = round(value, 3)
                                    else:
                                        value = round(value, 1)
                                var = var.replace(",stderr", " ±")
                                table[var].append(value)
                            table["Method"].append("MSD")
                        else:
                            table["Run"].append("")
                            if correction is not None:
                                table["Correction"].append("")
                            table["Method"].append("")
                        table["Species"].append(smiles)
                    alpha = self._tensor_labels[i][0]
                    table["D" + alpha].append(v)
                    table["e" + alpha].append("")

                add_msd_trace(
                    plot,
                    x_axis,
                    y_axis,
                    smiles,
                    self._msd[spec],
                    ts,
                    fit=fit,
                )

            figure.grid_plots("MSD")

            # Write to disk
            filename = "MSD.graph"
            path = Path(self.directory) / filename
            figure.dump(path)

            if "html" in self.options and self.options["html"]:
                path = path.with_suffix(".html")
                figure.template = "line.html_template"
                figure.dump(path)

        if self._use_velocity:
            # Fit the Helfand moments
            nframes, nalpha = self._M[0].shape
            ts = np.arange(nframes) * self._velocity_dt.m_as("ps")  # Scale to ps
            ts = ts.tolist()

            # Set up storage for each run the first time
            M_data = self._M_data
            if run == 1:
                # Space for the results from each run. Arrays indexed by smiles & alpha
                per_run = M_data["per_run"] = {}
                for smiles in self.species.keys():
                    per_run[smiles] = {}
                    for i in range(nalpha):
                        per_run[smiles][i] = []

                # Space for the mean and stderr so far
                mean = M_data["mean"] = {}
                stderr = M_data["stderr"] = {}
                for smiles in self.species.keys():
                    mean[smiles] = {}
                    stderr[smiles] = {}
                    for i in range(nalpha):
                        mean[smiles][i] = []
                        stderr[smiles][i] = []
                    mean["D"] = []
                    stderr["D"] = []
                if self._use_msd:
                    mean = M_data["combined mean"] = {}
                    stderr = M_data["combined stderr"] = {}
                    for smiles in self.species.keys():
                        mean[smiles] = {}
                        stderr[smiles] = {}
                        for i in range(nalpha):
                            mean[smiles][i] = []
                            stderr[smiles][i] = []
                        mean["D"] = []
                        stderr["D"] = []

            # Create the plot for the Helfand moments
            figure = self.create_figure(
                module_path=("seamm",),
                template="line.graph_template",
                title="Helfand Moments",
            )
            plot = figure.add_plot("HelfandMoments")

            x_axis = plot.add_axis("x", label="Time (ps)")
            y_axis = plot.add_axis("y", label="M (Å^2)", anchor=x_axis)
            x_axis.anchor = y_axis

            for spec, smiles in enumerate(self.species.keys()):
                coeffs = M_data["per_run"][smiles]
                # Fit the slopes
                fit = []
                # Convert units and remember the factor of 1/3 in the Einstein equation
                factor = Q_("Å^2/ps").m_as("m^2/s") / 3
                for i in range(nalpha):
                    if weighted:
                        slope, err, xs, ys = fit_slope(
                            self._M[spec][:, i],
                            ts,
                            sigma=self._M_err[spec][:, i],
                            start=P["helfand_fit_start"],
                            end=P["helfand_fit_end"],
                        )
                    else:
                        slope, err, xs, ys = fit_slope(
                            self._M[spec][:, i],
                            ts,
                            start=P["helfand_fit_start"],
                            end=P["helfand_fit_end"],
                        )
                    d_coeff = slope * factor
                    err = err * factor

                    # Set a scale factor to make the numbers managable
                    if self._scale is None:
                        self._scale = 10 ** floor(log10(d_coeff))

                    if correction is not None:
                        if i < 3:
                            d_coeff += correction / 3
                        else:
                            d_coeff += correction

                    M_data["per_run"][smiles][i].append(d_coeff)

                    if run == 1:
                        v = f"{d_coeff / self._scale:.2f}"
                        M_data["mean"][smiles][i].append(float(v))
                        M_data["stderr"][smiles][i].append(None)
                    else:
                        coeffs = M_data["per_run"][smiles][i]
                        tmp = np.array(coeffs)
                        mean = float(tmp.mean())
                        std = float(tmp.std())
                        stderr = std / sqrt(run - 1)
                        v_mean, e_mean = fmt_err(
                            mean / self._scale, 2 * stderr / self._scale, as_float=True
                        )

                        v, e = fmt_err(
                            d_coeff / self._scale,
                            2 * stderr / self._scale,
                            as_float=True,
                        )
                        M_data["mean"][smiles][i].append(v_mean)
                        M_data["stderr"][smiles][i].append(e_mean)

                    if self._use_msd:
                        tmp = np.array(
                            M_data["per_run"][smiles][i]
                            + msd_data["per_run"][smiles][i]
                        )
                        mean = float(tmp.mean())
                        std = float(tmp.std())
                        stderr = std / sqrt(2 * run - 1)
                        v_mean, e_mean = fmt_err(
                            mean / self._scale,
                            2 * stderr / self._scale,
                            as_float=True,
                        )

                        M_data["combined mean"][smiles][i].append(v_mean)
                        M_data["combined stderr"][smiles][i].append(e_mean)

                    fit.append(
                        {
                            "Species": smiles,
                            "D": d_coeff,
                            "stderr": err,
                            "xs": xs,
                            "ys": ys,
                            "scale": self._scale,
                            "D_s": v,
                        }
                    )

                    if i == 0:
                        if spec == 0:
                            if not self._use_msd:
                                table["Run"].append(run)
                            else:
                                table["Run"].append("")
                            for var in self._state_vars:
                                value = self._state[var][-1]
                                if value != "":
                                    if "density" in var:
                                        value = round(value, 3)
                                    else:
                                        value = round(value, 1)
                                var = var.replace(",stderr", " ±")
                                table[var].append(value)
                            table["Method"].append("Helfand Moments")
                        else:
                            table["Run"].append("")
                            table["Method"].append("")
                        table["Species"].append(smiles)
                    alpha = self._tensor_labels[i][0]
                    table["D" + alpha].append(v)
                    table["e" + alpha].append("")

                add_helfand_trace(
                    plot,
                    x_axis,
                    y_axis,
                    smiles,
                    self._M[spec],
                    ts,
                    fit=fit,
                )

            figure.grid_plots("HelfandMoments")

            # Write to disk
            filename = "HelfandMoments.graph"
            path = Path(self.directory) / filename
            figure.dump(path)

            if "html" in self.options and self.options["html"]:
                path = path.with_suffix(".html")
                figure.template = "line.html_template"
                figure.dump(path)

        # Add the current best estimates
        if run > 1:
            for spec, smiles in enumerate(self.species.keys()):
                for i in range(nalpha):
                    if self._use_msd:
                        if i == 0:
                            if spec == 0:
                                table["Run"].append("")
                                for var in self._state_vars:
                                    if ",stderr" in var:
                                        continue
                                    tmp = np.array(self._state[var])
                                    t_mean = float(tmp.mean())
                                    std = float(tmp.std())
                                    t_stderr = std / sqrt(run - 1)
                                    t_v, t_e = fmt_err(
                                        t_mean, 2 * t_stderr, as_float=True
                                    )
                                    table[var].append(t_v)
                                    table[var + " ±"].append(t_e)
                                if self._use_velocity:
                                    table["Method"].append("Current MSD Estimate")
                                else:
                                    table["Method"].append("Current Estimate")
                            else:
                                table["Run"].append("")
                                table["Method"].append("")
                            table["Species"].append(smiles)

                        mean = msd_data["mean"][smiles][i][-1]
                        stderr = msd_data["stderr"][smiles][i][-1]
                        alpha = self._tensor_labels[i][0]
                        table["D" + alpha].append(mean)
                        table["e" + alpha].append(stderr)
                    if self._use_velocity:
                        if i == 0:
                            if spec == 0:
                                table["Run"].append("")
                                for var in self._state_vars:
                                    if ",stderr" in var:
                                        continue
                                    tmp = np.array(self._state[var])
                                    t_mean = float(tmp.mean())
                                    std = float(tmp.std())
                                    t_stderr = std / sqrt(run - 1)
                                    t_v, t_e = fmt_err(
                                        t_mean, 2 * t_stderr, as_float=True
                                    )
                                    table[var].append(t_v)
                                    table[var + " ±"].append(t_e)
                                if self._use_msd:
                                    table["Method"].append("Current HM Estimate")
                                else:
                                    table["Method"].append("Current Estimate")
                            else:
                                table["Run"].append("")
                                table["Method"].append("")
                            table["Species"].append(smiles)

                        mean = M_data["mean"][smiles][i][-1]
                        stderr = M_data["stderr"][smiles][i][-1]
                        alpha = self._tensor_labels[i][0]
                        table["D" + alpha].append(mean)
                        table["e" + alpha].append(stderr)
                    if self._use_msd and self._use_velocity:
                        if i == 0:
                            if spec == 0:
                                table["Run"].append("")
                                for var in self._state_vars:
                                    if ",stderr" in var:
                                        continue
                                    tmp = np.array(self._state[var])
                                    t_mean = float(tmp.mean())
                                    std = float(tmp.std())
                                    t_stderr = std / sqrt(run - 1)
                                    t_v, t_e = fmt_err(
                                        t_mean, 2 * t_stderr, as_float=True
                                    )
                                    table[var].append(t_v)
                                    table[var + " ±"].append(t_e)
                                table["Method"].append("Combined Estimate")
                            else:
                                table["Run"].append("")
                                table["Method"].append("")
                            table["Species"].append(smiles)

                        mean = M_data["combined mean"][smiles][i][-1]
                        stderr = M_data["combined stderr"][smiles][i][-1]
                        alpha = self._tensor_labels[i][0]
                        table["D" + alpha].append(mean)
                        table["e" + alpha].append(stderr)

        # Print the table of results
        text = ""
        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            disable_numparse=True,
        )
        if run == 1:
            length = len(tmp.splitlines()[0])
            text += "\n"
            text += f"Diffusion Coefficients (* {self._scale:.1e} m^2/s)".center(length)
            text += "\n"
            text += "\n".join(tmp.splitlines()[0:-1])
        else:
            first = 3
            if run is not None and run == P["nruns"]:
                text += "\n".join(tmp.splitlines()[first:])
            else:
                text += "\n".join(tmp.splitlines()[first:-1])

        printer.normal(__(text, indent=8 * " ", wrap=False, dedent=False))

    def analyze(
        self,
        indent="",
        P=None,
        run=None,
        weighted=False,
        **kwargs,
    ):
        """Do any analysis of the output from this run.

        Also print important results to the local step.out file using
        "printer".

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        weighted : bool
            Whether to use the stderr to weight the fit, default False
        """

        # Print the final table of results
        table = {
            "Species": [],
            "Method": [],
            "Dir": [],
            "D": [],
            "±": [],
            "95%": [],
        }
        # First, the state variables
        for var in self._state_vars:
            if ",stderr" in var:
                continue
            n = len(self._state[var])
            tmp = np.array(self._state[var])
            if n == 1 or np.allclose(tmp, tmp[0], atol=0.0):
                t_v = round(self._state[var][0], 3)
                t_e = ""
            else:
                t_mean = float(tmp.mean())
                std = float(tmp.std())
                t_stderr = std / sqrt(n - 1)
                t_v, t_e = fmt_err(t_mean, 2 * t_stderr, as_float=True)
            table["Species"].append("")
            table["Method"].append(var)
            table["Dir"].append("")
            table["D"].append(t_v)
            table["95%"].append(t_e)
            self.results[var] = t_v
            try:
                self.results[var + ",stderr"] = float(t_e)
            except ValueError:
                pass
            if var == "a":
                tmp = 1.0 / tmp
                if n == 1 or np.allclose(tmp, tmp[0], atol=0.0):
                    t_v = round(1 / self._state[var][0], 4)
                    t_e = ""
                else:
                    t_mean = float(tmp.mean())
                    std = float(tmp.std())
                    t_stderr = std / sqrt(n - 1)
                    t_v, t_e = fmt_err(t_mean, 2 * t_stderr, as_float=True)
                self.results["1/L"] = t_v
                try:
                    self.results["1/L,stderr"] = float(t_e)
                except ValueError:
                    pass
                table["Species"].append("")
                table["Method"].append("1/L")
                table["Dir"].append("")
                table["D"].append(t_v)
                table["95%"].append(t_e)

        # Round to this many decimals to avoid unpleasant roundoff
        decimals = -ceil(log10(self._scale)) + 6

        if self._use_msd:
            msd_data = self._msd_data
            for smiles in self.species.keys():
                for i in range(4):
                    v = msd_data["mean"][smiles][i][-1]
                    e = msd_data["stderr"][smiles][i][-1]
                    d_coeff = round(v * self._scale, decimals)
                    if e is not None:
                        err = round(e * self._scale, decimals)
                    table["Species"].append(smiles if i == 0 else "")
                    table["Method"].append("MSD" if i == 0 else "")
                    if self._tensor_labels[i][0] == "":
                        table["Dir"].append("total")
                        if "D {key} (MSD)" in self.results:
                            self.results["D {key} (MSD)"][smiles] = d_coeff
                            if e is not None:
                                self.results["D {key} (MSD),stderr"][smiles] = err
                        else:
                            self.results["D {key} (MSD)"] = {smiles: d_coeff}
                            if e is not None:
                                self.results["D {key} (MSD),stderr"] = {smiles: err}
                        if "D {key}" in self.results:
                            self.results["D {key}"][smiles] = d_coeff
                            if e is not None:
                                self.results["D {key},stderr"][smiles] = err
                        else:
                            self.results["D {key}"] = {smiles: d_coeff}
                            if e is not None:
                                self.results["D {key},stderr"] = {smiles: err}
                    else:
                        table["Dir"].append(self._tensor_labels[i][0])
                        key = f"D{self._tensor_labels[i][0]}" + " {key} (MSD)"
                        if key in self.results:
                            self.results[key][smiles] = d_coeff
                            if e is not None:
                                self.results[key + ",stderr"][smiles] = err
                        else:
                            self.results[key] = {smiles: d_coeff}
                            if e is not None:
                                self.results[key + ",stderr"] = {smiles: err}
                        key = f"D{self._tensor_labels[i][0]}" + " {key}"
                        if key in self.results:
                            self.results[key][smiles] = d_coeff
                            if e is not None:
                                self.results[key + ",stderr"][smiles] = err
                        else:
                            self.results[key] = {smiles: d_coeff}
                            if e is not None:
                                self.results[key + ",stderr"] = {smiles: err}
                    table["D"].append(v)
                    table["±"].append("±")
                    table["95%"].append(e)
        if self._use_velocity:
            M_data = self._M_data
            for smiles in self.species.keys():
                for i in range(4):
                    v = M_data["mean"][smiles][i][-1]
                    e = M_data["stderr"][smiles][i][-1]
                    d_coeff = round(v * self._scale, decimals)
                    if e is not None:
                        err = round(e * self._scale, decimals)
                    table["Species"].append(smiles if i == 0 else "")
                    table["Method"].append("Helfand Moments" if i == 0 else "")
                    if self._tensor_labels[i][0] == "":
                        table["Dir"].append("total")
                        if "D {key} (HM)" in self.results:
                            self.results["D {key} (HM)"][smiles] = d_coeff
                            if e is not None:
                                self.results["D {key} (HM),stderr"][smiles] = err
                        else:
                            self.results["D {key} (HM)"] = {smiles: d_coeff}
                            if e is not None:
                                self.results["D {key} (HM),stderr"] = {smiles: err}
                        if "D {key}" in self.results:
                            self.results["D {key}"][smiles] = d_coeff
                            if e is not None:
                                self.results["D {key},stderr"][smiles] = err
                        else:
                            self.results["D {key}"] = {smiles: d_coeff}
                            if e is not None:
                                self.results["D {key},stderr"] = {smiles: err}
                    else:
                        table["Dir"].append(self._tensor_labels[i][0])
                        key = f"D{self._tensor_labels[i][0]}" + " {key} (HM)"
                        if key in self.results:
                            self.results[key][smiles] = d_coeff
                            if e is not None:
                                self.results[key + ",stderr"][smiles] = err
                        else:
                            self.results[key] = {smiles: d_coeff}
                            if e is not None:
                                self.results[key + ",stderr"] = {smiles: err}
                        key = f"D{self._tensor_labels[i][0]}" + " {key}"
                        if key in self.results:
                            self.results[key][smiles] = d_coeff
                            if e is not None:
                                self.results[key + ",stderr"][smiles] = err
                        else:
                            self.results[key] = {smiles: d_coeff}
                            if e is not None:
                                self.results[key + ",stderr"] = {smiles: err}
                    table["D"].append(v)
                    table["±"].append("±")
                    table["95%"].append(e)
        if self._use_msd and self._use_velocity:
            M_data = self._M_data
            for smiles in self.species.keys():
                for i in range(4):
                    v = M_data["combined mean"][smiles][i][-1]
                    e = M_data["combined stderr"][smiles][i][-1]
                    d_coeff = round(v * self._scale, decimals)
                    if e is not None:
                        err = round(e * self._scale, decimals)
                    table["Species"].append(smiles if i == 0 else "")
                    table["Method"].append("Combined" if i == 0 else "")
                    if self._tensor_labels[i][0] == "":
                        table["Dir"].append("total")
                        if "D {key}" in self.results:
                            self.results["D {key}"][smiles] = d_coeff
                            if e is not None:
                                if "D {key},stderr" in self.results:
                                    self.results["D {key},stderr"][smiles] = err
                                else:
                                    self.results["D {key},stderr"] = {smiles: err}
                        else:
                            self.results["D {key}"] = {smiles: d_coeff}
                            if e is not None:
                                self.results["D {key},stderr"] = {smiles: err}
                    else:
                        table["Dir"].append(self._tensor_labels[i][0])
                        key = f"D{self._tensor_labels[i][0]}" + " {key}"
                        if key in self.results:
                            self.results[key][smiles] = d_coeff
                            if e is not None:
                                if key + ",stderr" in self.results:
                                    self.results[key + ",stderr"][smiles] = err
                                else:
                                    self.results[key + ",stderr"] = {smiles: err}
                        else:
                            self.results[key] = {smiles: d_coeff}
                            if e is not None:
                                self.results[key + ",stderr"] = {smiles: err}
                    table["D"].append(v)
                    table["±"].append("±")
                    table["95%"].append(e)

        text = ""
        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=(
                "center",
                "center",
                "center",
                "decimal",
                "decimal",
            ),
        )
        length = len(tmp.splitlines()[0])
        text += "\n"
        text += f"Diffusion Coefficients (* {self._scale:.1e} m^2/s)".center(length)
        text += "\n"
        text += tmp
        text += "\n"

        printer.normal(__(text, indent=8 * " ", wrap=False, dedent=False))

        # And store results, only for the full output at the end
        ff = self.get_variable("_forcefield")
        if ff == "OpenKIM":
            self._model = "OpenKIM/" + self.get_variable("_OpenKIM_Potential")
        else:
            # Valence forcefield...
            self._model = ff.current_forcefield

        self.store_results(
            configuration=self.configuration,
            data=self.results,
        )

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = "diffusivity-step"
        parser = seamm_util.getParser()

        # Remember if the parser exists ... this type of step may have been
        # found before
        parser_exists = parser.exists(parser_name)

        # Create the standard options, e.g. log-level
        super().create_parser(name=parser_name)

        if not parser_exists:
            # Any options for diffusivity itself
            parser.add_argument(
                parser_name,
                "--html",
                action="store_true",
                help="whether to write out html files for graphs, etc.",
            )

        # Now need to walk through the steps in the subflowchart...
        self.subflowchart.reset_visited()
        node = self.subflowchart.get_node("1").next()
        while node is not None:
            node = node.create_parser()

        return self.next()

    def description_text(self, P=None, short=False):
        """Create the text description of what this step will do.
        The dictionary of control values is passed in as P so that
        the code can test values, etc.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        Returns
        -------
        str
            A description of the current step.
        """
        if P is None:
            P = self.parameters.values_to_dict()

        if P["approach"] == "both":
            text = (
                "Calculate the diffusivity using both the Helfand moments method and "
                "the mean square displacement (MSD)"
            )
        else:
            text = f"Calculate the diffusivity using the {P['approach']} approach"
        text += f", averaging over {P['nruns']} runs."

        if "Helfand" not in P["approach"]:
            text += (
                "The slope of the MSD vs time  will be fit from fraction "
                f"{P['msd_fit_start']} to {P['msd_fit_end']} of the total "
                "curve."
            )
        if "MSD" not in P["approach"]:
            text += (
                "The numerical integral for the Helfand moments will use a maximum of "
                "{maximum Helfand Integral length} steps. The cost scales as the "
                "square of the number of steps in the integration, so be cautious "
                "with large values. "
                "The slope of the Helfand moment vs time  will be fit from fraction "
                f"{P['helfand_fit_start']} to {P['helfand_fit_end']} of "
                "the total curve."
            )

        correction = P["hydrodynamic correction"]
        if type(correction) is not bool and self.is_expr(correction):
            text += (
                f" The expression '{correction}' will determine whether to use the "
                "Yeh-Hummer hydrodynamic correction for the finite cell size using "
                f"a viscosity of {P['viscosity']}."
            )
        if type(correction) is bool or "yes" in correction:
            text += (
                " The calculated diffusivity will be corrected for the finite cell "
                "size using the Yeh-Hummer hydrodynamic correction using "
                f"a viscosity of {P['viscosity']}."
            )

        text += "\n\n"

        # Make sure the subflowchart has the data from the parent flowchart
        self.subflowchart.root_directory = self.flowchart.root_directory
        self.subflowchart.executor = self.flowchart.executor
        self.subflowchart.in_jobserver = self.subflowchart.in_jobserver

        if not short:
            # Get the first real node
            node = self.subflowchart.get_node("1").next()
            node.all_options = self.all_options

            while node is not None:
                try:
                    text += __(node.description_text()).__str__()
                except Exception as e:
                    print(f"Error describing diffusivity flowchart: {e} in {node}")
                    logger.critical(
                        f"Error describing diffusivity flowchart: {e} in {node}"
                    )
                    raise
                except:  # noqa: E722
                    print(
                        "Unexpected error describing diffusivity flowchart: "
                        f"{sys.exc_info()[0]} in {str(node)}"
                    )
                    logger.critical(
                        "Unexpected error describing diffusivity flowchart: "
                        f"{sys.exc_info()[0]} in {str(node)}"
                    )
                    raise
                text += "\n"
                node = node.next()

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def run(self):
        """Run the diffusivity step.

        Parameters
        ----------
        None

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        next_node = super().run(printer)

        # Get the values of the parameters, dereferencing any variables
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Remember the configuration
        _, self._configuration = self.get_system_configuration()

        # Reset parameters if called in e.g. loop
        self._n_molecules = None
        self._species = None

        self.state_vars = {}
        self._results = {}

        self._use_msd = False
        self._msd_dt = None
        self._msds = None
        self._msd_errs = None
        self._msd = None
        self._msd_err = None
        self._scale = None  # Factor to make the results nice numbers
        self._msd_data = {}

        self._use_velocity = False
        self._velocity_dt = None
        self._Ms = None
        self._M_errs = None
        self._M = None
        self._M_err = None
        self._M_coeff = []
        self._M_stderr = []
        self._M_data = {}

        self._state = {}

        self._msd_samples = None
        self._msd_times = []
        self._helfand_integral_length = None
        self._helfand_integral_times = []

        # Decide what to do
        approach = P["approach"].lower()
        self._use_msd = approach == "both" or "msd" in approach
        self._use_velocity = approach == "both" or "helfand" in approach

        # Print what we are doing
        printer.important(__(self.description_text(P, short=True), indent=self.indent))

        # Transfer key parameters to results
        self._results = {}
        for key, value in P.items():
            self.results[key] = value

        # And the species...
        table = {
            "SMILES": [*self.species.keys()],
            "Count": [len(x) for x in self.species.values()],
            "mol %": [
                f"{100*len(x)/self.n_molecules:.2f}" for x in self.species.values()
            ],
        }
        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="simple_outline",
            disable_numparse=True,
            colalign=(
                "center",
                "decimal",
                "decimal",
            ),
        )
        length = len(tmp.splitlines()[0])
        text = "\n"
        text += "Species".center(length)
        text += "\n"
        text += tmp
        text += "\n"
        printer.important(__(text, indent=8 * " ", wrap=False, dedent=False))

        # Find the handler for job.out and set the level up
        job_handler = None
        out_handler = None
        for handler in job.handlers:
            if (
                isinstance(handler, logging.FileHandler)
                and "job.out" in handler.baseFilename
            ):
                job_handler = handler
                job_level = job_handler.level
                job_handler.setLevel(printing.JOB)
            elif isinstance(handler, logging.StreamHandler):
                out_handler = handler
                out_level = out_handler.level
                out_handler.setLevel(printing.JOB)

        # Make sure the subflowchart has the data from the parent flowchart
        self.subflowchart.root_directory = self.flowchart.root_directory
        self.subflowchart.executor = self.flowchart.executor
        self.subflowchart.in_jobserver = self.subflowchart.in_jobserver

        # Get the first real node
        first_node = self.subflowchart.get_node("1").next()

        # Ensure the nodes have their options
        node = first_node
        while node is not None:
            node.all_options = self.all_options
            node = node.next()

        # Loop over the runs
        nruns = P["nruns"]
        fmt = f"0{len(str(nruns))}"
        for run in range(1, nruns + 1):
            # Direct most output to iteration.out
            run_dir = Path(self.directory) / f"run_{run:{fmt}}"
            run_dir.mkdir(parents=True, exist_ok=True)

            # A handler for the file
            if self._file_handler is not None:
                self._file_handler.close()
                job.removeHandler(self._file_handler)
            path = run_dir / "Run.out"
            path.unlink(missing_ok=True)
            self._file_handler = logging.FileHandler(path)
            self._file_handler.setLevel(printing.NORMAL)
            formatter = logging.Formatter(fmt="{message:s}", style="{")
            self._file_handler.setFormatter(formatter)
            job.addHandler(self._file_handler)

            # Add the run to the ids so the directory structure is reasonable
            self.flowchart.reset_visited()
            self.set_subids((*self._id, f"run_{run:{fmt}}"))

            # Run through the steps in the loop body
            node = first_node
            try:
                while node is not None:
                    node = node.run()
            except DeprecationWarning as e:
                printer.normal("\nDeprecation warning: " + str(e))
                traceback.print_exc(file=sys.stderr)
                traceback.print_exc(file=sys.stdout)
            except Exception as e:
                printer.job(f"Caught exception in run {run}: {str(e)}")
                with open(run_dir / "stderr.out", "a") as fd:
                    traceback.print_exc(file=fd)
                if "continue" in P["errors"]:
                    continue
                elif "exit" in P["errors"]:
                    break
                else:
                    raise
            else:
                self.process_run(run, run_dir, P)
                if job_handler is not None:
                    job_handler.setLevel(job_level)
                if out_handler is not None:
                    out_handler.setLevel(out_level)

                self.analyze_run(P=P, run=run)

                if job_handler is not None:
                    job_handler.setLevel(printing.JOB)
                if out_handler is not None:
                    out_handler.setLevel(printing.JOB)

            self.logger.debug(f"End of run {run}")

        # Remove any redirection of printing.
        if self._file_handler is not None:
            self._file_handler.close()
            job.removeHandler(self._file_handler)
            self._file_handler = None
        if job_handler is not None:
            job_handler.setLevel(job_level)
        if out_handler is not None:
            out_handler.setLevel(out_level)

        # Analyze the results
        self.analyze(P=P, style="full")

        # Add other citations here or in the appropriate place in the code.
        # Add the bibtex to data/references.bib, and add a self.reference.cite
        # similar to the above to actually add the citation to the references.

        return next_node

    def parse_molecules(self):
        """Get the lists of molecules of each type.

        Returns
        -------
        {str, [int]}
            The list of molecule numbers (0-based) for each canonical SMILES string
            for a species.
        """
        # Create the subsets for the molecules
        SMILES = self.configuration.get_molecule_smiles()
        self._n_molecules = len(SMILES)
        result = {}
        for molecule, smiles in enumerate(SMILES):
            if smiles not in result:
                result[smiles] = [molecule]
            else:
                result[smiles].append(molecule)
        return result

    def process_run(self, run, run_dir, P=None):
        """Get the positions and velocities from the run and do initial processing.

        Parameters
        ----------
        run : int
            The run number
        run_dir : pathlib.Path
            The toplevel directory of the run.
        P : [str: any]
           Dictionary of control parameters, or default of None in which case the
           standard parameters are used.
        """
        if P is None:
            P = self.parameters.current_values_to_dict(
                context=seamm.flowchart_variables._data
            )

        # Get the state variables if saved in e.g. MD run
        if run == 1:
            # Reset any data of the system's state since this is a new calculation
            self._state = {}

        self._state_vars = []
        for variable in (
            "T",
            "T,stderr",
            "P",
            "P,stderr",
            "density",
            "density,stderr",
            "a",
            "a,stderr",
            "V",
            "V,stderr",
        ):
            if self.variable_exists(variable):
                if variable not in self._state:
                    self._state[variable] = []
                self._state[variable].append(self.get_variable(variable))
                self._state_vars.append(variable)
                # If they don't have standard errors, we still need space for them
                if ",stderr" not in variable:
                    var2 = variable + ",stderr"
                    if not self.variable_exists(var2):
                        if var2 not in self._state:
                            self._state[var2] = []
                        self._state[var2].append("")
                        self._state_vars.append(var2)

        n_species = len(self.species)
        if self._use_msd:
            tic = time.perf_counter_ns()
            paths = sorted(run_dir.glob("**/com_positions.trj"))

            if len(paths) == 0:
                raise RuntimeError(f"There is no com position data for run {run}.")
            elif len(paths) > 1:
                raise NotImplementedError(
                    f"Cannot handle multiple com position files from run {run}."
                )

            species = [x for x in self.species.values()]

            if self._msds is None:
                self._msds = [[] for i in range(n_species)]
                self._msd_errs = [[] for i in range(n_species)]

            metadata, result = read_vector_trajectory(paths[0])
            self._msd_dt = Q_(metadata["dt"], metadata["tunits"])
            msd, err = compute_msd(result, species)
            self._msd_samples = result.shape[0]

            for i in range(n_species):
                self._msds[i].append(msd[i])
                self._msd_errs[i].append(err[i])

            if True:
                # Analyze separately
                self._msd = msd
                self._msd_err = err
            else:
                # Average the MSD's then analyze
                if run == 1:
                    self._msd = msd
                    self._msd_err = err
                else:
                    for i in range(n_species):
                        tmp = np.stack(self._msds[i])
                        self._msd[i] = np.average(tmp, axis=0)
                        self._msd_err[i] = np.std(tmp, axis=0)
            toc = time.perf_counter_ns()
            _time = round((toc - tic) / 1e9, 2)
            self._msd_times.append(_time)
            if "t_msd" in self.results:
                self.results["t_msd"] += _time
            else:
                self.results["t_msd"] = _time
            self.logger.info(f"compute_msd: {_time:.3f}")

        if self._use_velocity:
            tic = time.perf_counter_ns()
            paths = sorted(run_dir.glob("**/com_velocities.trj"))

            if len(paths) == 0:
                raise RuntimeError(f"There is no com velocity data for run {run}.")
            elif len(paths) > 1:
                raise NotImplementedError(
                    f"Cannot handle multiple com velocity files from run {run}."
                )

            species = [x for x in self.species.values()]

            if self._Ms is None:
                self._Ms = [[] for i in range(n_species)]
                self._M_errs = [[] for i in range(n_species)]

            metadata, result = read_vector_trajectory(paths[0])
            self._velocity_dt = Q_(metadata["dt"], metadata["tunits"])

            # Limit the lengths of the data
            n = result.shape[0]
            if "maximum Helfand Integral length" in P:
                m = min(n // 2, P["maximum Helfand Integral length"])
            else:
                m = min(n // 2, 5000)
            self._helfand_integral_length = m

            # Convert units and remember the factor of 2 in the Helfand moments
            v_sq = Q_("Å^2/fs^2")
            constants = (v_sq * self._velocity_dt.to("fs") ** 2).m_as("Å^2") / 2

            M, err = create_helfand_moments(result, species, m=m)
            for i in range(n_species):
                M[i] *= constants
                err[i] *= constants
                self._Ms[i].append(M[i])
                self._M_errs[i].append(err[i])

            if True:
                # analyze separately
                self._M = M
                self._M_err = err
            else:
                # Average the Helfand moments, then analyze
                if run == 1:
                    self._M = M
                    self._M_err = err
                else:
                    for i in range(n_species):
                        tmp = np.stack(self._Ms[i])
                        self._M[i] = np.average(tmp, axis=0)
                        self._M_err[i] = np.std(tmp, axis=0)
            toc = time.perf_counter_ns()
            _time = round((toc - tic) / 1e9, 3)
            if "t_hm" in self.results:
                self.results["t_hm"] += _time
            else:
                self.results["t_hm"] = _time
            self._helfand_integral_times.append(_time)
            self.logger.info(f"create_helfand_moments: {_time:.3f}")

    def set_id(self, node_id=()):
        """Sequentially number the subnodes"""
        self.logger.debug("Setting ids for subflowchart {}".format(self))
        if self.visited:
            return None
        else:
            self.visited = True
            self._id = node_id
            self.set_subids(self._id)
            return self.next()

    def set_subids(self, node_id=()):
        """Set the ids of the nodes in the subflowchart"""
        self.subflowchart.reset_visited()
        node = self.subflowchart.get_node("1").next()
        n = 1
        while node is not None:
            node = node.set_id((*node_id, str(n)))
            n += 1
