# -*- coding: utf-8 -*-

"""Non-graphical part of the Diffusivity step in a SEAMM flowchart

The self-diffusion coefficients of the species in a fluid, from the mean square
displacement (MSD) of the molecules' centres of mass, from their single-molecule
Helfand moments (the square of the integral of the velocity), or from both, with the
Yeh-Hummer correction for the size of the periodic cell.

The MD is run before this step -- in a Loop for independent runs -- or read from
files; the step analyzes it.
"""

import importlib.resources
import logging
from math import floor, log10, sqrt
from pathlib import Path
import time

import numpy as np
from tabulate import tabulate

from seamm_md_analysis import (
    Estimate,
    component_check,
    fit_slope,
    is_cubic,
    isotropy,
    msd_by_species,
    open_md,
    run_state,
    species_helfand_moments,
    yeh_hummer_correction,
)
from seamm_md_analysis.plots import add_helfand_trace, add_msd_trace
from seamm_md_analysis.requirements import Need, Requirements
from seamm_md_analysis.step_support import (
    check_md_flowchart,
    md_inputs,
    preceding_viscosity,
    print_details,
    read_files,
    report_state,
    size_warnings,
    state_data,
)
from seamm_md_analysis.transport import YEH_HUMMER_XI

import diffusivity_step
import molsystem
import seamm
import seamm_util
from seamm_util import Q_
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
path = importlib.resources.files("diffusivity_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)

#: Å²/ps in m²/s
A2_PER_PS = Q_(1.0, "Å^2/ps").m_as("m^2/s")
#: The name of each method in the output
METHODS = {"MSD": "MSD", "HM": "Helfand moments", "both": "Combined"}
#: The suffix of each method's results; the best estimate has none
SUFFIX = {"MSD": " (MSD)", "HM": " (HM)", "both": ""}
#: The directions
DIRECTIONS = ("x", "y", "z")
#: The fewest frames of a block of a single run
MIN_FRAMES = 10
#: The number of frames read at a time when computing centres of mass
CHUNK = 100

_SUPERSCRIPTS = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def power_of_ten(exponent):
    """10 to an integer power, as text, e.g. "10⁻⁹"."""
    return "10" + str(exponent).translate(_SUPERSCRIPTS)


class Diffusivity(seamm.Node):
    """
    The non-graphical part of a Diffusivity step in a flowchart.

    Attributes
    ----------
    parameters : DiffusivityParameters
        The control parameters for Diffusivity.
    paths : [pathlib.Path]
        The md.h5md file of each run analyzed.

    See Also
    --------
    TkDiffusivity, DiffusivityParameters
    """

    def __init__(
        self,
        flowchart=None,
        title="Diffusivity",
        extension=None,
        logger=logger,
    ):
        """A step for Diffusivity in a SEAMM flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.
        title: str
            The name displayed in the flowchart.
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes
        """
        logger.debug(f"Creating Diffusivity {self}")
        super().__init__(
            flowchart=flowchart,
            title="Diffusivity",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = diffusivity_step.metadata
        self.parameters = diffusivity_step.DiffusivityParameters()
        self._reset()

    @property
    def version(self):
        """The semantic version of this module."""
        return diffusivity_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return diffusivity_step.__git_revision__

    def _reset(self):
        """Clear the data of earlier runs, e.g. from an earlier loop iteration."""
        self.paths = []  # The md.h5md file of each run
        self._viscosity = None  # The viscosity for the Yeh-Hummer correction

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = "diffusivity-step"
        parser = seamm_util.getParser()

        parser_exists = parser.exists(parser_name)
        super().create_parser(name=parser_name)
        if not parser_exists:
            parser.add_argument(
                parser_name,
                "--html",
                action="store_true",
                help="whether to write out html files for graphs, etc.",
            )
        return self.next()

    def description_text(self, P=None, short=False):
        """Create the text description of what this step will do.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        short : bool
            Not used: the description is always short.

        Returns
        -------
        str
            A description of the current step.
        """
        if P is None:
            P = self.parameters.values_to_dict()

        approach = str(P["approach"])
        use_msd, use_hm = self._approaches(approach)
        text = "Calculate the self-diffusion coefficients of the species from "
        if use_msd and use_hm:
            text += (
                "both the mean square displacement (MSD) and the Helfand moments of "
                "the molecules"
            )
        elif use_hm:
            text += "the Helfand moments of the molecules"
        else:
            text += "the mean square displacement (MSD) of the molecules"
        source = P["md source"]
        if source == read_files:
            text += f", for the MD in {P['md files']}."
        else:
            text += f", for {source} earlier in the flowchart."

        if use_msd:
            text += (
                f" The slope of the MSD is fit from {P['msd_fit_start']} to "
                f"{P['msd_fit_end']} of the curve."
            )
        if use_hm:
            text += (
                " The Helfand moments are integrated over at most "
                f"{P['maximum Helfand Integral length']} steps -- the cost grows as "
                "the square of the number -- and their slope is fit from "
                f"{P['helfand_fit_start']} to {P['helfand_fit_end']} of the curve."
            )

        correction = P["hydrodynamic correction"]
        if isinstance(correction, str) and self.is_expr(correction):
            text += (
                f" The expression '{correction}' will determine whether to use the "
                "Yeh-Hummer hydrodynamic correction for the finite cell size."
            )
        else:
            if isinstance(correction, str):
                correction = "yes" in correction
            if correction:
                text += (
                    " The diffusion coefficients are also given corrected for the "
                    "finite size of the cell with the Yeh-Hummer correction"
                )
                if P.get("viscosity source", "given value") == "given value":
                    text += f", using a viscosity of {P['viscosity']}."
                else:
                    text += (
                        ", using the viscosity from the preceding Viscosity step, if "
                        "there is one."
                    )

        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    @staticmethod
    def _approaches(approach):
        """Whether to use the MSD and the Helfand moments."""
        approach = str(approach).lower()
        both = approach == "both"
        return both or "msd" in approach, both or "helfand" in approach

    def run(self):
        """Run the diffusivity step.

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        next_node = super().run(printer)
        self._reset()

        # Get the values of the parameters, dereferencing any variables
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        printer.important(__(self.description_text(P), indent=self.indent))

        _, configuration = self.get_system_configuration()
        if configuration.periodicity != 3:
            raise RuntimeError(
                "The Diffusivity step needs a periodic system: the diffusion "
                "coefficient of a cluster or molecule is not a bulk property, and "
                "the finite-size correction needs the cell."
            )

        # The viscosity for the Yeh-Hummer correction, and where it came from
        self._viscosity, text = self._viscosity_for_correction(P)
        if text is not None:
            printer.important(__(text, indent=self.indent + 4 * " "))
            printer.important("")

        self.paths = self._md_paths(P)
        n = len(self.paths)
        printer.job(
            __(
                f"Analyzing {n} MD run{'s' if n > 1 else ''}: "
                + ", ".join(str(p) for p in self.paths),
                indent=self.indent + 4 * " ",
            )
        )

        self.analyze(P=P, configuration=configuration)
        return next_node

    def _md_paths(self, P):
        """The md.h5md files to analyze: those of the runs that sampled the first
        need, as checked before the job ran."""
        for name in self.md_requirements(P).needs[0].names:
            try:
                return md_inputs(self, P, name=name)
            except RuntimeError:
                continue
        return md_inputs(self, P)

    def _viscosity_for_correction(self, P):
        """The viscosity for the Yeh-Hummer correction, and where it came from.

        Parameters
        ----------
        P : dict
            The values of the parameters.

        Returns
        -------
        dict or None
            The viscosity: "value" and "stderr" (Pa s, the error may be None), "T"
            (K, or None) and "origin", for the output; None if the correction is
            not applied.
        str or None
            A note for the output.
        """
        if not P["hydrodynamic correction"]:
            return None, None
        if P.get("viscosity source", "given value") == "given value":
            eta = P["viscosity"]
            if eta.magnitude <= 0:
                return (
                    None,
                    "Note: the Yeh-Hummer correction is skipped because the "
                    "viscosity given is not positive.",
                )
            record = {
                "value": eta.m_as("Pa*s"),
                "stderr": None,
                "T": None,
                "origin": "given in this step",
            }
            text = f"The Yeh-Hummer correction uses the given viscosity, {eta:~P}."
            return record, text
        published = preceding_viscosity(self)
        if published is None:
            return (
                None,
                "Note: the Yeh-Hummer correction is skipped because no Viscosity "
                "step ran before this one. Add one, or give the viscosity.",
            )
        origin = (
            f"from step {published['id']} ({published['step']}, "
            f"{published['method']}"
        )
        if published.get("T") is not None:
            origin += f", at {published['T']:.1f} K"
        origin += ")"
        record = {
            "value": published["value"],
            "stderr": published.get("stderr"),
            "T": published.get("T"),
            "origin": origin,
        }
        eta = Q_(published["value"], "Pa*s").m_as("mPa*s")
        text = f"The Yeh-Hummer correction uses the viscosity, {eta:.4g}"
        if record["stderr"] is not None:
            error = Q_(record["stderr"], "Pa*s").m_as("mPa*s")
            text += f" ± {error:.2g} (standard error)"
        text += f" mPa*s, {origin}."
        return record, text

    # ------------------------------------------------------------------
    # Reading the runs
    # ------------------------------------------------------------------
    @staticmethod
    def species(configuration):
        """The molecules of each species.

        Returns
        -------
        {str: [int]}
            The molecules (0-based) of each species, by canonical SMILES, in the
            order of find_molecules.
        """
        result = {}
        for molecule, smiles in enumerate(configuration.get_molecule_smiles()):
            result.setdefault(smiles, []).append(molecule)
        return result

    def _read(self, path, molecules, masses, use_msd, use_hm):
        """What the analysis needs from one run.

        Returns
        -------
        dict
            "position" (Å) and "velocity" (Å/fs) of the molecules' centres of
            mass, (n_frames, n_molecules, 3), with "position dt" and "velocity dt"
            (fs); the mean temperature "T" (K), volume "V" (Å³) and "L", the cube
            root of the volume (Å); and "isotropic", (bool, reason), and "cubic".
        """
        result = {"path": path}
        with open_md(path) as md:
            if not md.periodic:
                raise RuntimeError(
                    f"The MD in {path} is not periodic, so it has no diffusivity."
                )
            state = run_state(md)
            result["T"] = state.temperature
            result["V"] = state.volume
            result["L"] = None if not state.volume else state.volume ** (1 / 3)
            result["isotropic"] = isotropy(md)
            result["cubic"] = is_cubic(md.cell())
            for name, units, use in (
                ("position", "Å", use_msd),
                ("velocity", "Å/fs", use_hm),
            ):
                if use:
                    values, dt = self._molecular(
                        md, name, units, molecules, masses, path
                    )
                    result[name] = values
                    result[name + " dt"] = dt
        return result

    @staticmethod
    def _molecular(md, name, units, molecules, masses=None, path=""):
        """The positions or velocities of the molecules' centres of mass.

        From the molecules' group if the MD sampled it, otherwise from the atoms':
        the positions are unwrapped, so the molecules are whole.

        Parameters
        ----------
        md : seamm_md_analysis.MDData
            The run.
        name : str
            "position" or "velocity".
        units : str
            The units to return them in.
        molecules : [[int]]
            The atoms (0-based) of each molecule.
        masses : array_like = None
            The masses of the atoms, if the file has none.
        path : str
            The file, for messages.

        Returns
        -------
        numpy.ndarray(n_frames, n_molecules, 3), float
            The values, and the interval between frames in fs (None for one frame).
        """
        n_molecules = len(molecules)
        if md.has_particles(name, "molecules"):
            series = md.particle_series(
                name, group="molecules", units=units, time_units="fs"
            )
            if series.values.shape[1] != n_molecules:
                raise RuntimeError(
                    f"The MD in {path} has the {name}s of {series.values.shape[1]} "
                    f"molecules, but the structure has {n_molecules}."
                )
            dt = series.dt
            return series.values, None if dt is None else float(dt.m_as("fs"))
        if not md.has_particles(name, "all"):
            raise RuntimeError(
                f"The MD in {path} did not sample the {name}s of the molecules or "
                "of the atoms."
            )
        if name == "position":
            attrs = md.first["particles/all/position/value"].attrs
            if not attrs.get("seamm_unwrapped", 0):
                raise RuntimeError(
                    f"The atomic positions in {path} are not unwrapped, so the "
                    "centres of mass of the molecules cannot be found. Sample the "
                    "positions of the molecules' centres of mass."
                )
        file_masses = md.masses("all")
        if file_masses is not None:
            masses = file_masses
        if masses is None:
            raise RuntimeError(f"The MD in {path} has no masses for the atoms.")
        masses = np.asarray(masses, dtype=float)
        n_atoms = len(md.ids("all"))
        sizes = [len(m) for m in molecules]
        if n_atoms != sum(sizes) or len(masses) != n_atoms:
            raise RuntimeError(
                f"The MD in {path} has {n_atoms} atoms, but the structure has "
                f"{sum(sizes)}."
            )
        atoms = np.concatenate([np.asarray(m, dtype=int) for m in molecules])
        starts = np.cumsum([0] + sizes[:-1])
        weights = masses[atoms]
        weights = weights / np.repeat(np.add.reduceat(weights, starts), sizes)
        values = []
        times = []
        for _, t, chunk in md.iter_frames(
            group="all",
            fields=(name,),
            chunk=CHUNK,
            units={name: units},
            time_units="fs",
        ):
            weighted = chunk[name][:, atoms, :] * weights[None, :, None]
            values.append(np.add.reduceat(weighted, starts, axis=1))
            times.append(np.asarray(t, dtype=float))
        values = np.concatenate(values)
        times = np.concatenate(times)
        dt = float(np.median(np.diff(times))) if len(times) > 1 else None
        return values, dt

    @staticmethod
    def _split_into_blocks(run, n_blocks):
        """Divide one run into consecutive blocks, analyzed as separate samples."""
        n_blocks = max(1, int(n_blocks))
        blocks = [dict(run) for _ in range(n_blocks)]
        for name in ("position", "velocity"):
            values = run.get(name)
            if values is None:
                continue
            n = len(values) // n_blocks
            if n < MIN_FRAMES:
                raise RuntimeError(
                    f"The run has only {len(values)} frames of the {name}s, too few "
                    f"for {n_blocks} blocks."
                )
            for i in range(n_blocks):
                blocks[i][name] = values[i * n : (i + 1) * n]
        return blocks

    # ------------------------------------------------------------------
    # The analysis
    # ------------------------------------------------------------------
    @staticmethod
    def _fit_msd(sample, species, start, end):
        """The diffusion coefficient in each direction from the MSD.

        Returns
        -------
        numpy.ndarray(n_species, 3)
            D in m^2/s.
        [numpy.ndarray(n, 4)]
            The MSD of each species, in Å^2.
        numpy.ndarray(n)
            The times, in ps.
        """
        positions = sample["position"]
        msd, _ = msd_by_species(positions, species)
        ts = np.arange(positions.shape[0]) * sample["position dt"] / 1000
        # The Einstein relation: MSD = 2 D t in each direction
        D = np.array(
            [
                [
                    fit_slope(curve[:, i], ts, start=start, end=end).slope
                    * A2_PER_PS
                    / 2
                    for i in range(3)
                ]
                for curve in msd
            ]
        )
        return D, msd, ts

    @staticmethod
    def _fit_helfand(sample, species, max_length, start, end):
        """The diffusion coefficient in each direction from the Helfand moments.

        Returns
        -------
        numpy.ndarray(n_species, 3)
            D in m^2/s.
        [numpy.ndarray(m, 4)]
            The Helfand moments of each species, in Å^2.
        numpy.ndarray(m)
            The times, in ps.
        int
            The length of the integration, in frames.
        """
        velocities = sample["velocity"]
        dt = sample["velocity dt"]
        m = min(velocities.shape[0] // 2, int(max_length))
        M, _ = species_helfand_moments(velocities, species, m=m)
        # Å²/fs² × fs², and the factor of 1/2 in the moments
        M = [moment * dt**2 / 2 for moment in M]
        ts = np.arange(m) * dt / 1000
        D = np.array(
            [
                [
                    fit_slope(curve[:, i], ts, start=start, end=end).slope * A2_PER_PS
                    for i in range(3)
                ]
                for curve in M
            ]
        )
        return D, M, ts, m

    def analyze(self, indent="", P=None, configuration=None, **kwargs):
        """Analyze, print, plot and store the diffusion coefficients.

        Parameters
        ----------
        indent : str
            An extra indentation for the output.
        P : dict
            The values of the parameters.
        configuration : molsystem._Configuration
            The structure that was simulated.
        """
        if configuration is None:
            raise RuntimeError("There is no structure to analyze.")
        species = self.species(configuration)
        names = list(species)
        species_lists = list(species.values())
        molecules = configuration.find_molecules(as_indices=True)
        masses = configuration.atoms.atomic_masses
        self._print_species(species, len(molecules))

        use_msd, use_hm = self._approaches(P["approach"])
        methods = (["MSD"] if use_msd else []) + (["HM"] if use_hm else [])
        nruns = len(self.paths)
        source = "runs" if nruns > 1 else "blocks"

        # One run at a time, keeping only the results, so that the memory needed
        # does not grow with the number of runs. Each run, or each block of a
        # single run, is a sample.
        runs = []
        samples = []
        D = {m: [] for m in methods}  # (n_samples, n_species, 3)
        curves = {m: [] for m in methods}
        timings = {m: 0.0 for m in methods}
        integral_length = None
        for path in self.paths:
            run = self._read(path, molecules, masses, use_msd, use_hm)
            if nruns > 1:
                pieces = [run]
            else:
                pieces = self._split_into_blocks(run, P["number of blocks"])
            for piece in pieces:
                if use_msd:
                    tic = time.perf_counter()
                    d, msd, ts = self._fit_msd(
                        piece, species_lists, P["msd_fit_start"], P["msd_fit_end"]
                    )
                    timings["MSD"] += time.perf_counter() - tic
                    D["MSD"].append(d)
                    curves["MSD"].append((msd, ts))
                if use_hm:
                    tic = time.perf_counter()
                    d, M, ts, integral_length = self._fit_helfand(
                        piece,
                        species_lists,
                        P["maximum Helfand Integral length"],
                        P["helfand_fit_start"],
                        P["helfand_fit_end"],
                    )
                    timings["HM"] += time.perf_counter() - tic
                    D["HM"].append(d)
                    curves["HM"].append((M, ts))
                samples.append(self._summary(piece))
            runs.append(self._summary(run))
            del run, pieces
        n_samples = len(samples)
        D = {m: np.array(values) for m, values in D.items()}
        if use_msd and use_hm:
            # The average of the two methods in each sample
            methods.append("both")
            D["both"] = (D["MSD"] + D["HM"]) / 2
        best = methods[-1]

        def estimate(values):
            values = np.asarray(values, dtype=float)
            if len(values) < 2:
                return Estimate(float(values[0]), None, 1, source=source)
            return Estimate.from_samples(values, source=source)

        # The Yeh-Hummer correction for each sample, and its error from that of η
        correction, notes = self._correction(samples)

        def corrected(values):
            if correction is None:
                return None
            result = estimate(np.asarray(values) + correction["samples"])
            if result.stderr is not None and correction["stderr"] is not None:
                result.stderr = sqrt(result.stderr**2 + correction["stderr"] ** 2)
            return result

        # The estimates: D is the average of the directions in each sample
        estimates = {}
        for m in methods:
            for s in range(len(names)):
                values = D[m][:, s, :]
                average = values.mean(axis=1)
                estimates[m, s, "D"] = (estimate(average), corrected(average))
                for i, c in enumerate(DIRECTIONS):
                    estimates[m, s, f"D{c}"] = (
                        estimate(values[:, i]),
                        corrected(values[:, i]),
                    )

        largest = max(
            (
                abs(estimates[m, s, "D"][0].value)
                for m in methods
                for s in range(len(names))
            ),
            default=0.0,
        )
        exponent = floor(log10(largest)) if largest > 0 else 0
        scale = 10.0**exponent
        unit = f"{power_of_ten(exponent)} m²/s"

        isotropic = True
        why = ""
        for run in runs:
            if not run["isotropic"][0]:
                isotropic, why = run["isotropic"]
                break

        self._print_table(
            names, methods, estimates, isotropic, correction is not None, scale, unit
        )

        # What the table shows
        lines = []
        what = (
            f"{nruns} independent runs"
            if source == "runs"
            else f"{n_samples} blocks of one run"
        )
        lengths = self._sample_lengths(samples)
        text = f"From {what}"
        if lengths:
            text += f" of {lengths} each"
        text += ". " + estimates[best, 0, "D"][0].describe_error()
        lines.append(text)
        if "both" in methods:
            lines.append(
                f"Combined is the average of the MSD and Helfand-moment values in "
                f"each {source[:-1]}."
            )
        if isotropic:
            lines.append(self._isotropy_line(names, methods, D))
        else:
            lines.append(
                f"The run is not treated as isotropic ({why}), so the components "
                "are given; D is their average."
            )
        if correction is not None:
            lines.append(self._correction_text(correction, scale, unit))
        lines.extend(notes)
        for line in lines:
            if line:
                printer.normal(__(line, indent=8 * " "))
        printer.normal("")

        self._print_samples(names, methods, D, samples, source, scale, unit)
        text = []
        if "MSD" in timings:
            text.append(
                f"The MSD, of {samples[0]['position frames']} frames per "
                f"{source[:-1]}, took {timings['MSD']:.2f} s in all"
            )
        if "HM" in timings:
            text.append(
                f"the Helfand moments, integrated over {integral_length} frames, "
                f"took {timings['HM']:.2f} s"
            )
        if text:
            text = "; ".join(text) + "."
            text = text[0].upper() + text[1:]
            print_details(printer, str(__(text, indent=8 * " ")) + "\n")

        state = report_state(self, printer, self.paths)

        self._plot(names, curves, P, scale)

        # The results. "D" is corrected for the cell size if the correction is
        # applied, as before; "D (box)" is the value in the periodic cell.
        data = {
            "approach": P["approach"],
            "nruns": nruns,
        }
        if use_msd:
            data["msd_fit_start"] = P["msd_fit_start"]
            data["msd_fit_end"] = P["msd_fit_end"]
            data["msd samples"] = samples[0]["position frames"]
            data["t_msd"] = round(timings["MSD"], 3)
        if use_hm:
            data["helfand_fit_start"] = P["helfand_fit_start"]
            data["helfand_fit_end"] = P["helfand_fit_end"]
            data["Helfand integral length"] = integral_length
            data["t_hm"] = round(timings["HM"], 3)
        keys = [(m, SUFFIX[m]) for m in methods if m != "both"] + [(best, "")]
        for m, suffix in keys:
            for s, name in enumerate(names):
                for component in ("D", "Dx", "Dy", "Dz"):
                    box, fixed = estimates[m, s, component]
                    value = box if fixed is None else fixed
                    key = f"{component} {{key}}{suffix}"
                    data.setdefault(key, {})[name] = value.value
                    if value.stderr is not None:
                        data.setdefault(key + ",stderr", {})[name] = value.stderr
        for s, name in enumerate(names):
            box = estimates[best, s, "D"][0]
            data.setdefault("D (box) {key}", {})[name] = box.value
            if box.stderr is not None:
                data.setdefault("D (box) {key},stderr", {})[name] = box.stderr
        if correction is not None:
            data["Yeh-Hummer correction"] = correction["value"]
            if correction["stderr"] is not None:
                data["Yeh-Hummer correction,stderr"] = correction["stderr"]
            data["viscosity for correction"] = Q_(correction["eta"], "Pa*s").m_as(
                "mPa*s"
            )
        inverse = [1 / r["L"] for r in runs if r["L"]]
        if inverse:
            one_over_L = Estimate.from_samples(inverse)
            data["1/L"] = one_over_L.value
            if one_over_L.stderr:
                data["1/L,stderr"] = one_over_L.stderr
        if state is not None:
            data.update(state_data(state))

        self.model = self._find_model()
        self.store_results(configuration=configuration, data=data)
        return data

    @staticmethod
    def _summary(sample):
        """A run or block without its trajectory: the number of frames of the
        positions and velocities, the interval between them, and the state."""
        result = {k: v for k, v in sample.items() if k not in ("position", "velocity")}
        for name in ("position", "velocity"):
            if sample.get(name) is not None:
                result[name + " frames"] = len(sample[name])
        return result

    @staticmethod
    def _sample_lengths(samples):
        """The length of the samples, e.g. "100 ps", or "" if not known."""
        for name in ("position", "velocity"):
            n = samples[0].get(name + " frames")
            dt = samples[0].get(name + " dt")
            if n is not None and dt is not None:
                return f"{n * dt / 1000:.4g} ps"
        return ""

    def _correction(self, samples):
        """The Yeh-Hummer correction of each sample, from its temperature and the
        size of the cell.

        Returns
        -------
        dict or None
            "samples" (m^2/s, one per sample), its mean "value", "stderr" from the
            error of the viscosity (or None), "eta" (Pa s), its "eta stderr" and
            "origin", and the mean "T" (K) and "L" (Å); None if not applied.
        [str]
            Notes for the output.
        """
        notes = []
        record = self._viscosity
        if record is None:
            return None, notes
        T = [s["T"] for s in samples]
        L = [s["L"] for s in samples]
        if any(t is None for t in T) or any(x is None for x in L):
            notes.append(
                "Note: the runs did not record the temperature or the cell, so the "
                "Yeh-Hummer correction is skipped."
            )
            return None, notes
        eta = record["value"]
        corrections = np.array(
            [yeh_hummer_correction(t, x * 1.0e-10, eta) for t, x in zip(T, L)]
        )
        value = float(corrections.mean())
        stderr = None
        if record.get("stderr") is not None:
            stderr = value * record["stderr"] / eta
        result = {
            "samples": corrections,
            "value": value,
            "stderr": stderr,
            "eta": eta,
            "eta stderr": record.get("stderr"),
            "origin": record["origin"],
            "T": float(np.mean(T)),
            "L": float(np.mean(L)),
            "varies": len(samples) > 1 and np.ptp(L) > 1.0e-6 * np.mean(L),
        }
        if not all(s["cubic"] for s in samples):
            notes.append(
                "Note: the cell is not cubic, so the Yeh-Hummer correction, which is "
                "for a cubic cell, is approximate; L is the cube root of the volume."
            )
        if record.get("T") is not None and abs(result["T"] - record["T"]) > 2.0:
            notes.append(
                f"Note: the viscosity is for {record['T']:.1f} K but the MD here is "
                f"at {result['T']:.1f} K, so the Yeh-Hummer correction is "
                "approximate."
            )
        self.references.cite(
            raw=self._bibliography["Yeh2004"],
            alias="Yeh2004",
            module="diffusivity_step",
            level=1,
            note="The Yeh-Hummer correction of the diffusivity for the cell size.",
        )
        return result, notes

    @staticmethod
    def _correction_text(correction, scale, unit):
        """A description of the correction, its error and the viscosity used."""
        value = correction["value"] / scale
        text = (
            "D corrected = D (box) + ΔD, the Yeh-Hummer correction for the "
            "hydrodynamic interaction of a molecule with its periodic images: "
            "ΔD = ξ kB T / (6π η L) = "
        )
        if correction["stderr"] is not None:
            error = correction["stderr"] / scale
            text += f"({value:.3g} ± {error:.2g} standard error) × {unit}"
        else:
            text += f"{value:.3g} × {unit}"
        T = correction["T"]
        L = correction["L"]
        mean = "the mean " if correction["varies"] else ""
        text += (
            f", with ξ = {YEH_HUMMER_XI}, {mean}T = {T:.1f} K, {mean}L = "
            f"{L:.3f} Å (1/L = {1 / L:.5f} Å⁻¹) and η = "
            f"{Q_(correction['eta'], 'Pa*s').m_as('mPa*s'):.4g}"
        )
        if correction["eta stderr"] is not None:
            error = Q_(correction["eta stderr"], "Pa*s").m_as("mPa*s")
            text += f" ± {error:.2g} (standard error)"
        text += f" mPa·s, {correction['origin']}."
        if correction["stderr"] is not None:
            text += " The ± of D corrected includes the error of η."
        return text

    @staticmethod
    def _isotropy_line(names, methods, D):
        """One line saying whether the diffusion is the same in every direction."""
        checks = []
        for m in methods:
            if m == "both":
                continue
            for s, name in enumerate(names):
                labels = [f"D{c} of {name} ({METHODS[m]})" for c in DIRECTIONS]
                checks.append(component_check(D[m][:, s, :], names=labels))
        usable = [c for c in checks if c.p is not None]
        if len(usable) == 0:
            return checks[0].text if checks else ""
        if len(usable) == 1:
            return usable[0].text
        differ = [c for c in usable if not c.agree]
        if len(differ) == 0:
            what = " and method" if len([m for m in methods if m != "both"]) > 1 else ""
            worst = max(usable, key=lambda c: c.chi2 / c.dof)
            return (
                f"Isotropy: the components agree for every species{what} "
                f"(largest χ²/dof = {worst.chi2 / worst.dof:.2f}, p = {worst.p:.2g})."
            )
        return " ".join(c.text for c in differ)

    def _print_species(self, species, n_molecules):
        """The species and their numbers."""
        table = {
            "Species": list(species),
            "Count": [len(x) for x in species.values()],
            "mol %": [f"{100 * len(x) / n_molecules:.2f}" for x in species.values()],
        }
        text = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=("left", "decimal", "decimal"),
        )
        length = len(text.splitlines()[0])
        printer.important(
            __(
                "\n" + "Species".center(length) + "\n" + text + "\n",
                indent=8 * " ",
                wrap=False,
                dedent=False,
            )
        )

    def _print_table(
        self, names, methods, estimates, isotropic, with_correction, scale, unit
    ):
        """The main table: D of each species, in the box and corrected."""
        headers = ["Species", "Method"]
        if not isotropic:
            headers.append("Component")
        headers += ["D (box)" if with_correction else "D", "± (95%)"]
        if with_correction:
            headers += ["D corrected", "± (95%)"]

        def scaled(estimate):
            return Estimate(
                estimate.value / scale,
                None if estimate.stderr is None else estimate.stderr / scale,
                estimate.n,
                source=estimate.source,
            ).formatted()

        components = ["D"] if isotropic else ["D", "Dx", "Dy", "Dz"]
        rows = []
        for s, name in enumerate(names):
            for k, m in enumerate(methods):
                for j, component in enumerate(components):
                    row = [
                        name if k == 0 and j == 0 else "",
                        METHODS[m] if j == 0 else "",
                    ]
                    if not isotropic:
                        row.append(
                            component if component == "D" else f"D_{component[1]}"
                        )
                    box, fixed = estimates[m, s, component]
                    row.extend(scaled(box))
                    if with_correction:
                        row.extend(scaled(fixed))
                    rows.append(row)
        colalign = ["left", "left"] + ([] if isotropic else ["left"])
        colalign += ["decimal"] * (len(headers) - len(colalign))
        text = tabulate(
            rows,
            headers=headers,
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=colalign,
        )
        length = len(text.splitlines()[0])
        title = f"Self-diffusion coefficients ({unit})"
        printer.normal(
            __(
                "\n" + title.center(length) + "\n" + text + "\n",
                indent=8 * " ",
                wrap=False,
                dedent=False,
            )
        )

    def _print_samples(self, names, methods, D, samples, source, scale, unit):
        """The values of each run or block, in the detailed output."""
        what = "Run" if source == "runs" else "Block"
        headers = [
            what,
            "Species",
            "Method",
            "D_x",
            "D_y",
            "D_z",
            "D",
            "T (K)",
            "L (Å)",
        ]
        rows = []
        for i, sample in enumerate(samples):
            first = True
            for s, name in enumerate(names):
                for m in methods:
                    values = D[m][i, s, :] / scale
                    row = [
                        str(i + 1) if first else "",
                        name,
                        METHODS[m],
                        *(f"{v:.3f}" for v in values),
                        f"{values.mean():.3f}",
                    ]
                    if first and sample["T"] is not None:
                        row.append(f"{sample['T']:.2f}")
                    else:
                        row.append("")
                    if first and sample["L"] is not None:
                        row.append(f"{sample['L']:.3f}")
                    else:
                        row.append("")
                    rows.append(row)
                    first = False
        text = tabulate(
            rows,
            headers=headers,
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=("right", "left", "left", *(["decimal"] * 6)),
        )
        length = len(text.splitlines()[0])
        title = f"By {what.lower()} ({unit})".center(length)
        print_details(
            printer,
            str(
                __(
                    "\n" + title + "\n" + text + "\n",
                    indent=8 * " ",
                    wrap=False,
                    dedent=False,
                )
            ),
        )

    def _plot(self, names, curves, P, scale):
        """The MSD and Helfand moments, averaged over the runs or blocks, with the
        fitted lines."""
        settings = {
            "MSD": (
                "Mean Square Displacement",
                "MSD",
                "MSD (Å^2)",
                ("msd_fit_start", "msd_fit_end"),
                (2, 6),
                add_msd_trace,
                "MSD.graph",
            ),
            "HM": (
                "Helfand Moments",
                "HelfandMoments",
                "M (Å^2)",
                ("helfand_fit_start", "helfand_fit_end"),
                (1, 3),
                add_helfand_trace,
                "HelfandMoments.graph",
            ),
        }
        for method, per_sample in curves.items():
            title, name, ylabel, fraction, divisors, add, filename = settings[method]
            start, end = P[fraction[0]], P[fraction[1]]
            figure = self.create_figure(
                module_path=("seamm",),
                template="line.graph_template",
                title=title,
            )
            plot = figure.add_plot(name)
            x_axis = plot.add_axis("x", label="Time (ps)")
            y_axis = plot.add_axis("y", label=ylabel, anchor=x_axis)
            x_axis.anchor = y_axis

            ts = per_sample[0][1]
            same = all(c[1].shape == ts.shape for c in per_sample)
            for s, species in enumerate(names):
                if same:
                    curve = np.mean([c[0][s] for c in per_sample], axis=0)
                else:
                    curve = per_sample[0][0][s]
                fit = []
                for i in range(curve.shape[1]):
                    line = fit_slope(curve[:, i], ts, start=start, end=end)
                    divisor = divisors[0] if i < 3 else divisors[1]
                    d = line.slope * A2_PER_PS / divisor
                    fit.append(
                        {
                            "xs": list(line.xs),
                            "ys": list(line.ys),
                            "D_s": round(d / scale, 3),
                            "scale": scale,
                        }
                    )
                add(plot, x_axis, y_axis, species, curve, ts.tolist(), fit=fit)

            figure.grid_plots(name)
            path = Path(self.directory) / filename
            figure.dump(path)
            if "html" in self.options and self.options["html"]:
                figure.template = "line.html_template"
                figure.dump(path.with_suffix(".html"))

    def _find_model(self):
        """The model (forcefield or model chemistry) of the MD, for properties."""
        try:
            if self.variable_exists("_model_chemistry"):
                return self.get_variable("_model_chemistry")["level"]
            ff = self.get_variable("_forcefield")
            if ff == "OpenKIM":
                return "OpenKIM/" + self.get_variable("_OpenKIM_Potential")
            if ff == "PyTorch":
                return "PyTorch/" + self.get_variable("_pytorch_model")
            return ff.current_forcefield
        except Exception:
            return "unknown"

    # ------------------------------------------------------------------
    # What the MD must provide
    # ------------------------------------------------------------------
    def md_requirements(self, P=None):
        """What the analysis needs from the MD.

        Parameters
        ----------
        P : dict = None
            The parameters, as values; by default the current ones.

        Returns
        -------
        seamm_md_analysis.Requirements
        """
        approach = self._static("approach") if P is None else P["approach"]
        use_msd, use_hm = self._approaches(approach)
        needs = []
        if use_msd:
            needs.append(
                Need(
                    ["molecules/position", "all/position"],
                    max_time="1 ps",
                    sample_time="100 fs",
                    why="the mean square displacement",
                )
            )
        if use_hm:
            needs.append(
                Need(
                    ["molecules/velocity", "all/velocity"],
                    max_time="20 fs",
                    sample_time="10 fs",
                    why="the Helfand moments",
                )
            )
        return Requirements(needs=needs, periodic=True, dynamics=True)

    def _static(self, key):
        """A parameter's value without evaluating variables, or None."""
        parameter = self.parameters[key]
        if parameter.is_expr:
            return None
        return parameter.get()

    def check_flowchart(self, preceding):
        """Check, before the job runs, that the MD will sample what is needed and is
        set up for the analysis.

        Parameters
        ----------
        preceding : [seamm.Node]
            The steps that run before this one.

        Returns
        -------
        [str]
            A description of each problem; empty if none.
        """
        source = self._static("md source")
        if source is None or self._static("approach") is None:
            return []
        return check_md_flowchart(self, preceding, self.md_requirements(), source)

    def check_configuration(self, configuration):
        """Called by an MD step before it runs MD that this step will analyze.

        Returns
        -------
        [str]
            Warnings about the structure, e.g. that it is small.

        Raises
        ------
        RuntimeError
            If the system is not periodic.
        """
        if configuration.periodicity != 3:
            raise RuntimeError(
                f"'{self.title}' will analyze the MD, but the system is not "
                "periodic: the diffusion coefficient of a cluster or molecule is not "
                "a bulk property, and the finite-size correction needs the cell."
            )
        return [f"'{self.title}': {w}" for w in size_warnings(configuration)]
