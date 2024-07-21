# -*- coding: utf-8 -*-

"""Routines to help do the analysis for diffusivity."""
import json

# import warnings

import numpy as np
from scipy.optimize import curve_fit  # , OptimizeWarning
from scipy.integrate import cumulative_trapezoid

# import statsmodels.tsa.stattools as stattools

tensor_labels = [
    ("x", "red", "rgba(255,0,0,0.1)"),
    ("y", "green", "rgba(0,255,0,0.1)"),
    ("z", "blue", "rgba(0,0,255,0.1)"),
    ("", "black", "rgba(0,0,0,0.1)"),
]


def axb(x, a, b):
    return a * x + b


def compute_msd(xyz, species, average=True):
    """Compute the mean square displacement vector, averaging over molecules.

    Parameters
    ----------
    xyz : numpy.ndarray(nsteps, nmolecules, 3)
        The coordinates at each step
    species : [[int]] or numpy.ndarray
        Indices of the molecules for each species
    average : bool = True
        Whether to return the average MSD as well as x, y, and z components.

    Returns
    -------
    [numpy.ndarray(nsteps, 3 or 4)] * nspecies
        The MSD plus the average MSD if <average> == True, as a list over species
    [numpy.ndarray(nsteps, 3 or 4)] * nspeceis
        The error of the MSD
    """

    nsteps, nmolecules, _ = xyz.shape

    tmp = np.zeros_like(xyz)
    for i in range(nsteps):
        j = nsteps - i
        tmp[:j] += (xyz[i:] - xyz[i]) ** 2

    msd = []
    err = []
    for i, molecules in enumerate(species):
        msd.append(np.average(tmp[:, molecules, :], axis=1))
        err.append(np.std(tmp[:, molecules, :], axis=1))

        if average:
            tmp_ave = np.sum(tmp, axis=2)
            msd_ave = np.average(tmp_ave[:, molecules], axis=1)
            err_ave = np.std(tmp_ave[:, molecules], axis=1)
            msd_ave = msd_ave.reshape((nsteps, 1))
            err_ave = err_ave.reshape((nsteps, 1))
            msd[i] = np.concatenate((msd[i], msd_ave), axis=1)
            err[i] = np.concatenate((err[i], err_ave), axis=1)

        # Normalize for origins
        norm = np.arange(nsteps, 0, -1)
        norm = norm.reshape((nsteps, 1))
        msd[i] /= norm
        err[i] /= norm

    return msd, err


def create_helfand_moments(v, species, m=None):
    """Create the Helfand moments from velocities

    Parameters
    ----------
    v : numpy.ndarray(nframes, nmols, 3)
        The heat fluxes in x, y, and z
    species : [[int]] or numpy.ndarray
        Indices of the molecules for each species
    m : int
        The length of the Helfand moments wanted

    Returns
    -------
    [numpy.ndarray(m, 3)] * nspecies
        The Helfand moments
    """

    n, nmols, _ = v.shape
    if m is None:
        m = min(n // 2, 5000)

    Ms = []
    M_errs = []
    for i, molecules in enumerate(species):
        M = np.zeros((m, 4))
        M_err = np.zeros((m, 4))

        for alpha in range(3):
            tmp = np.zeros((m, len(molecules)))
            for i in range(n - m):
                integral = cumulative_trapezoid(
                    v[i : m + i, molecules, alpha], initial=0.0, axis=0
                )
                tmp += integral * integral
            M[:, alpha] = np.average(tmp, axis=1)
            M_err[:, alpha] = np.std(tmp, axis=1)

        # and sum
        vsum = v[:, molecules, 0] + v[:, molecules, 1] + v[:, molecules, 2]
        tmp = np.zeros((m, len(molecules)))
        for i in range(n - m):
            integral = cumulative_trapezoid(vsum[i : m + i, :], initial=0.0, axis=0)
            tmp += integral * integral
        M[:, 3] = np.average(tmp, axis=1)
        M_err[:, 3] = np.std(tmp, axis=1)

        M /= n - m
        M_err /= n - m

        Ms.append(M)
        M_errs.append(M_err)
    return Ms, M_errs


def fit_slope(y, xs, sigma=None, start=0.1, end=0.9):
    """Find the best linear fit to longest possible segment.

    Parameters
    ----------
    y : [float] or numpy.ndarray()
        The MSD

    xs : [float]
        The time (x) coordinate

    sigma : [float] or numpy.ndarray()
        Optional standard error of y

    start : float
        Fraction of vector to ignore at the beginning

    end : float
        The fraction of vector to end at

    Returns
    -------
    slope : float
        The fit slope.
    stderr : float
        The 95% standard error of the slope
    xs : [float]
        The x values (time) for the fit curve
    ys : [float]
        The y values for the fit curve.
    """
    # We know the curves curve near the origin, so ignore the first part and last
    n = len(y)
    i = int(n * start)
    j = int(n * end)

    if sigma is None:
        popt, pcov, infodict, msg, ierr = curve_fit(
            axb,
            xs[i:j],
            y[i:j],
            full_output=True,
        )
    else:
        popt, pcov, infodict, msg, ierr = curve_fit(
            axb,
            xs[i:j],
            y[i:j],
            full_output=True,
            sigma=sigma[i:],
            absolute_sigma=True,
        )
    slope = float(popt[0])
    b = float(popt[1])
    err = float(np.sqrt(np.diag(pcov)[0]))

    ys = []
    for x in xs[i:j]:
        ys.append(axb(x, slope, b))

    return slope, err, xs[i:j], ys


def add_helfand_trace(
    plot, x_axis, y_axis, species, M, ts, err=None, fit=None, labels=tensor_labels
):
    """Add a trace for the Helfand moments.

    Parameters
    ----------
    plot : seamm_util.plot
        The plot that contains the traces

    x_axis :
        The x axis for the plot

    y_axis :
        The y axis for the plot

    species : str
        The label for the species

    M : numpy.mdarray(m, 3)
        The Helfand moments, in m^2

    ts : [float]
        The times associated with the moments, in ps

    err : numpy.mdarray(n, 3 or 4)
        The std error on the moments

    fit : {str: any}
        The information for the fit to the curve

    labels : {str: xxx}
        The labels for the directions
    """
    nsteps, nalpha = M.shape

    for i in range(nalpha):
        label, color, colora = tensor_labels[i]
        if fit is not None:
            hover = (
                f"{species} D{label} = {fit[i]['D_s']} * "
                f"{fit[i]['scale']:.1e} m^2/s"
            )
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"{species} fit{label}",
                hovertemplate=hover,
                x=fit[i]["xs"],
                xlabel="t",
                xunits="ps",
                y=fit[i]["ys"],
                ylabel=f"{species} fit{label}",
                yunits="Å^2",
                color=color,
                dash="dash",
                width=3,
            )
        if err is not None:
            errs = np.concatenate((M[:, i] + err[:, i], M[::-1, i] - err[::-1, i]))
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"{species} ±{label}",
                x=ts + ts[::-1],
                xlabel="t",
                xunits="ps",
                y=errs.tolist(),
                ylabel=f"{species} ±{label}",
                yunits="Å^2",
                color=colora,
                fill="toself",
                visible="legendonly",
            )
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name=f"{species} M{label}",
            x=ts,
            xlabel="t",
            xunits="ps",
            y=M[:, i].tolist(),
            ylabel=f"{species} M{label}",
            yunits="Å^2",
            color=color,
        )


def plot_helfand_moments(figure, M, ts, err=None, fit=None, labels=tensor_labels):
    """Create a plot for the Helfand moments.

    Parameters
    ----------
    figure : seamm_util.Figure
        The figure that contains the plots.

    M : numpy.mdarray(m, 3)
        The Helfand moments, in m^2

    ts : [float]
        The times associated with the moments, in ps
    """
    nsteps, nalpha = M.shape

    plot = figure.add_plot("HelfandMoments")

    x_axis = plot.add_axis("x", label="Time (ps)")
    y_axis = plot.add_axis("y", label="M (Å^2)", anchor=x_axis)
    x_axis.anchor = y_axis

    for i in range(nalpha):
        label, color, colora = tensor_labels[i]
        if fit is not None:
            hover = (
                f"D{label} = ({fit[i]['D_s']} ± {fit[i]['err_s']}) * "
                f"{fit[i]['scale']:.1e} m^2/s"
            )
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"fit{label}",
                hovertemplate=hover,
                x=fit[i]["xs"],
                xlabel="t",
                xunits="ps",
                y=fit[i]["ys"],
                ylabel=f"fit{label}",
                yunits="Å^2",
                color=color,
                dash="dash",
                width=3,
            )
        if err is not None:
            errs = np.concatenate((M[:, i] + err[:, i], M[::-1, i] - err[::-1, i]))
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"±{label}",
                x=ts + ts[::-1],
                xlabel="t",
                xunits="ps",
                y=errs.tolist(),
                ylabel=f"±{label}",
                yunits="Å^2",
                color=colora,
                fill="toself",
                visible="legendonly",
            )
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name=f"M{label}",
            x=ts,
            xlabel="t",
            xunits="ps",
            y=M[:, i].tolist(),
            ylabel=f"M{label}",
            yunits="Å^2",
            color=color,
        )
    return plot


def add_msd_trace(
    plot,
    x_axis,
    y_axis,
    species,
    msd,
    ts,
    err=None,
    fit=None,
    labels=tensor_labels,
):
    """Add a trace to the mean square deviation (MSD) plot.

    Parameters
    ----------
    plot : seamm_util.plot
        The plot that contains the traces

    x_axis :
        The x axis for the plot

    y_axis :
        The y axis for the plot

    species : str
        The label for the species

    msd : numpy.mdarray(n, 3 or 4)
        The MSD in x, y, and z, and optionally the average

    ts : [float]
        The times associated with the MSD, in ps

    err : numpy.mdarray(n, 3 or 4)
        The std error on the msd values

    fit : {str: any}
        The information for the fit to the curve

    labels : {str: xxx}
        The labels for the directions
    """
    nsteps, nalpha = msd.shape

    for i in range(nalpha):
        label, color, colora = tensor_labels[i]
        if fit is not None:
            hover = (
                f"{species} D{label} = {fit[i]['D_s']} * "
                f"{fit[i]['scale']:.1e} m^2/s"
            )
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"{species} fit{label}",
                hovertemplate=hover,
                x=fit[i]["xs"],
                xlabel="t",
                xunits="ps",
                y=fit[i]["ys"],
                ylabel=f"{species} fit{label}",
                yunits="Å^2",
                color=color,
                dash="dash",
                width=3,
            )
        if err is not None:
            errs = np.concatenate((msd[:, i] + err[:, i], msd[::-1, i] - err[::-1, i]))
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"{species} ±{label}",
                x=ts + ts[::-1],
                xlabel="t",
                xunits="ps",
                y=errs.tolist(),
                ylabel=f"{species} ±{label}",
                yunits="Å^2",
                color=colora,
                fill="toself",
                visible="legendonly",
            )
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name=f"{species} MSD{label}",
            x=ts,
            xlabel="t",
            xunits="ps",
            y=msd[:, i].tolist(),
            ylabel=f"{species} MSD{label}",
            yunits="Å",
            color=color,
        )


def plot_msd(figure, msd, ts, err=None, fit=None, labels=tensor_labels):
    """Create a plot for the mean square deviation (MSD).

    Parameters
    ----------
    figure : seamm_util.Figure
        The figure that contains the plots.

    msd : numpy.mdarray(n, 3 or 4)
        The MSD in x, y, and z, and optionally the average

    ts : [float]
        The times associated with the MSD, in ps
    """
    nsteps, nalpha = msd.shape

    plot = figure.add_plot("MSD")

    x_axis = plot.add_axis("x", label="Time (ps)")
    y_axis = plot.add_axis("y", label="MSD (Å^2)", anchor=x_axis)
    x_axis.anchor = y_axis

    nsteps, nalpha = msd.shape

    for i in range(nalpha):
        label, color, colora = tensor_labels[i]
        if fit is not None:
            hover = (
                f"D{label} = ({fit[i]['D_s']} ± {fit[i]['err_s']}) * "
                f"{fit[i]['scale']:.1e} m^2/s"
            )
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"fit{label}",
                hovertemplate=hover,
                x=fit[i]["xs"],
                xlabel="t",
                xunits="ps",
                y=fit[i]["ys"],
                ylabel=f"fit{label}",
                yunits="Å^2",
                color=color,
                dash="dash",
                width=3,
            )
        if err is not None:
            errs = np.concatenate((msd[:, i] + err[:, i], msd[::-1, i] - err[::-1, i]))
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"±{label}",
                x=ts + ts[::-1],
                xlabel="t",
                xunits="ps",
                y=errs.tolist(),
                ylabel=f"±{label}",
                yunits="Å^2",
                color=colora,
                fill="toself",
                visible="legendonly",
            )
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name=f"MSD{label}",
            x=ts,
            xlabel="t",
            xunits="ps",
            y=msd[:, i].tolist(),
            ylabel=f"MSD{label}",
            yunits="Å",
            color=color,
        )
    return plot


def read_vector_trajectory(path):
    """Read a standard vector trajectory from LAMMPS fix vector.

    Parameters
    ----------
    path : pathlib.Path
        The file as a string or path-like object.
    """
    data = []
    with open(path) as fd:
        lines = iter(fd)

        # Get the initial header line and check
        line = next(lines)
        tmp = line.split()
        if len(tmp) < 3:
            raise RuntimeError(f"Bad header for {path}: {line}")
        if tmp[0] != "!MolSSI":
            raise RuntimeError(f"Not a MolSSI file? {path}: {line}")
        if tmp[1] != "vector_trajectory":
            raise RuntimeError(f"Not a vector_trajectory file? {path}: {line}")
        if tmp[2][0] != "2":
            raise RuntimeError(
                f"Can only handle version 2 vector_trajectory files. {path}: {line}"
            )
        metadata = json.loads(" ".join(tmp[3:]))

        # Skip any commented header lines
        for line in lines:
            if line[0] != "!":
                break
        metadata["fields"] = line.split()[1:]

        for line in lines:
            _, n = line.split()
            frame = []
            for i in range(int(n)):
                frame.append([float(v) for v in next(lines).split()[1:]])
            data.append(frame)

    result = np.array(data)
    return metadata, result
