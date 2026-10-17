#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Regression test of the analysis on Brownian motion with a known D: the values,
the errors, the Yeh-Hummer correction and the stored results."""

import numpy as np
import pytest

from seamm_md_analysis import H5MDWriter, yeh_hummer_correction
from seamm_util import Q_

D_TRUE = 2.0e-9  # m^2/s
L = 20.0  # Å
N_ATOMS = 64
N_FRAMES = 400
DT = 500.0  # fs between frames
T = 300.0

P = {
    "md source": "files",
    "approach": "Mean Square Displacement (MSD)",
    "number of blocks": 4,
    "msd_fit_start": 0.1,
    "msd_fit_end": 0.6,
    "helfand_fit_start": 0.35,
    "helfand_fit_end": 1.0,
    "maximum Helfand Integral length": 1000,
}


def brownian_run(path, seed, cell=(L, L, L), velocities=False):
    """Write a run of independent random walks of argon atoms."""
    rng = np.random.default_rng(seed)
    # D in Å^2/fs; each step is Gaussian with variance 2 D dt in each direction
    sigma = np.sqrt(2 * D_TRUE * 1.0e5 * DT)
    steps = rng.normal(0.0, sigma, size=(N_FRAMES, N_ATOMS, 3))
    steps[0] = rng.random((N_ATOMS, 3)) * L
    x = np.cumsum(steps, axis=0)
    with H5MDWriter(path) as w:
        w.particles(
            ids=list(range(1, N_ATOMS + 1)),
            species=[18] * N_ATOMS,
            masses=[39.948] * N_ATOMS,
            cell=np.diag(cell),
        )
        appender = w.particle_appender(
            "position", N_ATOMS, "Angstrom", "fs", attrs={"seamm_unwrapped": 1}
        )
        for i, frame in enumerate(x):
            appender.append(i * 250, i * DT, frame)
        if velocities:
            # The mean velocity over each interval, so that it integrates to x
            v = np.concatenate([steps[1:], steps[-1:]]) / DT
            appender = w.particle_appender("velocity", N_ATOMS, "Angstrom fs-1", "fs")
            for i, frame in enumerate(v):
                appender.append(i * 250, i * DT, frame)
        w.observable(
            "temperature",
            np.arange(N_FRAMES) * 250,
            np.arange(N_FRAMES) * DT,
            np.full(N_FRAMES, T),
            "K",
        )
        w.parameters({"ensemble": "NVT", "temperature": T, "periodic": True})


@pytest.fixture()
def configuration():
    import molsystem

    db = molsystem.SystemDB(filename="file:seamm_db?mode=memory&cache=shared")
    system = db.create_system()
    configuration = system.create_configuration(
        periodicity=3, coordinate_system="Cartesian"
    )
    configuration.cell.parameters = (L, L, L, 90.0, 90.0, 90.0)
    xyz = np.random.default_rng(0).random((N_ATOMS, 3)) * L
    configuration.atoms.append(
        x=xyz[:, 0].tolist(),
        y=xyz[:, 1].tolist(),
        z=xyz[:, 2].tolist(),
        symbol=["Ar"] * N_ATOMS,
    )
    yield configuration
    db.close()


@pytest.fixture()
def node(tmp_path, monkeypatch):
    import seamm
    from seamm.variables import Variables

    monkeypatch.setattr(seamm, "flowchart_variables", Variables())
    flowchart = seamm.Flowchart(namespace="org.molssi.seamm", directory=str(tmp_path))
    node = flowchart.create_node("Diffusivity")
    flowchart.add_node(node)
    node._id = ("1",)
    (tmp_path / "1").mkdir()
    return node


def analyze(node, configuration, paths, viscosity=None):
    node.paths = paths
    node._viscosity = viscosity
    return node.analyze(P=dict(P), configuration=configuration)


def test_runs(node, configuration, tmp_path):
    paths = []
    for i in range(3):
        path = tmp_path / f"run_{i}.h5md"
        brownian_run(path, seed=i)
        paths.append(path)
    data = analyze(node, configuration, paths)

    assert data["nruns"] == 3
    D = data["D {key}"]["[Ar]"]
    stderr = data["D {key},stderr"]["[Ar]"]
    assert D == pytest.approx(D_TRUE, rel=0.1)
    assert 0 < stderr < 0.1 * D_TRUE
    # The best estimate is the MSD's, and without the correction D is in the box
    assert data["D {key} (MSD)"]["[Ar]"] == D
    assert data["D (box) {key}"]["[Ar]"] == D
    assert data["D (box) {key},stderr"]["[Ar]"] == stderr
    # D is the average of the components
    components = [data[f"D{c} {{key}}"]["[Ar]"] for c in "xyz"]
    assert np.mean(components) == pytest.approx(D)
    assert "Yeh-Hummer correction" not in data
    assert data["1/L"] == pytest.approx(1 / L)
    assert data["T"] == pytest.approx(T)
    assert (tmp_path / "1" / "MSD.graph").exists()


def test_one_run_in_blocks_with_the_correction(node, configuration, tmp_path):
    path = tmp_path / "run.h5md"
    brownian_run(path, seed=7)
    eta = 8.5e-4
    viscosity = {
        "value": eta,
        "stderr": 0.3e-4,
        "T": 300.0,
        "origin": "from step 4 (Viscosity, Green-Kubo, at 300.0 K)",
    }
    data = analyze(node, configuration, [path], viscosity)

    assert data["nruns"] == 1
    assert data["msd samples"] == N_FRAMES // 4
    correction = yeh_hummer_correction(T, L * 1.0e-10, eta)
    assert data["Yeh-Hummer correction"] == pytest.approx(correction)
    assert data["Yeh-Hummer correction,stderr"] == pytest.approx(
        correction * 0.3e-4 / eta
    )
    assert data["viscosity for correction"] == pytest.approx(
        Q_(eta, "Pa*s").m_as("mPa*s")
    )
    box = data["D (box) {key}"]["[Ar]"]
    D = data["D {key}"]["[Ar]"]
    assert D == pytest.approx(box + correction)
    assert box == pytest.approx(D_TRUE, rel=0.15)
    # The error of the corrected D includes that of the viscosity
    box_error = data["D (box) {key},stderr"]["[Ar]"]
    error = data["D {key},stderr"]["[Ar]"]
    assert error == pytest.approx(
        np.sqrt(box_error**2 + data["Yeh-Hummer correction,stderr"] ** 2)
    )


def test_both_methods_in_a_cell_that_is_not_cubic(node, configuration, tmp_path):
    paths = []
    for i in range(3):
        path = tmp_path / f"run_{i}.h5md"
        brownian_run(path, seed=10 + i, cell=(L, L, 1.2 * L), velocities=True)
        paths.append(path)
    node.paths = paths
    node._viscosity = None
    data = node.analyze(
        P={**P, "approach": "both", "maximum Helfand Integral length": 100},
        configuration=configuration,
    )
    msd = data["D {key} (MSD)"]["[Ar]"]
    hm = data["D {key} (HM)"]["[Ar]"]
    assert msd == pytest.approx(D_TRUE, rel=0.15)
    assert hm == pytest.approx(D_TRUE, rel=0.15)
    # The best estimate is the average of the two methods
    assert data["D {key}"]["[Ar]"] == pytest.approx((msd + hm) / 2)
    assert data["Helfand integral length"] == 100
    assert (tmp_path / "1" / "HelfandMoments.graph").exists()
