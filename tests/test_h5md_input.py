#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests that Diffusivity reads the molecules' centres of mass from md.h5md, from
the molecules' group or from the atoms'."""

import numpy as np
import pytest

from seamm_md_analysis import H5MDWriter, open_md

from diffusivity_step.diffusivity import Diffusivity

cell = np.diag([20.0, 20.0, 20.0])


def write(path, group, name, values, unit, every=10, masses=None, attrs=None):
    n = values.shape[1]
    with H5MDWriter(path) as w:
        if group == "all":
            w.particles(ids=list(range(1, n + 1)), masses=masses, cell=cell)
        else:
            w.particles(ids=[1, 2, 3, 4, 5, 6], cell=cell)
            w.particles(group=group, ids=list(range(n)), cell=cell)
        appender = w.particle_appender(name, n, unit, "fs", group=group, attrs=attrs)
        for i, frame in enumerate(values):
            appender.append(i * every, i * every * 0.5, frame)


def test_reads_molecular_positions(tmp_path):
    com = np.random.default_rng(1).random((5, 2, 3))
    path = tmp_path / "md.h5md"
    write(path, "molecules", "position", com * 10, "nm")
    with open_md(path) as md:
        values, dt = Diffusivity._molecular(md, "position", "Å", [[0, 1], [2, 3]])
    assert np.allclose(values, com * 100, rtol=1e-6)
    assert dt == pytest.approx(5.0)


def test_reads_velocities_in_angstrom_per_fs(tmp_path):
    path = tmp_path / "md.h5md"
    write(path, "molecules", "velocity", np.ones((4, 2, 3)), "Angstrom ps-1")
    with open_md(path) as md:
        values, _ = Diffusivity._molecular(md, "velocity", "Å/fs", [[0], [1]])
    assert np.allclose(values, 0.001)


def test_centres_of_mass_from_the_atoms(tmp_path):
    """Two diatomics, O-H like, from unwrapped atomic positions."""
    rng = np.random.default_rng(2)
    x = rng.random((7, 4, 3)) * 30  # unwrapped: may lie outside the cell
    masses = np.array([16.0, 1.0, 1.0, 16.0])
    molecules = [[0, 1], [3, 2]]
    path = tmp_path / "md.h5md"
    write(
        path,
        "all",
        "position",
        x,
        "Angstrom",
        masses=masses,
        attrs={"seamm_unwrapped": 1},
    )
    with open_md(path) as md:
        values, dt = Diffusivity._molecular(md, "position", "Å", molecules)
    x = x.astype(np.float32).astype(float)
    expected = np.stack(
        [
            (16 * x[:, 0] + x[:, 1]) / 17,
            (16 * x[:, 3] + x[:, 2]) / 17,
        ],
        axis=1,
    )
    assert values.shape == (7, 2, 3)
    assert np.allclose(values, expected, rtol=1e-6)
    assert dt == pytest.approx(5.0)


def test_single_atoms_are_molecules(tmp_path):
    x = np.random.default_rng(3).random((3, 3, 3))
    path = tmp_path / "md.h5md"
    write(
        path,
        "all",
        "velocity",
        x,
        "Angstrom fs-1",
        masses=[39.95] * 3,
    )
    with open_md(path) as md:
        values, _ = Diffusivity._molecular(md, "velocity", "Å/fs", [[2], [0], [1]])
    assert np.allclose(values, x[:, [2, 0, 1]], rtol=1e-6)


def test_wrapped_atomic_positions_are_refused(tmp_path):
    path = tmp_path / "md.h5md"
    write(path, "all", "position", np.zeros((3, 2, 3)), "Angstrom", masses=[1.0, 1.0])
    with open_md(path) as md:
        with pytest.raises(RuntimeError, match="not unwrapped"):
            Diffusivity._molecular(md, "position", "Å", [[0, 1]])


def test_the_atoms_must_match_the_structure(tmp_path):
    path = tmp_path / "md.h5md"
    write(
        path,
        "all",
        "position",
        np.zeros((3, 2, 3)),
        "Angstrom",
        masses=[1.0, 1.0],
        attrs={"seamm_unwrapped": 1},
    )
    with open_md(path) as md:
        with pytest.raises(RuntimeError, match="the structure has 3"):
            Diffusivity._molecular(md, "position", "Å", [[0, 1], [2]])


def test_requirements_follow_the_approach():
    d = Diffusivity.__new__(Diffusivity)
    req = d.md_requirements({"approach": "both"})
    assert [n.names[0] for n in req.needs] == [
        "molecules/position",
        "molecules/velocity",
    ]
    assert req.periodic
    assert req.dynamics
    req = d.md_requirements({"approach": "Mean Square Displacement (MSD)"})
    assert [n.names for n in req.needs] == [["molecules/position", "all/position"]]
    req = d.md_requirements({"approach": "Helfand moments"})
    assert [n.names[0] for n in req.needs] == ["molecules/velocity"]


def test_blocks_split_one_run():
    values = np.arange(100 * 2 * 3, dtype=float).reshape(100, 2, 3)
    run = {"position": values, "position dt": 10.0, "T": 300.0, "L": 20.0}
    blocks = Diffusivity._split_into_blocks(run, 4)
    assert len(blocks) == 4
    assert [len(b["position"]) for b in blocks] == [25] * 4
    assert np.array_equal(blocks[1]["position"], values[25:50])
    assert blocks[3]["T"] == 300.0
    assert blocks[3]["position dt"] == 10.0
    with pytest.raises(RuntimeError, match="too few"):
        Diffusivity._split_into_blocks(run, 20)
