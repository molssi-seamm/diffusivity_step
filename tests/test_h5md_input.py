#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests that Diffusivity reads md.h5md, and falls back to older files."""

import numpy as np
import pytest

from seamm_md_analysis import H5MDWriter

from diffusivity_step.diffusivity import Diffusivity

cell = np.diag([20.0, 20.0, 20.0])


def write(path, group, name, values, unit, every=10):
    with H5MDWriter(path) as w:
        w.particles(ids=[1, 2, 3, 4, 5, 6], cell=cell)
        if group != "all":
            w.particles(group=group, ids=list(range(values.shape[1])), cell=cell)
        appender = w.particle_appender(name, values.shape[1], unit, "fs", group=group)
        for i, frame in enumerate(values):
            appender.append(i * every, i * every * 0.5, frame)


def step(n_molecules=2):
    d = Diffusivity.__new__(Diffusivity)
    d._n_molecules = n_molecules
    return d


def test_reads_molecular_positions(tmp_path):
    run_dir = tmp_path / "run_1"
    (run_dir / "1" / "3").mkdir(parents=True)
    com = np.random.default_rng(1).random((5, 2, 3))
    write(run_dir / "1" / "3" / "md.h5md", "molecules", "position", com * 10, "nm")
    series = step()._read_h5md(1, run_dir, "position", "Å")
    assert np.allclose(series.values, com * 100, rtol=1e-6)
    assert series.dt.m_as("fs") == pytest.approx(5.0)


def test_reads_velocities_in_angstrom_per_fs(tmp_path):
    run_dir = tmp_path / "run_1"
    run_dir.mkdir()
    v = np.ones((4, 2, 3))
    write(run_dir / "md.h5md", "molecules", "velocity", v, "Angstrom ps-1")
    series = step()._read_h5md(1, run_dir, "velocity", "Å/fs")
    assert np.allclose(series.values, 0.001)


def test_no_h5md_returns_none(tmp_path):
    assert step()._read_h5md(1, tmp_path, "position", "Å") is None


def test_two_runs_with_data_are_refused(tmp_path):
    for sub in ("a", "b"):
        (tmp_path / sub).mkdir()
        write(
            tmp_path / sub / "md.h5md",
            "molecules",
            "position",
            np.zeros((2, 2, 3)),
            "Angstrom",
        )
    with pytest.raises(NotImplementedError):
        step()._read_h5md(1, tmp_path, "position", "Å")


def test_atomic_data_must_be_per_molecule():
    d = step(n_molecules=2)
    d._check_n_molecules(np.zeros((3, 2, 3)), 1, "positions")
    with pytest.raises(RuntimeError, match="centers of mass"):
        d._check_n_molecules(np.zeros((3, 6, 3)), 1, "positions")
