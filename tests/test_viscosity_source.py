#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests of where Diffusivity gets the viscosity for the Yeh-Hummer correction."""

import pytest

from seamm_util import Q_

import diffusivity_step
from diffusivity_step.diffusivity import Diffusivity


def test_old_flowcharts_keep_their_given_viscosity():
    data = {
        "hydrodynamic correction": {"value": "yes", "units": None},
        "viscosity": {"value": "0.89", "units": "cP"},
    }
    P = diffusivity_step.DiffusivityParameters(data=data)
    assert P["viscosity source"].value == "given value"
    assert P["viscosity"].value == "0.89"


def test_new_steps_use_the_viscosity_step():
    P = diffusivity_step.DiffusivityParameters()
    assert P["hydrodynamic correction"].value == "yes"
    assert P["viscosity source"].value == "the preceding Viscosity step"
    assert not P.applies("viscosity")


class Node(Diffusivity):
    def __init__(self, variables):
        self.variables = variables

    def variable_exists(self, name):
        return name in self.variables

    def get_variable(self, name):
        return self.variables[name]


P = {
    "hydrodynamic correction": True,
    "viscosity source": "the preceding Viscosity step",
    "viscosity": Q_(0.5, "cP"),
}


def test_viscosity_from_the_viscosity_step():
    record = {
        "value": 8.5e-4,
        "stderr": 3e-5,
        "units": "Pa*s",
        "T": 300.0,
        "method": "Green-Kubo",
        "step": "Viscosity",
        "id": "4",
    }
    eta, T, text = Node({"_viscosity": record})._viscosity_for_correction(P)
    assert eta == pytest.approx(8.5e-4)
    assert T == pytest.approx(300.0)
    assert "0.85 ± 0.03 mPa*s, from step 4 (Viscosity, Green-Kubo, at 300.0 K)" in text


def test_no_viscosity_step_skips_the_correction():
    eta, T, text = Node({})._viscosity_for_correction(P)
    assert eta is None
    assert "skipped" in text


def test_given_viscosity():
    eta, _, text = Node({})._viscosity_for_correction(
        {**P, "viscosity source": "given value"}
    )
    assert eta == pytest.approx(5e-4)
    assert "given viscosity" in text
    eta, _, text = Node({})._viscosity_for_correction(
        {**P, "viscosity source": "given value", "viscosity": Q_(0.0, "cP")}
    )
    assert eta is None


def test_no_correction():
    assert Node({})._viscosity_for_correction(
        {**P, "hydrodynamic correction": False}
    ) == (None, None, None)
