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
    eta, text = Node({"_viscosity": record})._viscosity_for_correction(P)
    assert eta["value"] == pytest.approx(8.5e-4)
    assert eta["stderr"] == pytest.approx(3e-5)
    assert eta["T"] == pytest.approx(300.0)
    assert eta["origin"] == "from step 4 (Viscosity, Green-Kubo, at 300.0 K)"
    assert (
        "0.85 ± 0.03 (standard error) mPa*s, from step 4 (Viscosity, Green-Kubo, "
        "at 300.0 K)" in text
    )


def test_no_viscosity_step_skips_the_correction():
    eta, text = Node({})._viscosity_for_correction(P)
    assert eta is None
    assert "skipped" in text


def test_given_viscosity():
    eta, text = Node({})._viscosity_for_correction(
        {**P, "viscosity source": "given value"}
    )
    assert eta["value"] == pytest.approx(5e-4)
    assert eta["stderr"] is None
    assert "given viscosity" in text
    eta, text = Node({})._viscosity_for_correction(
        {**P, "viscosity source": "given value", "viscosity": Q_(0.0, "cP")}
    )
    assert eta is None


def test_no_correction():
    assert Node({})._viscosity_for_correction(
        {**P, "hydrodynamic correction": False}
    ) == (None, None)


def test_correction_text_labels_every_error():
    correction = {
        "value": 0.6e-9,
        "stderr": 0.05e-9,
        "eta": 8.5e-4,
        "eta stderr": 3e-5,
        "origin": "from step 4 (Viscosity, Green-Kubo, at 300.0 K)",
        "T": 298.2,
        "L": 24.97,
        "varies": False,
    }
    text = Diffusivity._correction_text(correction, 1.0e-9, "10⁻⁹ m²/s")
    assert "ΔD = ξ kB T / (6π η L) = (0.6 ± 0.05 standard error) × 10⁻⁹ m²/s" in text
    assert "η = 0.85 ± 0.03 (standard error) mPa·s, from step 4" in text
    assert "L = 24.970 Å" in text
    assert "includes the error of η" in text


def test_old_flowcharts_still_read():
    """The settings for running the MD in the step are dropped; a step that ran its
    own MD says so, for the check of the flowchart."""
    import seamm

    if not hasattr(seamm.Parameters, "retired"):
        pytest.skip("needs seamm with Parameters.retired")
    assert diffusivity_step.DiffusivityParameters.retired == ("nruns", "errors")
    data = {
        "approach": {"value": "both", "units": None},
        "nruns": {"value": "$nruns", "units": None},
        "errors": {"value": "continue to next run", "units": None},
    }
    P = diffusivity_step.DiffusivityParameters(data=data)
    assert "nruns" not in P and "errors" not in P
    assert P["md source"].value == "run the MD in the sub-flowchart"
    data = {
        "md source": {
            "value": "the last MD run of each loop iteration",
            "units": None,
        },
        "nruns": {"value": "5", "units": None},
    }
    P = diffusivity_step.DiffusivityParameters(data=data)
    assert P["md source"].value == "the last MD run of each loop iteration"
