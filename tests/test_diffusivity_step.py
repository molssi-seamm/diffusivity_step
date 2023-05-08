#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `diffusivity_step` package."""

import pytest  # noqa: F401
import diffusivity_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = diffusivity_step.Diffusivity()
    assert str(type(result)) == "<class 'diffusivity_step.diffusivity.Diffusivity'>"
