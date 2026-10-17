# -*- coding: utf-8 -*-

"""Smoke test of the Tk dialog: create it and re-lay it out for every combination
of the choices that change it. Skipped when no display is available."""

import itertools

import pytest

from seamm_md_analysis.step_support import source_choices

from diffusivity_step import DiffusivityParameters


@pytest.fixture()
def tk_node():
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available for Tk")
    root.withdraw()
    import Pmw
    import seamm

    Pmw.initialise(root)
    flowchart = seamm.Flowchart(namespace="org.molssi.seamm", directory=".")
    tk_flowchart = seamm.TkFlowchart(
        master=root, flowchart=flowchart, namespace="org.molssi.seamm.tk"
    )
    node = flowchart.create_node("Diffusivity")
    flowchart.add_node(node)
    plugin = tk_flowchart.plugin_manager.get("Diffusivity")
    tk_node = plugin.create_tk_node(
        tk_flowchart=tk_flowchart, node=node, canvas=tk_flowchart.canvas, x=100, y=100
    )
    yield tk_node
    root.destroy()


def test_dialog_layouts(tk_node):
    tk_node.create_dialog()
    tk_node.reset_dialog()
    P = tk_node.node.parameters
    choices = {
        "md source": source_choices,
        "approach": P["approach"].enumeration,
        "hydrodynamic correction": ("no", "yes"),
        "viscosity source": P["viscosity source"].enumeration,
    }
    keys = [k for k in DiffusivityParameters.parameters if k != "results"]
    for combination in itertools.product(*choices.values()):
        settings = dict(zip(choices, combination))
        for key, value in settings.items():
            tk_node[key].set(value)
        tk_node.reset_dialog()
        values = {k: tk_node[k].get() for k in keys}
        for key in keys:
            shown = tk_node[key].grid_info() != {}
            assert shown == P.applies(key, values), (key, settings)


def test_no_own_md(tk_node):
    """The dialog has no sub-flowchart and no settings for running MD."""
    tk_node.create_dialog()
    notebook = tk_node["notebook"]
    tabs = [notebook.tab(t, "text") for t in notebook.tabs()]
    assert tabs == ["Parameters", "Results"]
    for key in ("nruns", "errors"):
        assert key not in DiffusivityParameters.parameters


def test_viscosity_layouts(tk_node):
    tk_node.create_dialog()
    shown = {}
    for correction in ("no", "yes"):
        for source in ("the preceding Viscosity step", "given value"):
            tk_node["hydrodynamic correction"].set(correction)
            tk_node["viscosity source"].set(source)
            tk_node.reset_dialog()
            shown[(correction, source)] = (
                tk_node["viscosity source"].grid_info() != {},
                tk_node["viscosity"].grid_info() != {},
            )
    assert shown[("no", "the preceding Viscosity step")] == (False, False)
    assert shown[("no", "given value")] == (False, False)
    assert shown[("yes", "the preceding Viscosity step")] == (True, False)
    assert shown[("yes", "given value")] == (True, True)
