.. _user-guide:

**********
User Guide
**********

The diffusivity step in SEAMM calculates the diffusion coefficients using molecular
dynamics. You will need to setup an appropriate molecular dynamics run in the
subflowchart and ensure that it produces the center-of-mass positions and/or velocities
of the particles/molecules.

The positions are used to calculate the diffusion coefficient from the
mean-sqaure-displacement (MSD) of the molecules. The velocities are needed for the
approach using Helfand moments, which is related to the more traditional Green-Kubo
method.

..
    <remove the dots above and this line and unindent the toctree to expose it>
    Contents:

    .. toctree::
       :glob:
       :maxdepth: 2
       :titlesonly:

       *

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
