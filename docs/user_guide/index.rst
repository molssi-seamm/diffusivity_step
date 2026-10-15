.. _user-guide:

**********
User Guide
**********

The diffusivity step in SEAMM calculates the diffusion coefficients using molecular
dynamics. You will need to setup an appropriate molecular dynamics run in the
subflowchart and ensure that it produces the center-of-mass positions and/or velocities
of the particles/molecules.

The positions are used to calculate the diffusion coefficient from the
mean-square-displacement (MSD) of the molecules. The velocities are needed for the
approach using Helfand moments, which is related to the more traditional Green-Kubo
method.

The step reports the diffusion coefficients in each direction, Dx, Dy and Dz, as
well as the total D, their average. For an isotropic fluid the four agree within
their errors. The errors are 95% intervals (twice the standard error) from the
spread of the values over the runs, so at least two runs are needed.

In a periodic cell the diffusion is slowed by the hydrodynamic interaction of a
molecule with its periodic images. The Yeh-Hummer correction, added to each of the
coefficients when requested, removes this effect for a cubic cell; it needs the
viscosity of the fluid.

The step needs a periodic system: the diffusion coefficient of a molecule or cluster
in vacuum is not a bulk property, and the correction needs the cell. It stops before
running any MD if the structure is not periodic.

The step reads the positions and velocities from the ``md.h5md`` file that the LAMMPS
step writes for each MD run, using those of the molecules' centres of mass. Atomic
positions or velocities can be used only when every molecule is a single atom, as in
argon; otherwise the step stops with a message asking for the centres of mass. Runs
made with older versions of the LAMMPS step, which have no ``md.h5md``, are read
from their text files as before.

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
