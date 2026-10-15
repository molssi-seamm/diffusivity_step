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
coefficients, removes this effect for a cubic cell. It needs the viscosity of the
fluid, which by default (**Viscosity from**: the preceding Viscosity step) comes
from the last Viscosity step run before this one; the output says which value was
used, and notes if it is for a different temperature. If there is no Viscosity step,
the correction is skipped with a note. A viscosity can also be given instead, and
**Correction for cell size** turns the correction off. Flowcharts written before the
Viscosity step existed keep the viscosity they give.

The step needs a periodic system: the diffusion coefficient of a molecule or cluster
in vacuum is not a bulk property, and the correction needs the cell. It stops before
running any MD if the structure is not periodic.

The step reads the positions and velocities from the ``md.h5md`` file that the LAMMPS
step writes for each MD run, using those of the molecules' centres of mass. Atomic
positions or velocities can be used only when every molecule is a single atom, as in
argon; otherwise the step stops with a message asking for the centres of mass. Runs
made with older versions of the LAMMPS step, which have no ``md.h5md``, are read
from their text files as before.


Where the MD comes from
=======================

**MD to analyze** chooses where the MD comes from:

run the MD in the sub-flowchart
    As before: the step runs its sub-flowchart, usually a LAMMPS step, the given
    number of times, and the errors come from the spread of the runs. Sampling that
    the analysis needs and the MD lacks is turned on in the last MD sub-step, with a
    note in the output.
the last MD run
    The last MD run earlier in the flowchart that sampled what the step needs. One
    production run can feed several analysis steps, and the sub-flowchart is not used.
the last MD run of each loop iteration
    For a LAMMPS step inside a Loop: the run of each iteration, analyzed as
    independent runs.
files
    The ``md.h5md`` files given, each an independent run; ``job://<n>/<path>`` reads a
    file of job ``n``.

A single run is divided into **Blocks for errors from one run** (5 by default),
analyzed separately to estimate the errors; each block must be much longer than the
correlation time. Before the job runs, the step checks that the MD before it samples
what it needs, and the job does not start if it does not.

The step needs the positions of the molecules' centres of mass sampled at least every 1
ps for the mean square displacement, and their velocities at least every 20 fs for the
Helfand moments.

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
