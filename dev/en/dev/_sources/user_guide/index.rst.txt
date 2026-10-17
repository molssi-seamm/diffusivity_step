.. _user-guide:

**********
User Guide
**********

The Diffusivity step calculates the self-diffusion coefficient of each species in a
fluid from molecular dynamics run before it. The species are the distinct molecules
(by their canonical SMILES); the motion of a molecule is that of its centre of mass.

Two methods are available, and by default both are used:

Mean square displacement (MSD)
    The Einstein relation: in each direction the MSD grows as 2 D t at long times.
    The slope is fit between the fractions of the curve given by **Start point of
    MSD fit** and **End point of MSD fit**. It needs the positions.
Helfand moments
    The square of the integral of the velocity of each molecule, which grows as
    2 D t in each direction; this is the Green-Kubo method in another form. The
    integral runs over at most **Maximum Helfand Integral length** frames -- the cost
    grows as the square of the number -- and the slope is fit between **Start point
    of Helfand fit** and **End point of Helfand fit**. It needs the velocities.

With both, "Combined" is the average of the two in each run, and is the value stored
as D.

Where the MD comes from
=======================

The step analyzes MD that ran before it in the flowchart, or the ``md.h5md`` files of
earlier runs. For independent runs, put the MD in a Loop -- in parallel, with each run
starting from the same structure ("Discard the structure" in the MD sub-steps) and new
velocities -- and analyze "the last MD run of each loop iteration", the default. A
typical flowchart:

.. code-block:: text

    Model Chemistry (FF@...), Packmol, minimization
    NPT equilibration; Set Cell to its mean volume; NVT equilibration
    Loop over the replicas (parallel)
        LAMMPS: velocities, NVT decorrelation, NVT production
    Viscosity
    Diffusivity

**MD to analyze** chooses the runs:

the last MD run of each loop iteration
    The runs of the MD step in the innermost Loop around it, as independent runs; with
    no Loop, the one run. The step must be after the Loop, not in it.
the last MD run
    The last MD run before the step that sampled what it needs.
files
    The ``md.h5md`` files given, each an independent run; ``job://<n>/<path>`` reads a
    file of job ``n``.

The step needs the positions sampled at least every 1 ps for the MSD and the
velocities at least every 20 fs for the Helfand moments -- of the molecules' centres
of mass, or of the atoms, from which the step finds the centres of mass (the atomic
positions are unwrapped, so the molecules are whole). Checking the flowchart before it
runs reports missing sampling, with a fix that turns it on, and also warns when the
analysis is inside the loop holding its MD, when the loop changes the temperature or
pressure (the runs are then not replicas of one state), and when the thermostat
perturbs the dynamics (Langevin, Andersen) or there is none (NVE): use Nose-Hoover or
CSVR for production. Before the MD runs, the step warns if the system is small, with
fewer than about 200 molecules or a cell narrower than three times the longest
molecule.

The step needs a periodic system: the diffusion coefficient of a molecule or cluster
in vacuum is not a bulk property, and the correction for the size of the cell needs
the cell. Self-diffusion averages over the molecules, so a modest box with several
runs of at least 500 ps each is better than one long run.

Errors and isotropy
===================

The error comes from the spread of the independent runs or, for a single run, of
**Blocks for errors from one run** (5 by default), each analyzed separately. The
tables give the 95% confidence interval -- Student's t for the number of runs or
blocks times the standard error -- and say so; the stored results are the standard
errors.

For an isotropic run -- a cubic cell with nothing directional in it -- D of each
species is one value, the average of the three directions in each run, and a line
checks that the directions agree within their errors (a χ² test). Otherwise the
components D_x, D_y and D_z are given too. The values of each run or block, with
their temperature and cell length, are in the step's own output (``step.out``), as
are the graphs of the MSD and Helfand moments (``MSD.graph``,
``HelfandMoments.graph``) with the fitted lines.

The size of the cell
====================

In a periodic cell the diffusion is slowed by the hydrodynamic interaction of a
molecule with its periodic images. The Yeh-Hummer correction,
ΔD = ξ k\ :sub:`B` T / (6π η L) with ξ = 2.837297 for a cubic cell of side L,
removes it. The table gives both D in the box and D corrected, and the text under it
gives ΔD, the temperature, the cell length and the viscosity η used, with its error
and where it came from. The error of D corrected includes that of η.

By default (**Viscosity from**: the preceding Viscosity step) η comes from the last
Viscosity step run before this one, which is the natural place for it: both analyze
the same runs. The output notes if that viscosity is for a different temperature. If
there is no Viscosity step, the correction is skipped with a note. A viscosity can
also be given instead, and **Correction for cell size** turns the correction off.

Results
=======

* D of each species (corrected for the cell size when the correction is applied) and
  its standard error, from the best method, and from each method (MSD, HM); the
  components D_x, D_y and D_z;
* D in the box, without the correction, and its standard error;
* the Yeh-Hummer correction with its standard error, the viscosity used, and 1/L,
  for extrapolating D against the size of the cell;
* the state of the runs: temperature, pressure, volume and density with their errors,
  and for runs at constant volume the equilibrium density at the target pressure when
  an NPT run before them gives the bulk modulus;
* the number of runs and frames, and the time the analysis took.

They can be stored as properties of the configuration, or in variables and tables.

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
