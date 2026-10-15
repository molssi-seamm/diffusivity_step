=======
History
=======
2026.10.14 -- Bugfix: crash with Helfand moments; per-direction coefficients
    * The step stopped with ``KeyError: 'Dx'`` after the second run whenever it used
      Helfand moments (the "Helfand moments" or "both" approach), while printing the
      current estimates. It now prints them.
    * The diffusion coefficients in x, y and z (Dx, Dy and Dz) were a third of their
      correct values, and had only a third of the Yeh–Hummer correction: they were
      the contributions of each direction to D rather than diffusion coefficients.
      For an isotropic liquid they now equal D. The total D is unchanged by this.
    * With Helfand moments, the total D is now computed from the sum of the moments
      in x, y and z, as with the MSD, rather than from the summed velocities, whose
      cross terms added noise. The value changes within its error.
    * The analysis now uses the shared ``seamm_md_analysis`` library, also used by the
      Thermal Conductivity step, instead of its own copy. The Helfand moments and mean
      square displacements are computed much faster (by FFT), which matters for long
      runs and many molecules.
    * The MSD graphs now label their y units as Å^2 rather than Å.

2026.10.12 -- Model chemistries
    * The results name the model chemistry (e.g. ``FF@oplsaa``) when LAMMPS uses one,
      rather than requiring a Forcefield step.

2026.3.26 -- Bugfix: Corrected handling of viscosity <= 0.0
    * Ignore the hydrostatic correction if the given viscosity is less than or equal
      to 0. This allows a flowchart to work more simply whether the viscosity os know or
      not.

2026.3.15 -- Fixed an error with MLFFs and protected the viscosity correction.
    * Fixed and error with the name of the MLFF, which caused the code to crash after
      printing the results but before ending.
    * Protected the hydrodynamic correction using the viscosity if the viscosity given
      is zero or less, in which case the correction is ignored.
      
2026.3.1 -- Internal: switching from deprecated library pkg_resources to importlib

2026.1.8 -- Enhancement to use atomic position and velocity trajectories
    * Add support for the atomic position and velocity trajectories for
      forcefields/potentials such as EAM/MEAM where there are no bonds, and hence the
      center-of-mass trajectories used from molecules are not allowed by LAMMPS.

2025.9.23 -- Bugfix: Fixed error storing diffusion constants as properties

2024.7.30.1 -- Bugfix: Reinitialization of data in loops
    * Fixed a bug where the data was not being correctly initialized if the step was in
      a loop. This caused the timing results to be incorrect.
    * Added actual number of samples for the MSD and length of the Helfand integral to
      the results.

2024.7.30 -- Added optional correction for cell size
    * Added an option to use the Yeh-Hummer hydrostatic correction for the effects of
      the finite cell size. The viscosity is required as an input, but the correction
      eliminates the need to extrapolate to 1/L = 0.
    * Added control parameters and timings to the available results.
      
2024.7.21 -- Significant improvements!
    * Simplified error analysis to safe approach of analyzing the diffusion constants
      over runs.
    * Improved fitting of the curves to focus on the central linear portion. There are
      reasonable defaults but the user can adjust as needed.
    * Provided a combined average and error bars when both the MSF approach and Helfand
      moments are used.
    * Capture temperature, pressure, and cell size from the MD step, providing 1/L as a
      result since the true diffusion constants are found by extrapolating to 1/L = 0.
    * Provided control over the number of steps for the expensive numerical integration
      in the Helfand moments, providing a reasonable default of 1000.
      
2024.7.15 -- Bugfix: Significant error in Helfand Moment approach
    * Now fixed and seems to be working.
      
2024.7.4 -- Improved fitting of curves
    * Removed weighting of the fit by the stdev since it is too biased to the beginning
    * Added control over the portion of the data to fit in order to avoid the initial
      curvature and poor data towards the end.
	
2024.6.3 -- Bugfix: handling of options for subflowchart
    * Fixed a bug where the options for the subflowchart were not being parsed
      correctly.

2024.5.26 -- Updated for new task handling
    * The new handling of running tasks such as LAMMPS required a small change in the
      code.
      
2023.9.5 -- Changed default to using only MSD
    * The Helfand moments approach seems give incorrect results if the sampling time is
      too long. It is not dramatic, but gives increasingly incorrect results as the
      sampling time is increased. Thus using the Helfand moments is dangerous because
      the results may be wrong, but not obviously so.

2023.8.30 -- Initial working version
    * A working version that has been tested somewhat. Further testing and documentation
      will follow

2023.5.8 -- Initial development version created
    * Plug-in created using the SEAMM plug-in cookiecutter.
