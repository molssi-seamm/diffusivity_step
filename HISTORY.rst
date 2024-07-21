=======
History
=======
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
