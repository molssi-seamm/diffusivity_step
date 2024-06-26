=======
History
=======
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
