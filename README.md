Input files for the labs I did.
Sorted into each lab.
For the first lab there are two files for the DFT theory with the basis set cc-pVTZ to show the difference in resource allocation that sped up SCF convergence/run time.
As mentioned in report also, the npt.mdp timestep was increased to 500 ps to allow a closer relaxation so that the pressure could be nearer to one bar rather than the consistent 8 bar with a 200 ps timestep. 
using   gmx energy -f md.edr -o md_pressure.xvg
the average pressure was:
Energy                      Average   Err.Est.       RMSD  Tot-Drift
-------------------------------------------------------------------------------
Pressure                   -2.95914        5.6    603.848   -11.5074  (bar)
which isn't as close to 1 bar, but closer than the 8 bar consistently scored. The timestep should be increased to the ns scale.
