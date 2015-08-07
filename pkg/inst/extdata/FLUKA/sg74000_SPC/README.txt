====================================================================
sg74000 - FLUKA/FLAIR files and R script to produce SPC and DDD data
====================================================================

0. General remarks
==================
Geometry according to Parodi et al., 2012, [no detailed BAMS simulation], not to the actually geometry (scattering is to some extend accounted for by having the RIFI 100 mm upstream). Thus, the actual isocenter position is (no matter the RIFI) at 2.89 mm depth!

Compile the executable using ldpm3qmd and the user mgdraw routine.


1. First run
============
Leave 'PeakPos' variable set in DEFINE card at some value, e.g. 1.0
Deactivate USERDUMP card
Set primary energy to desired value (negative!)
Check if ripple filter is turned on (or off, as desired), material for CELLxx should be PMMA (or air if off)
Run few histories to find approximate peak position (using detector DOSEGRID)


2. Second run
=============
Set the DOSEHRES detector to cover the approximate peak position with steps of 10 µm.
Re-run with somewhat higher number of primaries.


3. Third run
============
Set 'PeakPos' variable to the peak position found. 
Enable USERDUMP card.
Run with sufficient primaries (e.g. 5 cycles with 5 runs using 2500 primaries each). 
If you run on a cluster remember to provide the additional mgdraw routine.


4. Postprocessing
=================

Use the HITXML R package to use the resulting phase space files, and convert them into DDD / SPC (not yet implemented - TODO: SPC writer, dEdx output).
