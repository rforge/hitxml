=========================================================
sg74001 - FLUKA/FLAIR files and R script to get DEDX data

=========================================================

6 Fluka runs looping Z from 3 to 8 with only 1 primary each,
producing and dEdx output table for water (customize,
make sure the same setting as in the SPC/DDD etc. apply!).

This is because the particle HEAVYION is determined by
the settings for the primary, i.e. it is EITHER Li, or Be, or...
etc.

Proton and Helium have in contrast their on ids.
