# ACE-FTS-Handler

This program written in the Interactive Data Language (IDL). It defines
a class designed to streamline the workflow in handling data from the
Atmospheric Chemistry Experiment Fourier Transform Spectrometer (ACE-FTS).
In particular, this program aims to 
1) load vmr data from files and properly preprocess it for comparison between species.
2) Handle the removal of flagged data.
3) Handle the rebinning and accumulation of data based on user-inputted coordinates
4) Provide a set of standard plotting routines for the user.

As of v0.2.0, all of these features have been implemented.

v0.2.1 adds support for v4.1/4.2 ACE-FTS data.

TODO

Fix issues regarding calculation of time intervals, especially averaging
    a block of months over two or more years.
Complete the standard plotting 
Add support for parallel computations on multiprocessor systems. 
