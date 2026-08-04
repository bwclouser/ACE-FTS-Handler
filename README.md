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

v0.3.0 harmonizes an arbitrary number of molecules rather than two. See the
DATA LAYOUT section of the class header in `acedata__define.pro` for how the
molecules and ratios are indexed.

## Configuration

The DMP version aliases passed to `setParam,'dmp_path',...` resolve against the
`ACE_DATA_ROOT` environment variable:

    setenv,'ACE_DATA_ROOT=/media/benjamin/Data/ACE'
    myData.setParam,'dmp_path','v5'

If `ACE_DATA_ROOT` is unset the previous hardcoded root is used, so existing
scripts keep working. A value that is not a recognised alias is treated as a
literal path.

## Tests

    idl -IDL_PATH "+<repo>:<IDL_DEFAULT>" -e 'ace_run_tests'
    idl -IDL_PATH "+<repo>:<IDL_DEFAULT>" -e 'ace_run_tests,/batch'   # exit status
    idl -IDL_PATH "+<repo>:<IDL_DEFAULT>" -e 'ace_benchmark'          # time calcBlocks

The tests inject synthetic data through the `ACETest` subclass, so no ACE-FTS
files are needed. Expected values are computed from the generating formulae
rather than recorded from a previous run.

Note that the class structure changes between versions, so a fresh IDL session
or `.full_reset` is needed after pulling rather than just a recompile.

TODO

Fix the calculation of time intervals. Two related defects: the per-bin
    year/month/day/hour limits never reach the occultation selection, so
    nTimes > 1 produces that many identical copies; and params.jullims is
    never populated, so the solidBlock presets select nothing.
Complete the standard plotting. makeMeanSDTP, makeHist, makeZonalPlot,
    setAxes and errorByExtrema currently MESSAGE that they are unimplemented.
Decide whether bin edges should be half-open. They are inclusive at both
    ends, so an occultation exactly on an interior edge falls in both
    neighbouring bins.
