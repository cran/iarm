# News for package 'iarm'

## Changes in Version 0.4.1

* Unadjusted and adjusted p values included in several functions.
* Rest score calculation for LD detection in both possible ways. 
* fixed a bug in item-restscore().
* fixed a bug in ICCplot().

## Changes in Version 0.4.0

* added a plotting function for ICCs.

## Changes in Version 0.3.0

* Correction method for multiple testing included.
* person_estimates() accepts raschmodel objects.
* clr_tests() accepts single factors or a data frame with several variables.
* item_obsexp(), item_restscore(), out_infit() and boot_fit() use only complete data.  
* Package mRm is not maintained anymore, mrm() has been replaced by raschmodel().
* fixed a bug in item_obsexp().


## Changes in Version 0.2.0

* bug fix for input arguments of partgam
* added functions partgam_LD and partgam_DIF
* person_estimates accepts now items with different number of levels, calculates person measures for all persons in the data set.
* test_prop accepts pcmodel object as input.
