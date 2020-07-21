## filetObs.py
Returns a single vector of summary stats from a window-based file as either the
mean or median value. Mask option to remove statistics affected by missing data

## filet_clustering.py
Splits and concatenates result from produced by FILET into direction and status

## filet_clustering_mean.py
similar to above, but uses a rolling mean to merge windows instead of a strict cutoff

## filet_pipeline.sh
Hard coded example of running FILET

## make_filtet_mask.py
Creates masking file for use with FILET to mask training simulations

## make_trainvecs_filetMP.py
python wrapper to calculate the training vector of summary statistics used with FILET using multi-processing

## msMaskAllRows 
https://github.com/kr-colab/FILET
FILET exe to apply mask to ms-type file

## normalizeTwoPopnStats.py
https://github.com/kr-colab/FILET
Normalize the FILET stats vector by number of sites

## removeNedOutColumnsFromMsFileKeepStars.py
https://github.com/kr-colab/FILET
Remove masked columns from ms-type file, but retain `*` with msmove

## twoPopnStats_forML
https://github.com/kr-colab/FILET
FILET exe to calculate stats vector
