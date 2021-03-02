# geant4-scripts
Post-processing scripts for Athena WFI Geant4 simulations.

Includes:
* Script to convert MIT Geant4 text step files to uniform FITS format pixel
  energy deposittion lists (`convert_to_fits_mit.py`)
* Script to sort primaries into frames and find events and particle tracks
  (blobs) from those frames (`find_events.py`)
* Wrapper script to parallelize either of the above (`wrapit.pl`)
* Concatenation script to combine output from above into single files (`combit.py`).
* Test data generated from MIT Geant4 simulations (`test_data/`)

* Development directory containing deprecated scripts and testing (`devel/`)
* Original Perl script that did everything for OU data (`devel/doit.pl`)
* Converted Python version of doit.pl (`devel/doit.py`)
* Jupyter Notebook with doit.py code + some of the process (`devel/dot.ipynb`)
* Perl initial random values for use in Python comparison to Perl (`devel/framedist_vals.txt`)
* Tests to compare the two (`devel/testing_txt` folder)

Notes on usage:
* The pre-processing scripts only work on the Geant4 output format
  specified, in this case MIT and OU. A pre-processing script for 
  MPE format (ROOT files) is under construction.
* The pre-processing script is designed to process a list of files called as
  arguments, e.g.:
  ```
  ./convert_to_fits_mit.py test_data/20201006173157_StepLog_0000.gdat.gz
  ```
  and it will write a single FITS pixel file for each file. The
  multi-thread wrapper takes as arguments the script to run and the list of
  files to process. All must be specified by the path from the current
  working directory, e.g.:
  ```
  ./wrapit.pl ./convert_to_fits_mit.py test_data/*StepLog*.gdat.gz
  ```
  The wrapper will fork processes to run the wrapped script in parallel; by
  default it will run 20 parallel threads. This can be changed by editing
  the `wrapit.pl` script $MAXPROCS setting. 
* The post-processing script acts and can be wrapped in the same way. It
  takes as arguments a list of raw pixel FITS files output by the
  pre-processing script:
  ```
  ./wrapit.pl ./find_events.py test_data/rawpix_*.fits
  ```
  and outputs a single frame, pixel, event, and blob file for each input file.
  Those four sets of files can be concatenated into single large files, with running
  columns like FRAME and TIME and ID numbers incremented and header information updated:
  ```
  ./combit.py test_data/
  ```
  Note the different syntax; this takes only the directory path which contains the FITS
  files, and will combine whatever it finds there with filenames matching what's expected
  from `find_events.py`.
