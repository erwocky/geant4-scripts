# geant4-scripts
Post-processing scripts for Athena WFI Geant4 simulations.

Includes:
* Script to convert MIT Geant4 text step files to uniform FITS format pixel
  energy deposittion lists (`convert_to_fits_mit.py`)
* Script to sort primaries into frames and find events and particle tracks
  (blobs) from those frames (`find_events.py`)
* Wrapper script to parallelize either of the above (`wrapit.pl`)
* Original Perl script that did everything for OU data (`doit.pl`)
* Converted Python version of doit.pl (`doit.py`)
* Jupyter Notebook with doit.py code + some of the process (`dot.ipynb`)
* Perl initial random values for use in Python comparison to Perl (`framedist_vals.txt`)
* Tests to compare the two (`testing_txt` folder)

Notes on usage:
* The pre-processing script only works on the Geant4 output format
  specified, in this case MIT. Pre-processing scripts for OU format
  (comma-separated text files) and MPE format (ROOT files) are under
  construction.
* The pre-processing script is designed to process a single file called as
  an argument, e.g.:
  ```
  /path/to/geant4-scripts/convert_to_fits_mit.py input_50Mprotons/20201029073013_32358689_StepLog_0097.gdat
  ```
  and it will write a single FITS pixel files for that file. The
  multi-thread wrapper takes as arguments the script to run and a list of
  files to process. All must be specified by the path from the current
  working directory, e.g.:
  ```
  /path/to/geant4-scripts/wrapit.pl /path/to/geant4-scripts/convert_to_fits_mit.py input_50Mprotons/*StepLog*gdat
  ```
  The wrapper will fork processes to run the wrapped script in parallel; by
  default it will run 20 parallel threads. This can be changed by editing
  the `wrapit.pl` script $MAXPROCS setting. 
* The post-processing script acts and can be wrapped in the same way. It
  takes as argument a raw pixel FITS file and outputs a single frame,
  pixel, event, and blob file. Right now those files must be concatenated
  by hand; this will be automated in a future release.
