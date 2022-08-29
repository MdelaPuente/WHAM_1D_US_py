## Tests

You can test the `wham1D.py` and `errors_wham1D.py` on the CV files contained in the `COLVAR` folder with the `input.txt` file as input. In the `example_outputs` folder there are the outputs obtained by running (from the `tests` folder with Python 3.8.8):
```sh
python ../src/wham1D.py input.txt TEST
python ../src/errors_wham1D.py input.txt TEST
```

**Important:** the CV files provided are "real" CV files from a deprotonation US calculation (CV=coordination number), but they correspond to only the first 5 windows (out of 50) and have been reduced to the first 10000 frames (out of 280000), so the results are unphysical. They should only be used to check the scripts and understand the different file formats.
