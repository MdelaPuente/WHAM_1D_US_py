<div id="top"></div>

<!-- PROJECT SHIELDS -->

[![GNU AGPL v3.0 License][license-shield]][license-url]

<!-- TABLE OF CONTENTS -->

<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
<div id="about"></div>

## About The Project

A python implementation of the Weighted Histogram Analysis Method (WHAM) to obtain free-energy profiles from Umbrella Sampling simulations (`wham1D.py`) as well as a script for computing statistically significant errors (`errors_wham1D.py`) as prescribed by Zhu and Hummer (see below). 

It works only in 1D and makes use of files containing the collective variables (CV) (one per window) formatted as the default [Plumed](https://www.plumed.org/) COLVAR files (one line per MD step). The `wham1D.py` code writes an output file (in text format) with the grid of CV values and associated free-energy values and saves `.png` images of the profile and other computed quantities (biased and unbiased probability distributions, window free-energy segments, convergence info, etc.). If these images are unwanted or if you do not want to import `matplotlib.pyplot` you can simply comment the `import` line at the beginning of the code as well as the `MAKE FIGURES` section at the end and you should still get the text file. The errors are calculated separately by `errors_wham1D.py` and printed to a separate text file. (see `tests`folder)

**Important:** you should define sensible minimal CV value (`qmin`), CV spacing (`Δq`), number of bins for the CV grid (`Nb`) and `kT` value (in eV, if different fromt that at 300 K) directly inside the codes in the `Define parameters and recover data` section. It might be included in the input file in a later version of the code (work in progress).

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- GETTING STARTED -->
<div id="getting-started"></div>

## Getting Started

<div id="prerequisites"></div>

### Prerequisites

This Python script has been tested with Python 3.6 and 3.7.3 and is expected to work with later versions as well.

The following packages are required for the codes to work:
* `os` 
* `numpy`
* `scipy.integrate`
* `matplotlib.pyplot` (this is optional and can be removed with the corresponding section if only a text file is wanted)

<div id="installation"></div>

### Installation

1. Clone the repository
   ```sh
   git clone https://github.com/laagegroup/WHAM_1D_US_py.git
   ```
   ```
2. Go into the tests folder

3. Test the program following the `tests/REDME.md` guidelines. 

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
<div id="usage"></div>

## Usage

Both the `wham1D.py` and `errors_wham1D.py` codes take the same text file as input (see `tests/README.md`) and a "suffix" string for the naming of output files. You should be able to call both scripts right after download: 

   ```sh
   python wham1D.py input_file.txt SUFFIX
   ```

The input_file.txt is a file containing a line for each window simulation formatted as:

   ```sh
   # Possible comment line
   PATH_TO_CV_FILE1/CV_FILE_NAME1 TARGET_CV_VALUE1 KAPPA_VALUE_IN_EV1 
   PATH_TO_CV_FILE2/CV_FILE_NAME2 TARGET_CV_VALUE2 KAPPA_VALUE_IN_EV2 
   ```

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- LICENSE -->
<div id="license"></div>

## License

Distributed under the GNU Affero General Public License v3.0. See `LICENSE` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->
<div id="acknowledgments"></div>

## Acknowledgments & Sources

* Roux, B. The Calculation of the Potential of Mean Force Using Computer Simulations. Computer Physics Communications 1995, 91 (1–3), 275–282.
* Zhu, F.; Hummer, G. Convergence and Error Estimation in Free Energy Calculations Using the Weighted Histogram Analysis Method. J. Comput. Chem. 2012, 33 (4), 453–465.
* Tuckerman, M.E. Statistical Mechanics: Theory and Molecular Simulation, Chapter 8, pp.340-344, Oxford University Press, 2010.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[license-shield]: https://img.shields.io/github/license/laagegroup/0_Template.svg?style=for-the-badge
[license-url]: https://github.com/laagegroup/0_Template/blob/main/LICENSE
