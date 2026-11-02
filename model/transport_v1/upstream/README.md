Summary
=========================

This repository contains the source code to reproduce the numerical simulations presented in:

Carina Bringedal, Stefanie Kiemle, C. J. van Duijn and Rainer Helmig. Impact of saturation on evaporation-driven density instabilities1
in porous media: mathematical and numerical analysis. Submitted 2024.

The results of the numerical simulations used in the paper and produced by the code presented in this git-repository are available in the following Dataset on [DaRUS](https://darus.uni-stuttgart.de/): 
Kiemle, Stefanie, 2024, "DuMuX code for modelling evaporation-driven density instabilities in unsaturated porous media", https://doi.org/10.18419/darus-4610 

Installation
=========================
To install the module with its dependencies run:
```bash
python3 install_kiemle2023b.py
```

Compile and run the simulation
=========================
To compile the simulation use the following commands:
```bash
cd build-cmake/appl/simpleRichards
make test_model_salt
```
To then run the simulation use the following command:
```bash
./test_model_salt
```

To run all simulation presented in the paper, run the following bash script in the folder appl:

```bash
cd kiemle2023b/appl/simpleRichards/
chmod u+x simulation_analysis.sh
./simulation_analysis.sh
```
Evaluation
=========================

**Step 1:** Extract data of the first grid row of all .vtu file. For this we use the script extractlinedata.py. As a large amount of files are called, we run a script called run_extractlinedata.py where extractlinedata.py is used as a subprocess.

```bash
mkdir Results/PlotOverLineRichards/
cd bin/
python3 run_extractlinedata.py
```

**Step 2:** Determine mean concentration over time, the numerical onset and the dominant number of waves in the first row. The just extracted data is now processed with the script wavenumber.py.


```bash
mkdir Results/Wavelength_Richards/
mkdir Results/Signal_Plots_Richards/
cd bin/
python3 wavelength.py
```

**Step 3:** Create plots. The following scripts are used to create the Figures from the corresponding publications. Additionally, two further scripts have been added which have been used to analyse the system, but the results have not been used in the publication

```bash
cd bin/
```

To recreate Figure 4.1, run the following script:

```bash
python3 plot_perturbation_random.py
```
To recreate Figure 4.2, run the following script:
```bash
python3 plot_perturbation_periodic_WL025.py
```

To recreate Figure 5.4, run the following script:
```bash
python3 plot_periodic_eigenprofile_waterlevel.py
```
Additionally, there are two scripts available to create plots to compare different bottom pressures. The first script is used if the random perturabtion is used:

```bash
python3 plot_random_waterlevel.py
```
and the second one if the periodic perturbation only in the top row is used:

```bash
python3 plot_periodic_waterlevel.py
```




Preparing the Sources
=========================

Additional to the software mentioned in README you'll need the
following programs installed on your system:

  cmake >= 3.13

Getting started
---------------

If these preliminaries are met, you should run

  dunecontrol all

which will find all installed dune modules as well as all dune modules
(not installed) which sources reside in a subdirectory of the current
directory. Note that if dune is not installed properly you will either
have to add the directory where the dunecontrol script resides (probably
./dune-common/bin) to your path or specify the relative path of the script.

Most probably you'll have to provide additional information to dunecontrol
(e. g. compilers, configure options) and/or make options.

The most convenient way is to use options files in this case. The files
define four variables:

CMAKE_FLAGS      flags passed to cmake (during configure)

An example options file might look like this:

#use this options to configure and make if no other options are given
CMAKE_FLAGS=" \
-DCMAKE_CXX_COMPILER=g++-5 \
-DCMAKE_CXX_FLAGS='-Wall -pedantic' \
-DCMAKE_INSTALL_PREFIX=/install/path" #Force g++-5 and set compiler flags

If you save this information into example.opts you can pass the opts file to
dunecontrol via the --opts option, e. g.

  dunecontrol --opts=example.opts all

More info
---------

See

     dunecontrol --help

for further options.


The full build system is described in the dune-common/doc/buildsystem (Git version) or under share/doc/dune-common/buildsystem if you installed DUNE!
