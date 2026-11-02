Using the eigenprofiles from the linear stability analysis for A(z)-profiles 
========================

**Step 1:** Convert *.mat files into *.csv files 
The profiles from the linear stability analysis are saved using a *.mat format. For further analysis we prefer to use *.csv

```bash
python3 convert_mat_to_csv.py
```
Here, we transform profiles.mat into profiles.csv. These contains all cases from Figure 5.4. 

**Step 2:** Separate the combined data file into single files 

For better analysis, we split the combined data file into four single data files. As we have only four cases here, this is done manually. 

**Step 3:** Get grid resolution from numerical simulation 

The linear stability analysis uses a different vertical discretization than the numerical simulations. To match the profile with the discretization, we need information about the distribution of the vertical grid. 

For this extract the coordinates from the vertical direction from ParaView. In this extraction, some coordinates will occur multiple times. We remove any duplicates from the list. 

```bash
python3 filter_dublicates.py
```

**Step 4:** Interpolate eigenprofile to match vertical grid of the numerical simulations 
Specify within the script which case you want to evaluate. 

```bash
python3 plot_eigenprofile.py
```
As an output, you will get an interpolated profile of the the eigenprofile from the linear stability analysis. Depending on the case, the number of polynoms might have to be adjusted. 

**Step 5:** Transpose coefficients into an array 
In order to simply implement the coefficients into the numerical model, we can transpose the coefficients from the polynominal equation derived in Step 4 into an array. 

```bash
python3 transpose_coefficients.py
```