# A Python package for sccomp - Tests differences in cell type proportions and variability from single-cell data
Cellular omics such as single-cell genomics, proteomics, and microbiomics allow the characterization of tissue and microbial community composition, which can be compared between conditions to identify biological drivers. This strategy has been critical to unveiling markers of disease progression in conditions such as cancer and pathogen infections.

For cellular omic data, no method for differential variability analysis exists, and methods for differential composition analysis only take a few fundamental data properties into account. Here we introduce sccomp, a generalised method for differential composition and variability analyses capable of jointly modelling data count distribution, compositionality, group-specific variability, and proportion mean-variability association, while being robust to outliers.

sccomp is an extensive analysis framework that allows realistic data simulation and cross-study knowledge transfer. We demonstrate that mean-variability association is ubiquitous across technologies, highlighting the inadequacy of the very popular Dirichlet-multinomial modeling and providing essential principles for differential variability analysis.

## installation
Here is a demo of installation in editable mode


```python
%pip install -e .
```

    Defaulting to user installation because normal site-packages is not writeable
    Obtaining file:///home/chzhan1/Python/SAiGENCI/sccompPy
    Requirement already satisfied: pandas>=1.0.0 in /home/chzhan1/.local/lib/python3.9/site-packages (from sccompPy==0.1.0) (2.2.2)
    Requirement already satisfied: numpy>=1.18.0 in /home/chzhan1/.local/lib/python3.9/site-packages (from sccompPy==0.1.0) (1.26.4)
    Requirement already satisfied: patsy>=0.5.0 in /home/chzhan1/.local/lib/python3.9/site-packages (from sccompPy==0.1.0) (0.5.6)
    Requirement already satisfied: cmdstanpy>=1.0.0 in /home/chzhan1/.local/lib/python3.9/site-packages (from sccompPy==0.1.0) (1.2.4)
    Requirement already satisfied: stanio<2.0.0,>=0.4.0 in /home/chzhan1/.local/lib/python3.9/site-packages (from cmdstanpy>=1.0.0->sccompPy==0.1.0) (0.5.1)
    Requirement already satisfied: tqdm in /home/chzhan1/.local/lib/python3.9/site-packages (from cmdstanpy>=1.0.0->sccompPy==0.1.0) (4.66.2)
    Requirement already satisfied: tzdata>=2022.7 in /home/chzhan1/.local/lib/python3.9/site-packages (from pandas>=1.0.0->sccompPy==0.1.0) (2024.1)
    Requirement already satisfied: python-dateutil>=2.8.2 in /home/chzhan1/.local/lib/python3.9/site-packages (from pandas>=1.0.0->sccompPy==0.1.0) (2.9.0.post0)
    Requirement already satisfied: pytz>=2020.1 in /home/chzhan1/.local/lib/python3.9/site-packages (from pandas>=1.0.0->sccompPy==0.1.0) (2024.1)
    Requirement already satisfied: six in /home/chzhan1/.local/lib/python3.9/site-packages (from patsy>=0.5.0->sccompPy==0.1.0) (1.16.0)
    Installing collected packages: sccompPy
      Attempting uninstall: sccompPy
        Found existing installation: sccompPy 0.1.0
        Uninstalling sccompPy-0.1.0:
          Successfully uninstalled sccompPy-0.1.0
    [33m  WARNING: Value for scheme.platlib does not match. Please report this to <https://github.com/pypa/pip/issues/10151>
      distutils: /home/chzhan1/.local/lib/python3.9/site-packages
      sysconfig: /home/chzhan1/.local/lib64/python3.9/site-packages[0m
    [33m  WARNING: Additional context:
      user = True
      home = None
      root = None
      prefix = None[0m
      Running setup.py develop for sccompPy
    Successfully installed sccompPy-0.1.0
    Note: you may need to restart the kernel to use updated packages.


## Import `sccompy` package 


```python
import sccompPy
```

## Load embeded dataset


```python
import pkg_resources
import pandas as pd

data_file_path = pkg_resources.resource_filename("sccompPy", "data/count_obj.csv")

count_obj = pd.read_csv(data_file_path)
count_obj
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sample</th>
      <th>type</th>
      <th>phenotype</th>
      <th>count</th>
      <th>cell_group</th>
      <th>proportion</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>10x_6K</td>
      <td>benign</td>
      <td>b_cell_macrophage_precursor_or_follicular_LTB_...</td>
      <td>42</td>
      <td>BM</td>
      <td>0.008350</td>
    </tr>
    <tr>
      <th>1</th>
      <td>10x_6K</td>
      <td>benign</td>
      <td>B_cell:immature</td>
      <td>361</td>
      <td>B1</td>
      <td>0.071769</td>
    </tr>
    <tr>
      <th>2</th>
      <td>10x_6K</td>
      <td>benign</td>
      <td>B_cell:immature_IGLC3_IGLC2</td>
      <td>57</td>
      <td>B2</td>
      <td>0.011332</td>
    </tr>
    <tr>
      <th>3</th>
      <td>10x_6K</td>
      <td>benign</td>
      <td>B_cell:Memory_ITM2C_IGHA1_MZB1_JCHAIN</td>
      <td>40</td>
      <td>B3</td>
      <td>0.007952</td>
    </tr>
    <tr>
      <th>4</th>
      <td>10x_6K</td>
      <td>benign</td>
      <td>Dendritic_CD11_CD1_high_mito</td>
      <td>75</td>
      <td>Dm</td>
      <td>0.014911</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>715</th>
      <td>SRR7244582</td>
      <td>benign</td>
      <td>T_cell:CD8+_GZMK_DUSP2_LYAR_CCL5</td>
      <td>197</td>
      <td>CD8 2</td>
      <td>0.060727</td>
    </tr>
    <tr>
      <th>716</th>
      <td>SRR7244582</td>
      <td>benign</td>
      <td>T_cell:CD8+_non_activated</td>
      <td>320</td>
      <td>CD8 3</td>
      <td>0.098644</td>
    </tr>
    <tr>
      <th>717</th>
      <td>SRR7244582</td>
      <td>benign</td>
      <td>T_cell:CD8+_PPBP_SAT1</td>
      <td>39</td>
      <td>CD8 4</td>
      <td>0.012022</td>
    </tr>
    <tr>
      <th>718</th>
      <td>SRR7244582</td>
      <td>benign</td>
      <td>T_cell:CD8+_S100B</td>
      <td>88</td>
      <td>CD8 5</td>
      <td>0.027127</td>
    </tr>
    <tr>
      <th>719</th>
      <td>SRR7244582</td>
      <td>benign</td>
      <td>T_cell:CD8low_TIMP1_PPBP</td>
      <td>107</td>
      <td>CD8 6</td>
      <td>0.032984</td>
    </tr>
  </tbody>
</table>
<p>720 rows × 6 columns</p>
</div>



## `sccomp_estimate` function


```python
estimate_res = sccompPy.sccomp_estimate(
    data = count_obj,
    formula_composition = '~ 0 + type', 
    sample = 'sample',
    cell_group = 'cell_group',
    count = 'count',
    verbose = False
)
```

    16:22:09 - cmdstanpy - INFO - CmdStan start processing



    chain 1 |          | 00:00 Status



    chain 2 |          | 00:00 Status



    chain 3 |          | 00:00 Status



    chain 4 |          | 00:00 Status


                                                                                                                                                                                                                                                                                                                                    

    16:22:18 - cmdstanpy - INFO - CmdStan done processing.
    16:22:18 - cmdstanpy - WARNING - Non-fatal error during sampling:
    Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is inf, but must be positive finite! (in 'glm_multi_beta_binomial.stan', line 214, column 16 to line 219, column 19) (in 'glm_multi_beta_binomial.stan', line 653, column 3 to line 683, column 8)
    Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is inf, but must be positive finite! (in 'glm_multi_beta_binomial.stan', line 214, column 16 to line 219, column 19) (in 'glm_multi_beta_binomial.stan', line 653, column 3 to line 683, column 8)
    Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is inf, but must be positive finite! (in 'glm_multi_beta_binomial.stan', line 214, column 16 to line 219, column 19) (in 'glm_multi_beta_binomial.stan', line 653, column 3 to line 683, column 8)
    Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in 'glm_multi_beta_binomial.stan', line 753, column 2 to column 53)
    	Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in 'glm_multi_beta_binomial.stan', line 753, column 2 to column 53)
    	Exception: Exception: beta_binomial_lpmf: First prior sample size parameter[1] is inf, but must be positive finite! (in 'glm_multi_beta_binomial.stan', line 214, column 16 to line 219, column 19) (in 'glm_multi_beta_binomial.stan', line 653, column 3 to line 683, column 8)
    Consider re-running with show_console=True if the above output is unclear!


    



```python
estimate_res.keys()
```




    dict_keys(['fit', 'model_input', 'truncation_df2', 'sample', 'cell_group', 'count', 'formula_composition', 'formula_variability', 'noise_model'])



## `sccomp_test` function


```python
test_res = sccompPy.sccomp_test(estimate_res, contrasts= ['type[cancer] - type[benign]'])
```

    Warning: These elements require backquotes: ['type[cancer]', 'type[benign]'], sccompPy will auto quote them.
    Warning: These elements require backquotes: ['type[cancer]', 'type[benign]'], sccompPy will auto quote them.


    /home/chzhan1/Python/SAiGENCI/sccompPy/sccompPy/utilities.py:295: FutureWarning: The default of observed=False is deprecated and will be changed to True in a future version of pandas. Pass observed=False to retain current behavior or observed=True to adopt the future default and silence this warning.
      grouped = draws.groupby([cell_group, 'M', 'parameter'])
    /home/chzhan1/Python/SAiGENCI/sccompPy/sccompPy/utilities.py:312: FutureWarning: The default of observed=False is deprecated and will be changed to True in a future version of pandas. Pass observed=False to retain current behavior or observed=True to adopt the future default and silence this warning.
      summary["FDR"] = summary.groupby("parameter")[f"{prefix}pH0"].transform(lambda pH0: get_FDR(pH0))


### `sccomp_test` returns a `dict` where the first element - *result* contains the result table


```python
test_res
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>cell_group</th>
      <th>parameter</th>
      <th>c_lower</th>
      <th>c_effect</th>
      <th>c_upper</th>
      <th>c_pH0</th>
      <th>FDR</th>
      <th>N_Eff</th>
      <th>R_k_hat</th>
      <th>v_lower</th>
      <th>v_effect</th>
      <th>v_upper</th>
      <th>v_pH0</th>
      <th>factor</th>
      <th>design_matrix_col</th>
      <th>count_data</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>B1</td>
      <td>type[benign]</td>
      <td>0.900753</td>
      <td>1.155635</td>
      <td>1.386001</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>4821.08</td>
      <td>0.999812</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>type</td>
      <td>type[benign]</td>
      <td>sample    type        phenotype  co...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>B1</td>
      <td>type[cancer]</td>
      <td>0.150861</td>
      <td>0.494823</td>
      <td>0.798735</td>
      <td>0.01475</td>
      <td>0.001500</td>
      <td>5975.15</td>
      <td>0.999897</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>type</td>
      <td>type[cancer]</td>
      <td>sample    type        phenotype  co...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>B1</td>
      <td>type[cancer] - type[benign]</td>
      <td>-1.065256</td>
      <td>-0.663332</td>
      <td>-0.288341</td>
      <td>0.00325</td>
      <td>0.001000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>sample    type        phenotype  co...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>B2</td>
      <td>type[benign]</td>
      <td>0.412158</td>
      <td>0.724773</td>
      <td>1.017911</td>
      <td>0.00050</td>
      <td>0.000036</td>
      <td>7054.69</td>
      <td>0.999426</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>type</td>
      <td>type[benign]</td>
      <td>sample    type                    p...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>B2</td>
      <td>type[cancer]</td>
      <td>-0.411091</td>
      <td>0.028000</td>
      <td>0.418074</td>
      <td>0.63900</td>
      <td>0.087287</td>
      <td>7181.51</td>
      <td>0.999822</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>type</td>
      <td>type[cancer]</td>
      <td>sample    type                    p...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>103</th>
      <td>TM2</td>
      <td>type[cancer]</td>
      <td>-1.214596</td>
      <td>-0.892860</td>
      <td>-0.598306</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>5793.83</td>
      <td>0.999928</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>type</td>
      <td>type[cancer]</td>
      <td>sample    type                     ...</td>
    </tr>
    <tr>
      <th>104</th>
      <td>TM2</td>
      <td>type[cancer] - type[benign]</td>
      <td>-0.167554</td>
      <td>0.276814</td>
      <td>0.717797</td>
      <td>0.20400</td>
      <td>0.053104</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>sample    type                     ...</td>
    </tr>
    <tr>
      <th>105</th>
      <td>TM3</td>
      <td>type[benign]</td>
      <td>-1.764759</td>
      <td>-0.827435</td>
      <td>0.264891</td>
      <td>0.09350</td>
      <td>0.015492</td>
      <td>1817.40</td>
      <td>1.000410</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>type</td>
      <td>type[benign]</td>
      <td>sample    type                     ...</td>
    </tr>
    <tr>
      <th>106</th>
      <td>TM3</td>
      <td>type[cancer]</td>
      <td>-3.906174</td>
      <td>-2.667245</td>
      <td>-1.388576</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>1453.44</td>
      <td>0.999974</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>type</td>
      <td>type[cancer]</td>
      <td>sample    type                     ...</td>
    </tr>
    <tr>
      <th>107</th>
      <td>TM3</td>
      <td>type[cancer] - type[benign]</td>
      <td>-3.109568</td>
      <td>-1.838550</td>
      <td>-0.743080</td>
      <td>0.00050</td>
      <td>0.000200</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>sample    type                     ...</td>
    </tr>
  </tbody>
</table>
<p>108 rows × 16 columns</p>
</div>


