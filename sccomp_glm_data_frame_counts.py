import pandas as pd
import numpy as np
from patsy import dmatrices
import cmdstanpy

def sccomp_glm_data_frame_counts(
    _data,
    formula_composition = '~ 1',
    formula_variability = '~ 1',
    _sample = 'sample',
    _cell_group = 'cell_group',
    _count = 'count',

    # Secondary arguments, will fill more
    truncation_ajustment = 1.1,
    inference_method = "variational",
    bimodal_mean_variability_association = False,
    use_data = True,
    prior_mean = dict(intercept = [0,1], coefficients = [0,1]),                    
    prior_overdispersion_mean_association = dict(intercept = [5, 2], slope = [0,  0.6], standard_deviation = [20, 40]),
    exclude_priors = False
):
    
    # fill complete formula to include _count as a variable prefix.
    formula_composition = f"{_count} {formula_composition}"
    formula_variability = f"{_count} {formula_variability}"

    # Initialize dictionary to hold data for Stan model
    data_for_model = {}

    # N: Total number of unique samples
    data_for_model['N'] = _data[_sample].nunique()

    # M: Total number of unique cell groups
    data_for_model['M'] = _data[_cell_group].nunique()

    # exposure: Calculate exposure as the sum of counts for each sample
    data_for_model['exposure'] = np.array(_data.groupby(_sample)[_count].sum())

    # is_proportion: Boolean flag to indicate if data represents proportions
    data_for_model['is_proportion'] = _data[_count].max() <= 1

    # y and y_proportion 
    ## Pivot table to arrange counts by sample and cell group
    y = _data.pivot_table(
        index=_sample,          # Rows indexed by 'sample'
        columns=_cell_group,    # Columns are 'cell_group' values
        values=_count,          # Values in the table are from 'count' column
        fill_value=0             # Fill NaN values with 0
    )

    if data_for_model['is_proportion']:
        y_proportion = y
        y = pd.DataFrame(columns=y.columns)
    else:
        y_proportion = pd.DataFrame(columns=y.columns)
        y = y.astype(int)
    
    data_for_model['y'] = y
    data_for_model['y_proportion'] = y_proportion

    # X : Generate design matrices for composition and variability formulas
    y_composition, X_composition = dmatrices(formula_composition, _data)
    tmp_data_composition = pd.DataFrame(X_composition, columns=X_composition.design_info.column_names)
    tmp_data_composition[_sample] = _data[_sample]
    tmp_data_composition.drop_duplicates(inplace=True)
    tmp_data_composition.set_index(_sample, inplace=True)

    data_for_model['X'] = tmp_data_composition

    # Xa and XA
    y_variability, X_variability = dmatrices(formula_variability, _data)
    tmp_data_variability = pd.DataFrame(X_variability, columns=X_variability.design_info.column_names)
    tmp_data_variability[_sample] = _data[_sample]
    tmp_data_variability.drop_duplicates(inplace=True)
    tmp_data_variability.set_index(_sample, inplace=True)

    data_for_model['Xa'] = tmp_data_variability
    data_for_model['XA'] = tmp_data_variability.drop_duplicates().reset_index(drop=True)

    data_for_model['C'] = data_for_model['X'].shape[1]
    data_for_model['A'] = data_for_model['XA'].shape[1]
    data_for_model['Ar'] = data_for_model['XA'].shape[0]
    data_for_model['truncation_ajustment'] = truncation_ajustment
    data_for_model['is_vb'] = int(inference_method in ["variational", "pathfinder"])
    data_for_model['bimodal_mean_variability_association'] = bimodal_mean_variability_association
    data_for_model['use_data'] = use_data
    data_for_model['is_random_effect'] = 0

    # random effect coef --- NOT POLISHED
    data_for_model["ncol_X_random_eff"] = [0,0]
    data_for_model['n_random_eff'] = 0
    data_for_model["n_groups"] = [0,0]
    data_for_model["X_random_effect"] = pd.DataFrame(index=_data[_sample].unique())
    data_for_model["X_random_effect_2"] = pd.DataFrame(index=_data[_sample].unique())
    data_for_model['group_factor_indexes_for_covariance'] = pd.DataFrame()
    data_for_model['group_factor_indexes_for_covariance_2'] = pd.DataFrame()
    data_for_model['how_many_factors_in_random_design'] = [0,0]

    # subject to change
    data_for_model['grainsize'] = 1
    data_for_model['enable_loo'] = False
    data_for_model['is_truncated'] = 0

    data_for_model['truncation_up'] = pd.DataFrame(-1, index=_data[_sample].unique(), columns=_data[_cell_group].unique())
    data_for_model['truncation_down'] = pd.DataFrame(-1, index=_data[_sample].unique(), columns=_data[_cell_group].unique())

    # check this later
    data_for_model['truncation_not_idx'] = list(range(1, _data.shape[0]+1))

    data_for_model['TNS'] =  _data.shape[0]
    data_for_model['truncation_not_idx_minimal'] = np.empty((0, 2))
    data_for_model['TNIM'] = 0

    # intercept
    data_for_model['intercept_in_design'] = data_for_model['X'].iloc[:,0].unique().tolist() == [1]

    if data_for_model['intercept_in_design'] or X_variability.design_info.term_names == [] or X_variability.design_info.term_names == ['Intercept']:
        data_for_model['A_intercept_columns'] = 1
    else:
        data_for_model['A_intercept_columns'] = _data[[col for col in X_composition.design_info.term_names if col in _data.columns]].drop_duplicates().shape[0]

    if data_for_model['intercept_in_design']:
        data_for_model['B_intercept_columns'] = 1
    else:
        data_for_model['B_intercept_columns'] = _data[[col for col in X_composition.design_info.term_names if col in _data.columns]].drop_duplicates().shape[0]
    
    # check this later
    data_for_model['user_forced_truncation_not_idx'] = list(range(1, _data.shape[0]+1))

    # fill params from args
    data_for_model['prior_prec_intercept'] = prior_overdispersion_mean_association['intercept']
    data_for_model['prior_prec_slope'] = prior_overdispersion_mean_association['slope']
    data_for_model['prior_prec_sd'] = prior_overdispersion_mean_association['standard_deviation']

    data_for_model['prior_mean_intercept'] = prior_mean['intercept']
    data_for_model['prior_mean_coefficients'] = prior_mean['coefficients']

    data_for_model['exclude_priors'] = exclude_priors

    # Compile Stan model
    stan_model = cmdstanpy.CmdStanModel(stan_file='./stan/glm_multi_beta_binomial.stan')

    # Run Stan model sampling
    fit = stan_model.sample(
        data=data_for_model,
        show_console=True
    )

    return fit, data_for_model # keep outputing data_for_model for debugging at dev stage