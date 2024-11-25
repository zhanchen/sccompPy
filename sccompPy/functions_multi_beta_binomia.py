import pandas as pd
import numpy as np
from patsy import dmatrices
import cmdstanpy
import os
import pkg_resources

def sccomp_glm_data_frame_counts(
    data,
    formula_composition = '~ 1',
    formula_variability = '~ 1',
    sample = 'sample',
    cell_group = 'cell_group',
    count = 'count',

    # Secondary arguments, will fill more
    contrasts=None,
    prior_mean = dict(intercept = [0,1], coefficients = [0,1]),     
    prior_overdispersion_mean_association = dict(intercept = [5, 2], slope = [0,  0.6], standard_deviation = [20, 40]),
    percent_false_positive=5,
    # check_outliers = True, # not completed yet
    variational_inference = None,
    inference_method = "variational",
    truncation_ajustment = 1.1,
    test_composition_above_logit_fold_change = 0.1,
    sample_cell_group_pairs_to_exclude = None,
    output_directory = "sccomp_draws_files",
    verbose = False,
    exclude_priors = False,
    bimodal_mean_variability_association = False,
    enable_loo = False,  
    use_data = True,
    cores = 4,
    mcmc_seed = None,
    max_sampling_iterations=20000,
    pass_fit=True,
    **kwargs
):
    
    # fill complete formula to include count as a variable prefix.
    formula_composition = f"{count} {formula_composition}"
    formula_variability = f"{count} {formula_variability}"

    # Initialize dictionary to hold data for Stan model
    data_for_model = {}

    # N: Total number of unique samples
    data_for_model['N'] = data[sample].nunique()

    # M: Total number of unique cell groups
    data_for_model['M'] = data[cell_group].nunique()

    # exposure: Calculate exposure as the sum of counts for each sample
    data_for_model['exposure'] = np.array(data.groupby(sample)[count].sum())

    # is_proportion: Boolean flag to indicate if data represents proportions
    data_for_model['is_proportion'] = data[count].max() <= 1

    # y and y_proportion 
    ## Pivot table to arrange counts by sample and cell group
    y = data.pivot_table(
        index=sample,          # Rows indexed by 'sample'
        columns=cell_group,    # Columns are 'cell_group' values
        values=count,          # Values in the table are from 'count' column
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
    y_composition, X_composition = dmatrices(formula_composition, data)
    tmp_data_composition = pd.DataFrame(X_composition, columns=X_composition.design_info.column_names)
    tmp_data_composition[sample] = data[sample]
    tmp_data_composition.drop_duplicates(inplace=True)
    tmp_data_composition.set_index(sample, inplace=True)

    data_for_model['X'] = tmp_data_composition

    # Xa and XA
    y_variability, X_variability = dmatrices(formula_variability, data)
    tmp_data_variability = pd.DataFrame(X_variability, columns=X_variability.design_info.column_names)
    tmp_data_variability[sample] = data[sample]
    tmp_data_variability.drop_duplicates(inplace=True)
    tmp_data_variability.set_index(sample, inplace=True)

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
    data_for_model["X_random_effect"] = pd.DataFrame(index=data[sample].unique())
    data_for_model["X_random_effect_2"] = pd.DataFrame(index=data[sample].unique())
    data_for_model['group_factor_indexes_for_covariance'] = pd.DataFrame()
    data_for_model['group_factor_indexes_for_covariance_2'] = pd.DataFrame()
    data_for_model['how_many_factors_in_random_design'] = [0,0]

    # subject to change
    data_for_model['grainsize'] = 1
    data_for_model['enable_loo'] = enable_loo
    data_for_model['is_truncated'] = 0

    data_for_model['truncation_up'] = pd.DataFrame(-1, index=data[sample].unique(), columns=data[cell_group].unique())
    data_for_model['truncation_down'] = pd.DataFrame(-1, index=data[sample].unique(), columns=data[cell_group].unique())

    # check this later
    data_for_model['truncation_not_idx'] = list(range(1, data.shape[0]+1))

    data_for_model['TNS'] =  data.shape[0]
    data_for_model['truncation_not_idx_minimal'] = np.empty((0, 2))
    data_for_model['TNIM'] = 0

    factor_parameter_dictionary = []
    for term, slice_obj in X_composition.design_info.term_slices.items():
        term_columns =  X_composition.design_info.column_names[slice_obj]
        factor_parameter_dictionary.extend([(col, term.factors[0].name()) for col in term_columns])

    factor_parameter_dictionary = pd.DataFrame(factor_parameter_dictionary, columns=["design_matrix_col", "factor"])
        
    # data_for_model['factor_parameter_dictionary'] = factor_parameter_dictionary

    # intercept
    data_for_model['intercept_in_design'] = data_for_model['X'].iloc[:,0].unique().tolist() == [1]

    if data_for_model['intercept_in_design'] or X_variability.design_info.term_names == [] or X_variability.design_info.term_names == ['Intercept']:
        data_for_model['A_intercept_columns'] = 1
    else:
        data_for_model['A_intercept_columns'] = data[[col for col in X_composition.design_info.term_names if col in data.columns]].drop_duplicates().shape[0]

    if data_for_model['intercept_in_design']:
        data_for_model['B_intercept_columns'] = 1
    else:
        data_for_model['B_intercept_columns'] = data[[col for col in X_composition.design_info.term_names if col in data.columns]].drop_duplicates().shape[0]
    
    # check this later
    data_for_model['user_forced_truncation_not_idx'] = list(range(1, data.shape[0]+1))

    # fill params from args
    data_for_model['prior_prec_intercept'] = prior_overdispersion_mean_association['intercept']
    data_for_model['prior_prec_slope'] = prior_overdispersion_mean_association['slope']
    data_for_model['prior_prec_sd'] = prior_overdispersion_mean_association['standard_deviation']

    data_for_model['prior_mean_intercept'] = prior_mean['intercept']
    data_for_model['prior_mean_coefficients'] = prior_mean['coefficients']

    data_for_model['exclude_priors'] = exclude_priors

    # Resolve the path to the .stan file dynamically
    stan_file_path = pkg_resources.resource_filename("sccompPy", "stan/glm_multi_beta_binomial.stan")

    # Compile Stan model
    stan_model = cmdstanpy.CmdStanModel(stan_file=stan_file_path)

    # Credible interval
    CI = 1 - (percent_false_positive / 100)
    if mcmc_seed is None:
        mcmc_seed = np.random.randint(0, int(1e5))

    # Run Stan model sampling
    fit = stan_model.sample(
        data=data_for_model,
        chains = cores,
        # quantile=CI, didnt find matching args
        # inference_method=inference_method,
        output_dir=output_directory,
        seed=mcmc_seed,
        # max_sampling_iterations=max_sampling_iterations,
        # pars=["beta", "alpha", "prec_coeff", "prec_sd", "alpha_normalised", "random_effect", "random_effect_2", "log_lik"],
        show_console=verbose
    )

    # move it here for debug
    data_for_model['factor_parameter_dictionary'] = factor_parameter_dictionary

    output = {
        'fit': fit,
        'model_input': data_for_model,
        'truncation_df2': data,
        'sample': sample,
        'cell_group': cell_group,
        'count': count,
        'formula_composition': formula_composition,
        'formula_variability': formula_variability
    }

    return output


def get_mean_precision_association(fit):
    # Extract the summary DataFrame
    summary_df = fit.summary()
    
    # Filter rows whose index starts with "prec_coeff"
    prec_coeff_values = summary_df[summary_df.index.str.startswith("prec_coeff")].Mean
        
    # Filter rows whose index starts with "prec_sd"
    prec_sd_values = summary_df[summary_df.index.str.startswith("prec_sd")].Mean

    output = {
        'prec_coeff': prec_coeff_values,
        'prec_sd' : prec_sd_values
    }
        
    return output