import pandas as pd
import numpy as np
from .utilities import *
from .functions_multi_beta_binomia import *

def sccomp_estimate(
    data,
    formula_composition="~1",
    formula_variability="~1",
    sample=None,
    cell_group=None,
    count=None,

    # Secondary arguments
    cores=4,  # Replace `detectCores()` with a default number of cores
    bimodal_mean_variability_association=False,
    percent_false_positive=5,
    inference_method="pathfinder",
    prior_mean={"intercept": [0, 1], "coefficients": [0, 1]},
    prior_overdispersion_mean_association={
        "intercept": [5, 2],
        "slope": [0, 0.6],
        "standard_deviation": [10, 20]
    },
    sample_cell_group_pairs_to_exclude=None,
    output_directory="sccomp_draws_files",
    verbose=True,
    enable_loo=False,
    noise_model="multi_beta_binomial",
    exclude_priors=False,
    use_data=True,
    mcmc_seed=None,
    max_sampling_iterations=20000,
    pass_fit=True,
    count_col=None,
    approximate_posterior_inference=None,
    variational_inference=None,
    **kwargs
):
    res = sccomp_glm_data_frame_counts(
        data,
        formula_composition=formula_composition,
        formula_variability=formula_variability,
        sample=sample,
        cell_group=cell_group,
        count=count,
        cores=cores,
        bimodal_mean_variability_association=bimodal_mean_variability_association,
        # percent_false_positive=percent_false_positive,
        inference_method=inference_method,
        prior_mean=prior_mean,
        prior_overdispersion_mean_association=prior_overdispersion_mean_association,
        sample_cell_group_pairs_to_exclude=sample_cell_group_pairs_to_exclude,
        output_directory=output_directory,
        verbose=verbose,
        enable_loo=enable_loo,
        exclude_priors=exclude_priors,
        use_data=use_data,
        mcmc_seed=mcmc_seed,
        max_sampling_iterations=max_sampling_iterations,
        pass_fit=pass_fit,
        **kwargs
        )

    # Informational message
    if verbose:
        print(
            "sccomp says: to do hypothesis testing run `sccomp_test()`. "
            "The `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, "
            "and 0.7 equates to ~100% increase, if the baseline is ~0.1 proportion. "
            "Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear)."
        )

    # Add noise model as an attribute to the result
    res["noise_model"] = noise_model

    return res



def sccomp_test(data, contrasts=None, percent_false_positive=5, test_composition_above_logit_fold_change=0.1, pass_fit=True):
    # Extract attributes
    fit = data.get('fit')
    sample = data.get("sample")
    cell_group = data.get("cell_group")
    count = data.get("count")
    model_input = data.get("model_input")
    truncation_df2 = data.get("truncation_df2")
    formula_composition = data.get("formula_composition")
    formula_variability = data.get("formula_variability")


    # Abundance
    abundance_CI = get_abundance_contrast_draws(data, contrasts)

    # If contrasts exist in abundance_CI
    if "parameter" in abundance_CI.columns:
        abundance_CI = draws_to_statistics(
            abundance_CI,
            percent_false_positive / 100,
            test_composition_above_logit_fold_change,
            cell_group,
            "c_"
        )

    # Variability
    variability_CI = get_variability_contrast_draws(data, contrasts)
    if "parameter" in variability_CI.columns:
        variability_CI = draws_to_statistics(
            variability_CI,
            percent_false_positive / 100,
            test_composition_above_logit_fold_change,
            cell_group,
            "v_"
        )

    # Factor parameter dictionary
    if "factor" not in model_input.get("factor_parameter_dictionary", {}).columns:
        factor_parameter_dictionary = pd.DataFrame(columns=["factor", "design_matrix_col"])
    else:
        factor_parameter_dictionary = model_input["factor_parameter_dictionary"][["factor", "design_matrix_col"]]

    # Merge results
    result = abundance_CI.merge(variability_CI, how="left")

    # Add factor labels
    result = result.merge(
        factor_parameter_dictionary,
        left_on="parameter",
        right_on="design_matrix_col",
        how="left"
    ).drop(columns=["M"], errors="ignore")

    # Add truncation data
    truncation_data = truncation_df2.drop(
        columns=["M", "N", ".variable", "mean", "se_mean", "sd", "N_Eff", "R_hat", "k_hat", "Rhat", ".lower", ".median", ".upper"],
        errors="ignore"
    ).groupby(cell_group).apply(
        lambda df: pd.Series({"count_data": df})
    )
    result = result.merge(
        truncation_data, 
        left_on=cell_group,  
        right_index=True,
        how="left"
    )

    output = {
        'result': result,
        'mean_concentration_association': get_mean_precision_association(fit),
        'test_composition_above_logit_fold_change': test_composition_above_logit_fold_change,
        'model_input': model_input,
        'truncation_df2': truncation_df2,
        'noise_model': data.get("noise_model"),
        'sample': sample,
        'cell_group': cell_group,
        'count': count,
        'formula_composition': formula_composition,
        'formula_variability': formula_variability
    }

    if pass_fit:
        output['fit'] = fit

    return output
