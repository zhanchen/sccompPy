import pandas as pd
import numpy as np

def sccomp_test(data, contrast=None, percent_false_positive=5, test_composition_above_logit_fold_change=0.1, pass_fit=True):
    # Extract attributes
    sample = data.get("sample")
    cell_group = data.get("cell_group")
    count = data.get("count")
    model_input = data.get("model_input")
    truncation_df2 = data.get("truncation_df2")

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

    