import pandas as pd
import numpy as np

def get_abundance_contrast_draws(
    data, 
    contrasts=None
):
    # Retrieve attributes from data
    cell_group = data.get(".cell_group", None)
    model_input = data.get("model_input", {})
    fit = data.get("fit", {})

    # Extract parameter names that match the specified pattern
    beta_factor_of_interest = data.get("X", {}).columns if data.get("X", None) is not None else []