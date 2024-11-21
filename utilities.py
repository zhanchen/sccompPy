import pandas as pd
import re
import numpy as np
from functools import reduce

#' draws_to_tibble_x_y
#'
#'
#' @param fit A fit object
#' @param par A character vector. The parameters to extract.
#' @param x A character. The first index.
#' @param y A character. The first index.
#'
#' @keywords internal
#' @noRd
def draws_to_tibble_x_y(fit, par, x, y, number_of_draws = None):

    # Extract parameter names that match the specified pattern
    par_names = [var for var in fit.stan_variables().keys() if re.search(par, var)]

    # Extract draws for the specified parameter and convert to DataFrame format
    draws_df = fit.draws_pd([par, 'chain__', 'iter__', 'draw__'])

    # Pivot longer to reshape the DataFrame, renaming and selecting relevant columns
    draws_long = draws_df.melt(var_name="parameter", value_name="value", value_vars=[col for col in draws_df.columns if par in col], id_vars = ['chain__', 'iter__', 'draw__'])

    # Extract chain, variable, x, and y indices from the parameter string
    pattern = r"([1-9]+)?\.?([a-zA-Z0-9_\.]+)\[([0-9]+),([0-9]+)"
    draws_long[['chain', 'variable', x, y]] = draws_long['parameter'].str.extract(pattern)

    # Convert extracted x and y values to integers
    draws_long[x] = draws_long[x].astype(int)
    draws_long[y] = draws_long[y].astype(int)

    # Sort and prepare the DataFrame
    draws_long = draws_long.sort_values(by=['variable', x, y, 'chain__'])
    draws_long = draws_long.groupby(['variable', x, y]).apply(lambda df: df.assign(draw__=range(1, len(df) + 1))).reset_index(drop=True)

    # Select relevant columns and filter by the parameter of interest
    draws_long = draws_long[['chain__', 'iter__', 'draw__', 'variable', x, y, 'value']]
    draws_long = draws_long[draws_long['variable'] == par]

    return draws_long

def summary_to_tibble(fit, par, x, y = None, probs = (5, 25, 50, 75, 95)):
    
    # Extract parameter names matching 'par', DOES NOT compatible with cmdstanpy
    # par_names = [name for name in fit.column_names if re.search(par, name)]

    # Get the summary DataFrame
    summary = fit.summary(probs)

    # Select only the desired columns
    filtered_summary = summary.loc[summary.index.str.contains(par)]

    # Convert probs to column names like "5%", "25%", etc.
    prob_cols = [f"{p}%" for p in probs] 
    columns_to_keep = ["Mean"] + prob_cols + ['N_Eff', 'R_hat']

    filtered_summary = filtered_summary[columns_to_keep]

    # Reset the index and store the old index in a new column called 'variable'
    filtered_summary = filtered_summary.reset_index().rename(columns={"index": "variable"})

    # Split variable names into '.variable', x, and (optional) y
    def split_variable_name(var_name):
        parts = re.split(r"\[|\]|,|\s", var_name)
        parts = [p for p in parts if p]  # Remove empty strings
        if y:
            return parts[:3]  # Keep .variable, x, and y
        else:
            return parts[:2]  # Keep .variable and x
    
    # Apply splitting logic
    split_cols = filtered_summary["variable"].apply(split_variable_name)
    split_df = pd.DataFrame(split_cols.tolist(), columns=["variable", x] + ([y] if y else []))

    # Ensure x and y are integers
    split_df[x] = split_df[x].astype(int)
    if y:
        split_df[y] = split_df[y].astype(int)

    filtered_summary = pd.concat([split_df, filtered_summary.drop(columns=['variable'])], axis=1)

    # Add missing columns if not present
    if "N_Eff" not in filtered_summary.columns:
        filtered_summary["N_Eff"] = np.nan
    if "R_hat" not in filtered_summary.columns:
        filtered_summary["R_hat"] = np.nan
    
    return filtered_summary



def mutate_from_expr_list(x, formula_expr, ignore_errors = True):

    # Ensure formula_expr keys have names
    if not formula_expr or not isinstance(formula_expr, dict):
        formula_expr = {k: k for k in formula_expr}
    
    # Check if all elements of contrasts are in the parameter
    parameter_names = x.columns.tolist()

    # Process contrast elements by removing fractions, decimals, and splitting expressions
    def clean_contrast_elements(expr):
        expr = re.sub(r"[0-9]+/[0-9]+\s?\*", "", expr)  # Remove fractions
        expr = re.sub(r"[-+]?[0-9]+\.[0-9]+\s?\*", "", expr)  # Remove decimals
        elements = re.split(r"[+\-*]", expr)  # Split by operators
        elements = [re.sub(r"[\(\)\s]", "", e) for e in elements]  # Remove parentheses and spaces
        return elements

    contrast_elements = reduce(lambda a, b: a + b, [clean_contrast_elements(f) for f in formula_expr.values()], [])

    # Check if backquotes are required (columns with special characters)
    def requires_backquotes(element):
        return not re.match(r"^[a-zA-Z0-9_]+$", element)  # Only valid chars for column names

    invalid_contrasts = [e for e in contrast_elements if requires_backquotes(e) and e not in parameter_names]
    if invalid_contrasts:
        print(f"Warning: These elements require backquotes: {invalid_contrasts}")

    # Check if contrast columns exist in the DataFrame
    missing_contrasts = [e for e in contrast_elements if e not in parameter_names]
    if missing_contrasts and not ignore_errors:
        raise ValueError(f"These contrasts are not present in the DataFrame: {missing_contrasts}")

    # Function to escape column names with backticks if they contain special characters
    def escape_column_names(formula, columns):
        for col in columns:
            if not re.match(r"^[a-zA-Z_][a-zA-Z0-9_]*$", col):  # Matches valid Python identifiers
                escaped_col = f"`{col}`"
                formula = formula.replace(col, escaped_col)
        return formula

    # Apply formulas to mutate the DataFrame
    def mutate_with_formula(df, column_name, formula):
        try:
            # Escape column names with special characters
            formula = escape_column_names(formula, parameter_names)
            df[column_name] = df.eval(formula)
        except Exception as e:
            if not ignore_errors:
                raise e
            print(f"Warning: Error while processing formula '{formula}': {e}")
        return df

    for new_col, formula in formula_expr.items():
        x = mutate_with_formula(x, new_col, formula)

    # # Add columns not in formula_expr to the final DataFrame
    # remaining_columns = [col for col in x.columns if col not in formula_expr.keys()]
    # x = x[formula_expr.keys() + remaining_columns]

    return x

def get_abundance_contrast_draws(
    data, 
    contrasts=None
):
    # Retrieve attributes from data
    cell_group = data.get("cell_group", None)
    model_input = data.get("model_input", {})
    fit = data.get("fit", {})

    # Beta
    beta_factor_of_interest = data.get('model_input').get("X").columns.tolist()
    beta = draws_to_tibble_x_y(fit, "beta", "C", "M")
    beta = beta.pivot(index = ['chain__', 'iter__', 'draw__', 'variable', 'M'], columns='C', values='value')
    beta.columns = beta_factor_of_interest
    beta.reset_index(inplace=True)

    # Abundance
    draws = beta.drop(columns=['variable'], errors='ignore')

    # Random effect
    n_random_eff = model_input.get('n_random_eff', 0)
    ##### will test this later when developing random effects ###
    if n_random_eff > 0:
        beta_random_effect_factor_of_interest = model_input.get("X_random_effect", {}).columns if model_input.get("X_random_effect", None) is not None else []
        beta_random_effect = draws_to_tibble_x_y(fit, "random_effect", "C", "M")
        
        # Add last component
        beta_random_effect_sum = beta_random_effect.groupby(['C', 'chain__', 'iter__', 'draw__', 'variable'])['.value'].sum().reset_index()
        beta_random_effect_sum['.value'] = -beta_random_effect_sum['.value']
        beta_random_effect_sum['M'] = beta_random_effect['M'].max() + 1
        beta_random_effect = pd.concat([beta_random_effect, beta_random_effect_sum])
        
        # Reshape and merge
        beta_random_effect = beta_random_effect.pivot(index=['chain__', 'iter__', 'draw__'], columns='C', values='.value')
        beta_random_effect.columns = beta_random_effect_factor_of_interest
        draws = draws.merge(beta_random_effect, left_index=True, right_index=True, how='left')
    
    # Second random effect
    if n_random_eff > 1:
        beta_random_effect_factor_of_interest_2 = model_input.get("X_random_effect_2", {}).columns if model_input.get("X_random_effect_2", None) is not None else []
        beta_random_effect_2 = draws_to_tibble_x_y(fit, "random_effect_2", "C", "M")
        
        # Add last component
        beta_random_effect_2_sum = beta_random_effect_2.groupby(['C', 'chain__', 'iter__', 'draw__', '.variable'])['value'].sum().reset_index()
        beta_random_effect_2_sum['.value'] = -beta_random_effect_2_sum['.value']
        beta_random_effect_2_sum['M'] = beta_random_effect_2['M'].max() + 1
        beta_random_effect_2 = pd.concat([beta_random_effect_2, beta_random_effect_2_sum])
        
        # Reshape and merge
        beta_random_effect_2 = beta_random_effect_2.pivot(index=['chain__', 'iter__', 'draw__'], columns='C', values='value')
        beta_random_effect_2.columns = beta_random_effect_factor_of_interest_2
        draws = draws.merge(beta_random_effect_2, left_index=True, right_index=True, how='left')
    
    # Apply contrasts if specified
    # develop this later, require an example of contrasts
    if contrasts:
        draws = mutate_from_expr_list(draws, contrasts, ignore_errors = True)

    # Convergence diagnostics
    convergence_df = summary_to_tibble(fit, "beta", "C", "M") 

    # Attach cell names
    y = model_input.get("y", None)
    
    cell_names = pd.Series(y.columns, name=cell_group)
    cell_names.index += 1  # Adjust for 1-based indexing in R
    
    draws = draws.merge(cell_names, left_on ='M', right_index=True, how='left')
    draws = draws[[cell_group] + [col for col in draws.columns if col != cell_group]]

    factor_names = pd.Series(beta_factor_of_interest, name="parameter")
    factor_names.index += 1
    
    convergence_df = convergence_df.merge(cell_names, left_on ='M', right_index=True, how='left')
    convergence_df = convergence_df.merge(factor_names, left_on ='C', right_index=True, how='left')
    
    # Check if no contrasts of interest; if so, return minimal output
    if len(draws.columns) <= 5:
        return draws[[cell_group, 'M']].drop_duplicates()
     
    if 'R_hat' in convergence_df.columns:
        convergence_df.rename(columns={'R_hat': 'R_k_hat'}, inplace=True)
    elif 'k_hat' in convergence_df.columns:
        convergence_df.rename(columns={'k_hat': 'R_k_hat'}, inplace=True)
    
    # Attach convergence info and reformat
    draws = draws.melt(id_vars=[cell_group, 'M', 'chain__', 'iter__', 'draw__'], var_name='parameter', value_name='value')
    draws = draws.merge(convergence_df[[cell_group, 'parameter', 'N_Eff', 'R_k_hat']], on=[cell_group, 'parameter'], how='left')

    return draws.sort_values('parameter')

def get_FDR(x):
    """
    Calculate the False Discovery Rate (FDR) for a series of values.

    Parameters:
        x (list or pd.Series): Input values.

    Returns:
        pd.Series: The calculated FDR values.
    """
    # Convert input to a DataFrame
    df = pd.DataFrame({'value': x})
    
    # Sort by value
    df = df.sort_values(by='value').reset_index()
    
    # Calculate cumulative mean (FDR)
    df['FDR'] = df['value'].expanding().mean()
    
    # Restore original order
    df = df.sort_values(by='index')
    
    return df['FDR'].values

def draws_to_statistics(draws, false_positive_rate, test_composition_above_logit_fold_change, cell_group, prefix=""):
    
    """
    Calculate summary statistics for draws.

    Parameters:
        draws (pd.DataFrame): The input DataFrame with draws.
        false_positive_rate (float): False positive rate for confidence intervals.
        test_composition_above_logit_fold_change (float): Threshold for determining significance.
        cell_group (str): The column representing the cell group.
        prefix (str): Prefix for the resulting column names.

    Returns:
        pd.DataFrame: The DataFrame with calculated statistics.
    """
    # Summarize the statistics
    grouped = draws.groupby([cell_group, 'M', 'parameter'])

    summary = grouped.apply(lambda group: pd.Series({
        f"{prefix}lower": group["value"].quantile(false_positive_rate / 2),
        f"{prefix}effect": group["value"].quantile(0.5),
        f"{prefix}upper": group["value"].quantile(1 - (false_positive_rate / 2)),
        f"{prefix}bigger_zero": (group["value"] > test_composition_above_logit_fold_change).sum(),
        f"{prefix}smaller_zero": (group["value"] < -test_composition_above_logit_fold_change).sum(),
        "R_k_hat": group["R_k_hat"].unique()[0] if "R_k_hat" in group else None,
        "N_Eff": group["N_Eff"].unique()[0] if "N_Eff" in group else None,
        "n": len(group)
    })).reset_index()

    # Calculate pH0
    summary[f"{prefix}pH0"] = 1 - summary[[f"{prefix}bigger_zero", f"{prefix}smaller_zero"]].max(axis=1) / summary["n"]

    # Calculate FDR for each parameter
    summary["FDR"] = summary.groupby("parameter")[f"{prefix}pH0"].transform(lambda pH0: get_FDR(pH0))

    # Select relevant columns
    relevant_cols = [cell_group, "M", "parameter", f"{prefix}lower", f"{prefix}effect", f"{prefix}upper", f"{prefix}pH0", "FDR"]
    if "N_Eff" in summary.columns:
        relevant_cols.append("N_Eff")
    if "R_k_hat" in summary.columns:
        relevant_cols.append("R_k_hat")

     # Return the result
    return summary[relevant_cols]

def get_variability_contrast_draws(data, contrasts=None):
    """
    Extract variability contrast draws from data.

    Parameters:
        data (dict): Input data containing attributes and model fits.
        contrasts (list or None): List of contrast expressions to apply.

    Returns:
        pd.DataFrame: A DataFrame with variability contrast draws.
    """
    # Extract attributes
    cell_group = data.get("cell_group")
    model_input = data.get("model_input")
    fit = data.get("fit")

    # Variability factors
    variability_factor_of_interest = model_input.get("XA").columns if "XA" in model_input else []

    # Get initial draws
    draws = draws_to_tibble_x_y(fit, "alpha_normalised", "C", "M")
    draws["value"] = -draws["value"]  # Invert values for variability

    # Pivot wider and rename columns
    draws = draws.pivot(index=["M", "chain__", "iter__", "draw__"], columns="C", values="value")
    draws.columns = variability_factor_of_interest  # Rename columns
    draws = draws.reset_index().drop(columns=["variable"], errors="ignore") 

    # Apply contrasts if provided
    if contrasts:
        draws = mutate_from_expr_list(draws, contrasts, ignore_errors=True)

    # Add cell group names
    cell_names = pd.DataFrame({
        "M": range(1, len(model_input["y"].columns) + 1),
        cell_group: model_input["y"].columns
    })
    draws = draws.merge(cell_names, on="M", how="left")
    draws = draws[[cell_group] + [col for col in draws.columns if col != cell_group]]

    # If no contrasts of interest, return minimal output
    if draws.shape[1] == 5:  # Includes M and cell_group
        draws.drop_duplicates(subset=["M", cell_group])

    # Get convergence
    convergence_df = summary_to_tibble(fit, "alpha_normalised", "C", "M")

    # Add cell names to convergence_df
    convergence_df = convergence_df.merge(cell_names, on="M", how="left")

    # Add factor names to convergence_df
    factor_names = pd.DataFrame({
        "C": range(1, len(variability_factor_of_interest) + 1),
        "parameter": variability_factor_of_interest
    })
    convergence_df = convergence_df.merge(factor_names, on="C", how="left")

    # Rename convergence columns
    if "R_hat" in convergence_df.columns:
        convergence_df.rename(columns={"R_hat": "R_k_hat"}, inplace=True)
    elif "k_hat" in convergence_df.columns:
        convergence_df.rename(columns={"k_hat": "R_k_hat"}, inplace=True)

    convergence_df = convergence_df[[cell_group, "parameter", "N_Eff", "R_k_hat"]].dropna()

    # Reshape draws to long format
    draws_long = draws.melt(
       id_vars=["M", "chain__", "iter__", "draw__", cell_group],
       var_name="parameter",
       value_name="value"
    )

    # Attach convergence
    draws_long = draws_long.merge(convergence_df, on=[cell_group, "parameter"], how="left")

    # Reorder parameters
    parameter_order = draws.columns[5:]  # Columns beyond first 5 are parameters
    draws_long["parameter"] = pd.Categorical(draws_long["parameter"], categories=parameter_order, ordered=True)
    draws_long = draws_long.sort_values(["parameter", 'M'])

    return draws_long


