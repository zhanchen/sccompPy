import pandas as pd
import re

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
    draws_df = fit.draws_pd(['beta', 'chain__', 'iter__', 'draw__'])

    # Pivot longer to reshape the DataFrame, renaming and selecting relevant columns
    draws_long = draws_df.melt(var_name="parameter", value_name="value", value_vars=[col for col in draws_df.columns if 'beta' in col], id_vars = ['chain__', 'iter__', 'draw__'])

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
