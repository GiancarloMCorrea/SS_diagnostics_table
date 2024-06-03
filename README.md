# SS_diagnostics_table

Get diagnostics from a SS model grid. Based on [Carvalho et al. (2021)](https://doi.org/10.1016/j.fishres.2021.105959). 

The diagnostics that are obtained are:

- Hessian invertible?
- Negative log-likelihood
- Number of parameters estimated
- Number of parameters on bounds
- Maximum gradient
- Run test (CPUE)
- Run test (mean length)
- Root mean square error (CPUE)
- Root mean square error (mean length)
- Mean absolute scaled error (CPUE)
- Mean absolute scaled error (mean length)
- Mohn's rho (SSB)
- Mohn's rho (F)
- Trends in recruitment deviates (see [Merino et al. 2022](https://doi.org/10.1016/j.fishres.2022.106478))

The *run_diagnostics.R* is the main script. Follow the instructions therein. 

If you want to generate two tables in Word to paste it in your report, use the *Diagnostic-table.qmd* file (it uses Quarto). 
