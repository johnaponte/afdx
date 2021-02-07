# Globals used within functions without explicit external binding
# 20210201 by JJAV

# logitexp function
utils::globalVariables(c(".", "rr"))

# sens_spec functions
utils::globalVariables(c("v.fever", "v.density", "tnc", "rrc"))

# make_cutoff function
utils::globalVariables(c("0", "denscat", "density", "fever", "mincat", "todel"))

# make_n_cutoff function
utils::globalVariables(
  c(
    "chg2fev",
    "chg2nofev",
    "chgfev",
    "chgfev_lag",
    "chgnofev",
    "chgnofev_lag",
    "cumfev",
    "cumfev_lag",
    "cumnofev",
    "cumnofev_lag",
    "cumtotal",
    "difpart",
    "intpart",
    "k",
    "lag",
    "lagpart",
    "no_fever",
    "total"
  )
)