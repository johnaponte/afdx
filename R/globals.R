# Globals used within functions without explicit external binding
# 20210201 by JJAV

# logitexp function
utils::globalVariables(c(".","rr"))

# sens_spec functions
utils::globalVariables(c("v.fever","v.density","tnc","rrc"))

# make_cutoff function
utils::globalVariables(c("0", "denscat", "density", "fever", "mincat", "todel"))