# afdx: Diagnosis performance indicators from attributable fraction estimates.

The afdx package provides functions to estimate the attributable
fraction using logit exponential model or bayesian latent class model.

## The logit exponential model

The logitexp function estimated the logit exponential function fitting a
maximum likelihood model. The senspec() function estimate the
sensitivity, specificity, positive predicted value and negative
predicted values for the specified cut-off points.

## The bayesian latent class model

The get_latent_model() provides an rjags model template to estimate the
attributable fraction and the sensitivity, specificity, positive
predicted value and negative predicted value of the latent class model.

@docType package @name afdx

## See also

Useful links:

- <https://github.com/johnaponte/afdx>

## Author

**Maintainer**: John J. Aponte <john.j.aponte@gmail.com>
([ORCID](https://orcid.org/0000-0002-3014-3673))

Authors:

- Orvalho Augusto <caveman@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-0005-3968))
