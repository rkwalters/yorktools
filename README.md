# YorkTools

R package for York regression with extended options (e.g. fixed intercept)

# Installation:

```
require(devtools)
devtools::install_github("rkwalters/yorktools", ref="main")
```

Currently requires a PAT (personal access token) to install from private repo. Instructions for creating a PAT can be found [here](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token). You can then either supply it to `install_github` with the argument `auth_token` or make it available with environment variable `GITHUB_PAT`, e.g. by adding it to your `~/.Rprofile`:

```
Sys.setenv(GITHUB_PAT="YOUR_PAT_HERE")
```

# Usage

See `?york`

# To do:

* Clean up git history before release
* add ML SEs
* add predict function
* add resid function
* add coef function
* add normal (not fixed) intercept

