# penguin_watch
Repo for Penguin Watch survivorship model.

All model runs being performed on high performance computing cluster.

## Interpretation
- Site (eta) and year (gamma) effect estimates
    - Does BS differ in time and space?
- Krill (rho) and SIC (pi) effect estimates
    - Is BS impacted by abiotic/anthropogenic factors?
- Survival latent state (mean of phi)
    - When in the season is survival lowest?
    - What do these estimates mean in terms of 'chicks per nest'?

**General model structure:**

$z_{t,i,j,k} \sim Binom(\phi_{t,i,j,k}, z_{t-1,i,j,k,})$
$y_{t,i,j,k} \sim Binom(p_{t,i,j,k} * w_{t,i,j,k}, z__{t,i,j,k})$
$logit(phi_{t,i,j,k}) = \mu_{\phi} + \gamma_{\phi_{j}} + \eta_{\phi_{k}} + \pi_{\phi} + \rho_{\phi} * KRILL_{j,k}$
$logit(p_{t,i,j,k}) <- \mu_{p} + \beta_{p} * x_{t} + \nu_{p_{j,k}}$
