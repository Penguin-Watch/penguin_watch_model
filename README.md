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

<a href="https://www.codecogs.com/eqnedit.php?latex=z_{t,i,j,k}&space;\sim&space;Binom(\phi_{t,i,j,k},&space;z_{t-1,i,j,k,})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?z_{t,i,j,k}&space;\sim&space;Binom(\phi_{t,i,j,k},&space;z_{t-1,i,j,k,})" title="z_{t,i,j,k} \sim Binom(\phi_{t,i,j,k}, z_{t-1,i,j,k,})" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=y_{t,i,j,k}&space;\sim&space;Binom(p_{t,i,j,k}&space;*&space;w_{t,i,j,k},&space;z_{t,i,j,k})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y_{t,i,j,k}&space;\sim&space;Binom(p_{t,i,j,k}&space;*&space;w_{t,i,j,k},&space;z_{t,i,j,k})" title="y_{t,i,j,k} \sim Binom(p_{t,i,j,k} * w_{t,i,j,k}, z_{t,i,j,k})" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=logit(\phi_{t,i,j,k})&space;=&space;\mu^{\phi}&space;&plus;&space;\gamma^{\phi}_{j}&space;&plus;&space;\eta^{\phi}_{k}&space;&plus;&space;\pi^{\phi}&space;*&space;SIC_{j,k}&space;&plus;&space;\rho^{\phi}&space;*&space;KRILL_{j,k}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?logit(\phi_{t,i,j,k})&space;=&space;\mu^{\phi}&space;&plus;&space;\gamma^{\phi}_{j}&space;&plus;&space;\eta^{\phi}_{k}&space;&plus;&space;\pi^{\phi}&space;*&space;SIC_{j,k}&space;&plus;&space;\rho^{\phi}&space;*&space;KRILL_{j,k}" title="logit(\phi_{t,i,j,k}) = \mu^{\phi} + \gamma^{\phi}_{j} + \eta^{\phi}_{k} + \pi^{\phi} * SIC_{j,k} + \rho^{\phi} * KRILL_{j,k}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=logit(p_{t,i,j,k})&space;=&space;\mu^{p}&space;&plus;&space;\beta^{p}&space;*&space;x_{t}&space;&plus;&space;\nu^{p}_{j,k}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?logit(p_{t,i,j,k})&space;=&space;\mu^{p}&space;&plus;&space;\beta^{p}&space;*&space;x_{t}&space;&plus;&space;\nu^{p}_{j,k}" title="logit(p_{t,i,j,k}) = \mu^{p} + \beta^{p} * x_{t} + \nu^{p}_{j,k}" /></a>
