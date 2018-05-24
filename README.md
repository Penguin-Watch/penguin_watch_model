# penguin_watch_model

**General model structure:**

<a href="https://www.codecogs.com/eqnedit.php?latex=z_{t,i,j,k}&space;\sim&space;Binom(\phi_{t,i,j,k},&space;z_{t-1,i,j,k,})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?z_{t,i,j,k}&space;\sim&space;Binom(\phi_{t,i,j,k},&space;z_{t-1,i,j,k,})" title="z_{t,i,j,k} \sim Binom(\phi_{t,i,j,k}, z_{t-1,i,j,k,})" /></a>

<a href="http://www.codecogs.com/eqnedit.php?latex=y_{t,i,j,k}&space;\sim&space;Binom(p_{t,i,j,k},&space;z_{t,i,j,k})" target="_blank"><img src="http://latex.codecogs.com/gif.latex?y_{t,i,j,k}&space;\sim&space;Binom(p_{t,i,j,k},&space;z_{t,i,j,k})" title="y_{t,i,j,k} \sim Binom(p_{t,i,j,k}, z_{t,i,j,k})" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=logit(\phi_{t,i,j,k})&space;=&space;\mu^{\phi}&space;&plus;&space;\gamma^{\phi}_{j}&space;&plus;&space;\eta^{\phi}_{k}&space;&plus;&space;\pi^{\phi}&space;*&space;SIC_{j,k}&space;&plus;&space;\rho^{\phi}&space;*&space;KRILL_{j,k}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?logit(\phi_{t,i,j,k})&space;=&space;\mu^{\phi}&space;&plus;&space;\gamma^{\phi}_{j}&space;&plus;&space;\eta^{\phi}_{k}&space;&plus;&space;\pi^{\phi}&space;*&space;SIC_{j,k}&space;&plus;&space;\rho^{\phi}&space;*&space;KRILL_{j,k}" title="logit(\phi_{t,i,j,k}) = \mu^{\phi} + \gamma^{\phi}_{j} + \eta^{\phi}_{k} + \pi^{\phi} * SIC_{j,k} + \rho^{\phi} * KRILL_{j,k}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=logit(p_{t,i,j,k})&space;=&space;\mu^{p}&space;&plus;&space;\beta^{p}&space;*&space;x_{t}&space;&plus;&space;\nu^{p}_{j,k}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?logit(p_{t,i,j,k})&space;=&space;\mu^{p}&space;&plus;&space;\beta^{p}&space;*&space;x_{t}&space;&plus;&space;\nu^{p}_{j,k}" title="logit(p_{t,i,j,k}) = \mu^{p} + \beta^{p} * x_{t} + \nu^{p}_{j,k}" /></a>

<br>

<a href="https://www.codecogs.com/eqnedit.php?latex=z" target="_blank"><img src="https://latex.codecogs.com/gif.latex?z" title="z" /></a> = true state

<a href="https://www.codecogs.com/eqnedit.php?latex=y" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y" title="y" /></a> = observed state

<a href="https://www.codecogs.com/eqnedit.php?latex=\phi" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\phi" title="\phi" /></a> = survival probability

<a href="https://www.codecogs.com/eqnedit.php?latex=p" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p" title="p" /></a> = detection probability

<a href="https://www.codecogs.com/eqnedit.php?latex=\mu" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu" title="\mu" /></a> = intercept term

<a href="https://www.codecogs.com/eqnedit.php?latex=\gamma" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma" title="\gamma" /></a> = year effect

<a href="https://www.codecogs.com/eqnedit.php?latex=\eta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\eta" title="\eta" /></a> = site effect

<a href="https://www.codecogs.com/eqnedit.php?latex=\pi" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\pi" title="\pi" /></a> = SIC effect

<a href="https://www.codecogs.com/eqnedit.php?latex=\rho" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\rho" title="\rho" /></a> = krill catch effect

<a href="https://www.codecogs.com/eqnedit.php?latex=\beta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /></a> = slope (change in detection over time)

<a href="http://www.codecogs.com/eqnedit.php?latex=x" target="_blank"><img src="http://latex.codecogs.com/gif.latex?x" title="x" /></a> = time step within season

<a href="https://www.codecogs.com/eqnedit.php?latex=\nu" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\nu" title="\nu" /></a> = site/year effect on detection

for time step <a href="https://www.codecogs.com/eqnedit.php?latex=t" target="_blank"><img src="https://latex.codecogs.com/gif.latex?t" title="t" /></a>, nest <a href="https://www.codecogs.com/eqnedit.php?latex=i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?i" title="i" /></a>, year <a href="https://www.codecogs.com/eqnedit.php?latex=j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?j" title="j" /></a>, site <a href="https://www.codecogs.com/eqnedit.php?latex=k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k" title="k" /></a>




## Interpretation
- Site (eta) and year (gamma) effect estimates
    - Does BS differ in time and space?
- Krill (rho) and SIC (pi) effect estimates
    - Is BS impacted by abiotic/anthropogenic factors?
- Survival latent state (mean of phi)
    - When in the season is survival lowest?
    - What do these estimates mean in terms of 'chicks per nest'?
