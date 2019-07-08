# penguin_watch_model

**General model structure:**

<a href="https://www.codecogs.com/eqnedit.php?latex=z_{t,i,j,k}&space;\sim&space;Binom(\phi_{t,i,j,k},&space;z_{t-1,i,j,k})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?z_{t,i,j,k}&space;\sim&space;Binom(\phi_{t,i,j,k},&space;z_{t-1,i,j,k})" title="z_{t,i,j,k} \sim Binom(\phi_{t,i,j,k}, z_{t-1,i,j,k})" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=y_{t,i,j,k}&space;\sim&space;Binom(p_{t,i,j,k},&space;z_{t,i,j,k})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y_{t,i,j,k}&space;\sim&space;Binom(p_{t,i,j,k},&space;z_{t,i,j,k})" title="y_{t,i,j,k} \sim Binom(p_{t,i,j,k}, z_{t,i,j,k})" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=logit(\phi_{t,i,j,k})&space;=&space;\mu_{\phi_{j,k}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?logit(\phi_{t,i,j,k})&space;=&space;\mu_{\phi_{j,k}}" title="logit(\phi_{t,i,j,k}) = \mu_{\phi_{j,k}}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=logit(p_{t,i,j,k})&space;=&space;\mu_{p}&space;&plus;&space;\nu_{p_{i,j,k}}&space;&plus;&space;\beta_{p_{j,k}}&space;\times&space;x_{t}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?logit(p_{t,i,j,k})&space;=&space;\mu_{p}&space;&plus;&space;\nu_{p_{i,j,k}}&space;&plus;&space;\beta_{p_{j,k}}&space;\times&space;x_{t}" title="logit(p_{t,i,j,k}) = \mu_{p} + \nu_{p_{i,j,k}} + \beta_{p_{j,k}} \times x_{t}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\mu_{\phi}&space;\sim&space;N(\theta_{\phi},&space;\sigma_{\mu_{\phi}})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu_{\phi}&space;\sim&space;N(\theta_{\phi},&space;\sigma_{\mu_{\phi}})" title="\mu_{\phi} \sim N(\theta_{\phi}, \sigma_{\mu_{\phi}})" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\beta_{p_{j,k}}&space;\sim&space;N(\mu_{\beta_{p}},&space;\sigma_{\beta_{p}})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta_{p_{j,k}}&space;\sim&space;N(\mu_{\beta_{p}},&space;\sigma_{\beta_{p}})" title="\beta_{p_{j,k}} \sim N(\mu_{\beta_{p}}, \sigma_{\beta_{p}})" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\nu_{p_{i,j,k}}&space;\sim&space;N(0,&space;\sigma_{\nu_{p}})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\nu_{p_{i,j,k}}&space;\sim&space;N(0,&space;\sigma_{\nu_{p}})" title="\nu_{p_{i,j,k}} \sim N(0, \sigma_{\nu_{p}})" /></a>


**Derived quantities:**

<a href="https://www.codecogs.com/eqnedit.php?latex=Z_{t,j,k}&space;=&space;\forall&space;\&space;t,j,k;&space;\sum_{i}{}&space;z_{t,i,j,k}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Z_{t,j,k}&space;=&space;\forall&space;\&space;t,j,k;&space;\sum_{i}{}&space;z_{t,i,j,k}" title="Z_{t,j,k} = \forall \ t,j,k; \sum_{i}{} z_{t,i,j,k}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\overline{p_{t,j,k}}&space;=&space;\forall&space;\&space;t,j,k;&space;\frac{1}{I}&space;\sum_{i}{}&space;p_{t,i,j,k}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\overline{p_{t,j,k}}&space;=&space;\forall&space;\&space;t,j,k;&space;\frac{1}{I}&space;\sum_{i}{}&space;p_{t,i,j,k}" title="\overline{p_{t,j,k}} = \forall \ t,j,k; \frac{1}{I} \sum_{i}{} p_{t,i,j,k}" /></a>


<br>

<a href="https://www.codecogs.com/eqnedit.php?latex=z" target="_blank"><img src="https://latex.codecogs.com/gif.latex?z" title="z" /></a> = true state

<a href="https://www.codecogs.com/eqnedit.php?latex=y" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y" title="y" /></a> = observed state

<a href="https://www.codecogs.com/eqnedit.php?latex=\phi" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\phi" title="\phi" /></a> = survival probability

<a href="https://www.codecogs.com/eqnedit.php?latex=p" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p" title="p" /></a> = detection probability

<a href="https://www.codecogs.com/eqnedit.php?latex=\mu_{\phi}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu_{\phi}" title="\mu_{\phi}" /></a> = intercept for survival probability (varies by site/year)

<a href="https://www.codecogs.com/eqnedit.php?latex=\mu_{p}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu_{p}" title="\mu_{p}" /></a> = intercept term for detection probability (mean detection for all nests/sites/years at mean x)

<a href="https://www.codecogs.com/eqnedit.php?latex=\beta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /></a> = slope (change in detection over time); varies by site/year

<a href="http://www.codecogs.com/eqnedit.php?latex=x" target="_blank"><img src="http://latex.codecogs.com/gif.latex?x" title="x" /></a> = time step within season (daily)

<a href="https://www.codecogs.com/eqnedit.php?latex=\nu" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\nu" title="\nu" /></a> = nest/site/year effect on detection

<a href="https://www.codecogs.com/eqnedit.php?latex=Z" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Z" title="Z" /></a> = number of chicks (true state) for each time step/site/year (across all nests)

<a href="https://www.codecogs.com/eqnedit.php?latex=\bar{p}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\bar{p}" title="\bar{p}" /></a> = mean detection probability for each time step/site/year (across all nests)

<a href="https://www.codecogs.com/eqnedit.php?latex=t" target="_blank"><img src="https://latex.codecogs.com/gif.latex?t" title="t" /></a> = time step

<a href="https://www.codecogs.com/eqnedit.php?latex=i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?i" title="i" /></a> = nest

<a href="https://www.codecogs.com/eqnedit.php?latex=j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?j" title="j" /></a> = year

<a href="https://www.codecogs.com/eqnedit.php?latex=k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k" title="k" /></a> = site
