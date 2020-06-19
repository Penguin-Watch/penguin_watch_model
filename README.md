# penguin_watch_model

Code for Youngflesh et al. *Accepted.* *Remote Sensing in Ecology and Conservation*


Associated publication:

Youngflesh, C.\*, Jones, F.M.\*, Lynch, H.J., Arthur, J., Macháčková, Z., Tosey, H.R., Hart, T. *Accepted.* Large-scale assessment of intra- and inter-annual nesting success using a remote camera network. __*Remote Sensing in Ecology and Conservation*__. [link](URL_HERE)

\* _Authors contributed equally_


Associated data available here:

DRYAD CITATION [link](URL_HERE)


Repository structure:

* `Scripts/`
  * `1-process-pw-data.R` - Prepare data for modeling
  * `2-model.R` - R script for capture-recapture model
  * `2-run-model.sh` - Shell script to run model on HPC resources
  * `3-process-output.R` - Process model output and precipitation data
  * `4-process-krill-data.R` - Process krill data
  * `5-process-tourism-data.R` - Process tourism data (available upon request from IAATO)
  * `6-analyze-output.R ` - Posterior predictive checks, create breeding success plots, breeding success as a function of covariates
  * `7-mortality-timing.R` - Quantify the timing of mortality events
  * `8-krill-plots.R` - Create krill catch plots
  * `9-gif-precip.R` - Create gifs with cameras images and model results
