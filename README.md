# penguin_watch
Repo for Penguin Watch mark-recapture survivorship model.

The most recent working script for model is `HPC/model_script.R`. Tweaks are continually being made. 

All model runs being performed on high performance computing cluster at SBU. PBS script for cluster is `HPC/pwatch.pbs`. Model output is located in `HPC/Archive`, though some data objects are too large to be commited to git and were transfered to local machine from cluster via `scp`.

**General model structure:**

![](model_notes_SESYNC.JPG)
