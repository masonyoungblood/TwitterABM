# TwitterABM
*Agent-Based Model of Biased Cultural Transmission on Twitter*

This package includes an agent-based model (ABM) of cultural transmission on Twitter that simulates the three main forms of transmission bias: content bias, frequency bias, and demonstrator bias. The ABM has elements from [Carrignon et al. (2019)](https://www.nature.com/articles/s41599-019-0295-9), [Lachlan et al. (2018)](https://www.nature.com/articles/s41467-018-04728-1), and [Youngblood & Lahti (2021)](https://www.biorxiv.org/content/10.1101/2021.03.05.434109v1). 

Install the package using the following code in R:

`devtools::install_github("masonyoungblood/TwitterABM")`

You can access the documentation of the ABM function by entering `?twitter_ABM` in R.

The `example` folder includes the analysis scripts and data for the [study that this model was developed for](https://psyarxiv.com/2jksg/). The three analysis scripts are labelled in the order in which they should be run, and the output files from each script are available in respective subfolders in `data`. The only saved model output that we excluded from the repository (too large) is the full best fitting model, which can be run locally using code from `3_glmm.R`.
