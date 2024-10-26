# PhD projects code
[Dissertation](https://github.com/nirvolo/PhD-Project/blob/main/TopicsInPredictionForSpatiallyVaryingFunctionalData.pdf)

## Introduction
This repository contains the code for the Dissertation titled "Topics in Prediction for Spatially Varying Functional Data" by Nir Voloshin. All the necessary data, models and model fitting are included.

### Dissertation abstract
Functional data are inherently infinite-dimensional because they represent continuous functions. Since the GSVFM suffers from the curse of dimensionality, functional models can not be estimated directly. To address both the infinite-dimensionality and the spatial varying components of the data, a novel two-step procedure is introduced. The first step of the procedure is to reduce the dimension of the SVFM through the method of Functional Principal Components Analysis (FPCA). This reduces the GSVFM to a Generalized Spatial Varying Coefficient Model (GSVCM) which is the second step in the procedure. The GSVCM considers the spatial locations in the data. The proposed two-step procedure is able to capture location-specific effects that previous functional regression models can't. 

This research is motivated by a crop-yield prediction application in agriculture. The agriculture data is collected at the county-level from five Midwest states, Kansas, Iowa, Illinois, Indiana and Missouri. For each county, we observe daily minimum and maximum temperature time series data. The temperature time series data can be viewed as functions, where the temperature is indexed by the day. Since the temperature data varies across the Midwest counties, this represents the multivariate spatially varying functional data. The precipitation, irrigated land and crop-yield are collected at the county level. The goal is to apply the GSVFM to predict the spatially varying crop yield through the scalar predictor variables and the multivariate spatially varying functional data. Existing functional models are used to compare performance with the GSVFM. 
    
The dissertation consists of two projects that use the novel two-step procedure to estimate the GSVFM and the SVFQM. The first project aims at predicting the conditional mean and the second project extends the GSVFM model to the SVFQM that predicts the conditional quantile. This research addresses the current gap in functional models that do not consider the spatial component.

## Directory structure

### GSVFM
This folder contains all the code to run the GSVFM for both simulations and real data. The necessary functions to run the GSVFM can be found in "util_funcs" and "app_funcs". 

The GSVFM is compared to the PLFAM and FLM which have their own separate code under the folder "util_funcs". The "util_funcs" consist of code that is used throughout all the code files. There is one main file that is used for the simulations which can be found under the folder "sim_funcs". The GSVCM code developed by Kim and Wang (2021), can be found under the folder "GSVCM_Kim_Wang". The model fitting the Gaussian case including simulations is performed separately in its own file titled accordingly. Since the GSVFM is a generalized model, it can also be used for a Poisson response variable. The corresponding simulation results can be found in "PoissonSimulations.R".

The results for the bootstrap hypothesis test can be found in the files titled "KansasBootstrapTest.R" and "MidwestBootstrapTest.R". 

### SVFQM
The SVFQM folder contains all the necessary code to run the SVFQM for both simulations and real data. The main functions can be found in the folders "sim_funcs","app_funcs", and "util_funcs". The rqPen package is used to fit the SVFQM but there are some changes that are made to the original functions. These updated functions can be found under the "rqPen" folder. The files "KansasQuantRegRes.R", "MidwestQuantRegResults.R", and "SimTrialsQuantReg.R" contain the results for the real data and simulations.

### common_util
This folder contains files that are used throughout both projects. One of the main files which is necessary for both projects is the "preprocessing_function.R". This contains the code that processes the functional and non-functional data before model fitting. The files "KansasMaps.Rmd" and "MidwestMaps.Rmd" contain the code for obtaining the county centroids and construction of the triangulation. The file "AgricultureDataEDA.Rmd" includes some exploratory data analysis that is performed for the agriculture data.

### Data_Extraction
This folder contains the code that is used for extracting the temperature data from the servers and the preprocessing of this data and the agriculture data.












