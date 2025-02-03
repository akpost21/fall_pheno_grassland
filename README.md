# Modeling fall (autumn) phenology in grasslands

Code used for analyses in this publication: 

*Post, A.K. & Richardson, A.D. (Accepted). Predicting End-of-Season Timing Across Diverse North American Grasslands. Oecologia.*


**For code to work, use my forked version of the "phenor" R package, available here:**
[https://github.com/akpost21/phenor](https://github.com/akpost21/phenor)


## Contents:

#### Code folder
- New_fall_models_final.R: Code for new fall phenology models
- model_fitting.R: Code to fit new models to PhenoCam data

#### Files folder
- AllSites_fall_1day_90_trans50_noOut_v4.rds: Full dataset of PhenoCam fall transition dates ("AllSites") used in paper
- fall_parameter_ranges_final.csv: contains parameter ranges used to fit models

*We used the same steps to prepare and process the fall transition date (EOS) dataset as previously used for the spring transition date (SOS) dataset [(Post et al., 2022)](https://www.sciencedirect.com/science/article/pii/S0168192322003914). Please see the ["spring_pheno_grassland" repository](https://github.com/akpost21/Spring_pheno_grassland) for additional details.
