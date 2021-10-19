## Data 
### Obtain data from NICHD DASH
The raw data are available in the form of sas7bdat files at https://dash.nichd.nih.gov/study/13016. The descriptive documents could be found at the same link. The variables we used in the paper including:
- 4 response variables: PxDFR(parent-reported family responsibility), CxDFR (child-reported familiy responsibility), PxPCC (parent-reported parent-child conflict) and CxPCC (child-reported parent-child relationship) at 0-, 6-, 12-, 18- and 24-month
- 1 predictor: GROUP (Control/treatment group)
- 2 adjusting covariates: AGE (age at 1st home visit) and SEX (gender of the child).
### Make working datasets
- The file containing 4 response variables `FMOD_imputed.csv`: The first column is the family's id and the last for columns are the 4 manifesto variables, the dataset is ordered by id
- The file containing the predictor 'A1C.csv': One family one row ordered by id. The id's order is the same as the one for 'FMOD_imputed.csv'
- The file containing the covariates `age_sex.csv`: Five rows for one family ordered by id. The id's order is the same as the one for 'FMOD_imputed.csv'
### Obtain working datasets
The working datasets for regenerating results in the paper could be obtained for the purpose of reviewing upon request.
