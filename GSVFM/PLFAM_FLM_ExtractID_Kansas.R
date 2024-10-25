##################################################################################
# This file determines the train and test indices for PLFAM and FLM model fitting
##################################################################################

#### load in data file ####
load("Data_Extraction_cleaning/RdataFiles/NewKansasFromMidw_PreProcRegdat.RData")
# Defining the nonfd data objects
nonfd = kns_nonfd_new
n_iters = 9 # number of repetitions
n_folds <- 5 # number of folds in validations
test_inds = matrix(NA, nr = 300, nc = n_folds*n_iters)

# define expression for tryCatch
oneiteration_test_inds_nosmooth <- expression({
  set.seed(12876+i*123)
  # Percentage of data to use for training
  ratio = 1 - (1/n_folds)
  
  # For each fold I am doing the following
  for (j in (1:n_folds)){
    cat("Fold", j, "\n")
    cat('------------', "\n")
    # Indices to index the nonfd matrices
    nfd_tst_ind = 0
    # Vector to store cnty ids that have less than a certain number of years
    excl_cnty_ids = c()
    for(cnty_id in unique(nonfd$CountyI)){  
      n_obs = nrow(nonfd[nonfd$CountyI == cnty_id,])
      # Not including counties that have less than 8 years of data
      if(n_obs < 5){ 
        cat("yes",cnty_id,"\n") 
        excl_cnty_ids = append(excl_cnty_ids, cnty_id)
        next
      }
      # Number of training and testing observations for this county
      n_train = floor(ratio*n_obs)
      n_test = n_obs - n_train
      # The train and test indices for the current county
      # In each fold we will get a different set of indices for training with slight
      # probability for repetition.
      trn_inds = sample(1:n_obs, n_train)
      tst_inds = (1:n_obs)[-trn_inds]
      # The nonfd data for the current county and window
      nonfd_cnty = nonfd[nonfd$CountyI == cnty_id,]
      # Adding the test Id's for the current county, fold and iteration
      test_inds[(nfd_tst_ind + 1):(nfd_tst_ind + n_test),(i-1)*n_folds + j] = nonfd_cnty[tst_inds,"Id"]
      # Updating the indices for the nonfd matrices
      nfd_tst_ind = nfd_tst_ind + n_test
    } # inner loop
  } # outer loop
})

for (i in 1:n_iters){
  eval(eval(oneiteration_test_inds_nosmooth))
}



