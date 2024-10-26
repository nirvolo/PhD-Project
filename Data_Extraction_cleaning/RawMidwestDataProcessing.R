# Loading the midwest yield data
il_yield = read.csv("Data_Extraction_cleaning/ExcelDataFiles/IllinoisRawYield_RmvD_RmvEmpty.csv")
in_yield = read.csv("Data_Extraction_cleaning/ExcelDataFiles/IndianaRawYield_RmvD_RmvEmpty.csv")
ia_yield = read.csv("Data_Extraction_cleaning/ExcelDataFiles/IowaRawYield_RmvD_RmvEmpty.csv")
mo_yield = read.csv("Data_Extraction_cleaning/ExcelDataFiles/MissouriRawYield_RmvD_RmvEmpty.csv")
ks_yield = read.csv("Data_Extraction_cleaning/ExcelDataFiles/KansasRawYield_RmvD_RmvEmpty.csv")
# Combining the yield data from all states
midw_yield = rbind(il_yield, in_yield, ia_yield, mo_yield, ks_yield)
# Changing the column names for the area, yield and irrigation
colnames(midw_yield)[19:24] = c("Area", "Area_CV", "Yield", "Yield_CV", "Area.Irrigated", "Area.Irrigated_CV")
# Removing , from the integers in Area and Area.Irrigated
midw_yield$Area = gsub(",","", midw_yield$Area)
midw_yield$Area.Irrigated = gsub(",","", midw_yield$Area.Irrigated)
# Converting Area and Area.Irrigated to numeric variables
midw_yield |> lapply(\(x) class(x))
midw_yield$Area = as.numeric(midw_yield$Area)
midw_yield$Area.Irrigated = as.numeric(midw_yield$Area.Irrigated)
midw_yield |> lapply(\(x) class(x))

# Filtering out the years in midw_yield that with NA yield values
nrow(midw_yield) # 11,729 years of data before removing NA
midw_yield_noNA = midw_yield[!is.na(midw_yield$Yield),]
sum(!is.na(midw_yield_noNA$Area.Irrigated))
sum(is.na(midw_yield_noNA$Area.Irrigated))
nrow(midw_yield_noNA) # 9,652 years of data after removing NA
head(midw_yield_noNA)

# Removing the county called "OTHER (COMBINED) COUNTIES"
midw_yield_noNA = midw_yield_noNA[!midw_yield_noNA$County %in% c("OTHER (COMBINED) COUNTIES"),]

# Loading the midwest temperature data where the NA values are filled
midw_temp_NAfilled = read.csv("Data_Extraction_cleaning/ExcelDataFiles/MidwTempDataNA_filled.csv") 
# Removing any rows that still have NA values for max or min temp
# NA values weren't filled if there were over 100 NAs or over 20 consecutive NAs
midw_temp_NAfilled = midw_temp_NAfilled[!is.na(midw_temp_NAfilled$mint),]
midw_temp_NAfilled = midw_temp_NAfilled[!is.na(midw_temp_NAfilled$maxt),]

# Looking at how many NA irrigation values there are for each county
summary(midw_yield_noNA$Area.Irrigated) # After removing yield NA, there are 884 values of Area.Irrigated
num_irrig_per_cnty = c()  
irrig_mn_per_cnty = c() # To store the mean of irrigation per county
unq_cntys = unique(midw_yield_noNA[,c("State","County")])
for(i in 1:nrow(unq_cntys)){
  st = unq_cntys[i,1]
  cnty = unq_cntys[i,2]
  crnt_dat = midw_yield_noNA[midw_yield_noNA$State == st & midw_yield_noNA$County == cnty, "Area.Irrigated"]
  num_irrig_per_cnty[i] = sum(is.na(crnt_dat))
  # Calculating the cnty_mean for irrigation
  irrig_cnty_mean = mean(crnt_dat, na.rm = T)
  indx_na = which(is.na(crnt_dat))
  # Replacing the NA values in the given state and county with the irrigation county mean
  midw_yield_noNA[midw_yield_noNA$State == st & midw_yield_noNA$County == cnty, 
                  "Area.Irrigated"][indx_na]  = irrig_cnty_mean
}

# Creating a new column called irrig_actualprop. This is the ratio between area.irrig and area.
midw_yield_noNA$Irrigated.actualprop  = midw_yield_noNA$Area.Irrigated/midw_yield_noNA$Area
# If there are still NA values for irrigated.actualprop, simply replacing it with 0
midw_yield_noNA[is.na(midw_yield_noNA$Irrigated.actualprop),"Irrigated.actualprop"] = 0 


# Converting the county and state column names in the yield data to lower case so it matches the temp data.
# Also removing punctation from county names
midw_yield_noNA$County = tolower(midw_yield_noNA$County)
midw_yield_noNA$State = tolower(midw_yield_noNA$State)
midw_yield_noNA$County = gsub(" |'|\\.", "", midw_yield_noNA$County)

# Calculating the average precipitation for each year using the prcp in the midwest temp data
avg_prcp = c()
for (i in 1:nrow(midw_yield_noNA)){
  prcp.st = midw_yield_noNA[i, "State"]
  prcp.yr = midw_yield_noNA[i, "Year"]
  prcp.cnty = midw_yield_noNA[i,"County"]
  prcp.set = subset(midw_temp_NAfilled, State == prcp.st & Year==prcp.yr & County == prcp.cnty)
  # The units for prcp is in millimeters (converted in python already)
  avg_prcp[i] = mean(prcp.set$prcp, na.rm = TRUE)
}
# Adding precipitation to midw_yield data
midw_yield_noNA$avgPRCP = avg_prcp
# Removing years which don't have prcp and columns that I don't need
midw_regdat = midw_yield_noNA[!is.na(midw_yield_noNA$avgPRCP),c("Year","State","County","Area","Yield",
                                                                "Area.Irrigated","Irrigated.actualprop","avgPRCP")]

# Adding County ID to midw_regdat and midw_temp_noNA. Some states have the same county names
# so I need to iterate over states as well
midw_temp_NAfilled$CountyI = NA
midw_regdat$CountyI = NA
unq_cnty_st = unique(midw_regdat[,c("State","County")])
for(i in 1:nrow(unq_cnty_st)){
  cnty = unq_cnty_st[i,"County"]
  st = unq_cnty_st[i, "State"]
  # For each combination of County and State, it gets its own county id i
  midw_regdat[midw_regdat$State == st & midw_regdat$County == cnty,"CountyI"] = i
  midw_temp_NAfilled[midw_temp_NAfilled$State == st & midw_temp_NAfilled$County == cnty,"CountyI"] = i
}

# Removing years which aren't complete (doesn't include all the days in the year)
# Creating a vector to store the indices
rmv_yr_ind = c()
for(yr in unique(midw_temp_NAfilled$Year)){
  # Iterate over each county and check how many days of temp data there is.
  # Using the state and county since there are counties with the same name
  unq_cnty_st_temp = unique(midw_temp_NAfilled[midw_temp_NAfilled$Year == yr, c("State","County")])
  for(j in 1:nrow(unq_cnty_st_temp)){
    cnty = unq_cnty_st_temp[j,"County"]
    st = unq_cnty_st_temp[j, "State"]    
    # The indices corresponding to the current year and county
    crnt_ind = as.double(which(midw_temp_NAfilled$Year == yr & midw_temp_NAfilled$County == cnty 
                               & midw_temp_NAfilled$State == st))
    # Checking if the number of days is less than 365
    if(length(crnt_ind) < 365 & length(crnt_ind) != 0){
      rmv_yr_ind = append(rmv_yr_ind, crnt_ind)
    }
    if(length(crnt_ind) > 366) print("leap year")
  }
}
# Removing the incomplete years from midw_temp_NAfilled. Also excluding the prcp variable since I don't need it for temp
midw_fd = midw_temp_NAfilled[-rmv_yr_ind,-7] # 9743 years of data
# Changing the column names of midw_fd
colnames(midw_fd) = c("State", "County", "Year","Day", "TMAX","TMIN", "CountyI")

# Checking if there is any county/year in regdat and not in fundat
rm_indx_regdat = c()
for(yr in unique(midw_regdat$Year)){
  # Iterate over each county and check how many days of temp data there is.
  # Using the county ID since there are counties with the same name
  for(id in unique(midw_regdat[midw_regdat$Year == yr,"CountyI"])){
    # The indices corresponding to the current year and county
    crnt_ind = as.double(which(midw_fd$Year == yr & midw_fd$CountyI == id))
    # If crnt_ind is empty, I want to remove that county and year from regdat
    if(length(crnt_ind) == 0){
      # Storing the indices I need to remove from regdat so that it matches with fundat
      rm_indx_regdat = append(rm_indx_regdat, which(midw_regdat$Year == yr & midw_regdat$CountyI == id))
    }
  }
}
# Subset of regdat that matches all combinations of year and county id in fundat
midw_regdat = midw_regdat[-rm_indx_regdat,]

################################# Running the preprocessing function on the midwest data ##################################
source("common_utils/preprocessing_function.R") 
# The midwest boundary
state = map_data("state")
midwest = state[state$region %in% c("kansas", "missouri", "iowa", "illinois", "indiana"),]
kansas <- state[state$region %in% c("kansas"),]
missouri <- state[state$region %in% c("missouri"),]
iowa <- state[state$region %in% c("iowa"),]
illinois <- state[state$region %in% c("illinois"),]
indiana <- state[state$region %in% c("indiana"),]
midw_bnd = list("kns" = kansas, "iln" = illinois, "ind" = indiana, "iowa" = iowa, "miss" = missouri)

# Saving the fd and nonfd data before preprocessing so that I can use it for Kansas. I need to run preprocessing for Kansas
# again since the centroid is different
# save(midw_fd, file = "RawMidwestFundatBeforePreProcFxn_incl_irrig.RData")
# save(midw_regdat, file = "RawMidwestRegdatBeforePreProcFxn_incl_irrig.RData")

# Preprocessing the data including irrigation
midw_raw_preprocess = preprocessing_nosmooth(fun_dat = midw_fd, non_fd = midw_regdat, states = c("Kansas", "Illinois", "Indiana", "Iowa", "Missouri"), 
                                             vars = c("Yield", "Area", "Irrigated.actualprop","avgPRCP", "CountyI", "Year", "State"), st_bnd = midw_bnd)
midw_fd_new = midw_raw_preprocess$fd_dat
midw_regdat_new = midw_raw_preprocess$non_fd

# Adding the Id's for both datasets since it's needed for training and test
midw_regdat_new[,"Id"] = seq(1, nrow(midw_regdat_new))
midw_fd_new[,"Id"] = NA
# Iterating over all the rows since I still need to add county Id's for counties
# that aren't in regdat. Look at the issue below. The way I coded this is fine since
# regdat and fd are in the same order so setting the obs Id equal to i will work. It'll 
# go up to 8405 (num obs in regdat) and then the rest of the Id's will be 8406-9545 which
# don't exist in regdat
for(i in 1:(nrow(midw_fd_new)/365)){
  midw_fd_new[((i-1)*365 + 1):(365*i),"Id"] = i
}

# Changing column names in fd from Day to DateI so that it works for the PLFAM code 
colnames(midw_fd_new) = c("State", "County", "Year",  "DateI",  "TMAX", "TMIN", "CountyI", "county", "long", "lat",    
                          "Id")
# save(midw_fd_new, file = "RawMidwestFundatPreProc_incl_irrig.RData")
# save(midw_regdat_new, file = "RawMidwestRegdatPreProc_incl_irrig.RData")





