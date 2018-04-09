library(tidyverse)
library(lubridate)
library(here)
# Wed Nov 15 09:52:14 2017 ------------------------------

## Notes: 1. Only have data from 11 pentads
#         2. There are a number of records that have no associated 
#            track ID (point records) - confirm they are ad hoc records  
#         3. Import spatial data and add covariates (survey duration, distance)
#         4. Import temperature data and add as detection covariate
#         5. Given the different searching abilities, observer ID should be added as a detection covariate (factor)

# Import specimen data and clean ------------------------------------------
data <- read_csv(here("data input","scorpions_all.csv"))

scorp <- data %>% 
  rename_all(tolower) %>% 
  rename_all(make.names) %>%                            # Cool function (removes spaces)
  rename_all(funs(str_replace_all(.,"\\.","_"))) %>% 
  mutate(datetime = ymd_hms(str_c(track_date,track_time))) %>%
  mutate(spp = paste(genus, species, sep = " ")) %>% 
  select (pentad,datetime, family,spp, no_spec,observer) %>% 
  filter(!is.na(datetime)) %>% 
  group_by(pentad) %>% 
  mutate(transect = as.numeric(as.factor(datetime))) %>% 
  ungroup()
  
scorp

# Import temperature data and assign --------------------------------------
airtemp <- read_csv(here("data input","iButton_data.csv"))

temp1 <- airtemp %>% 
  filter(row_number() %in% 1:2011) %>% 
  mutate(datetime = mdy_hm(datetime))

temp2 <- airtemp %>% 
  filter(row_number() %in% 2012:n()) %>% 
  mutate(datetime = dmy_hm(datetime))

airtemp <- rbind(temp1,temp2)
rm(temp1);rm(temp2)

## Match closest data point
scorp$temperature <- airtemp$temp[sapply(scorp$datetime, function(z) which.min(abs(airtemp$datetime - z)))]

# Transect summaries ------------------------------------------------------
rep_df <- scorp %>% group_by(pentad) %>% 
          summarise(n_transects = length(unique((transect))))

n_reps <- rep_df %>% select(n_transects) %>% max()
rep_ref <- c(1:n_reps)

n_spp <- scorp %>% distinct(spp) %>% nrow()
spp_ref <- scorp %>% arrange(spp) %>% distinct(spp) %>% arrange() %>% pull()

n_pen <- scorp %>% distinct(pentad) %>% nrow()
pen_ref <- scorp %>% arrange(pentad) %>% distinct(pentad) %>% pull() 

airtemp <- scorp %>% 
  mutate(transid = paste0(pentad,"_",transect)) %>% 
  select(transid,temperature) %>%
  group_by(transid) %>% 
  summarise(meantemp = mean(temperature)) %>% 
  mutate(pentad = substr(transid,1,9))

observer <- scorp %>% 
  mutate(transid = paste0(pentad,"_",transect)) %>% 
  select(transid,observer) %>%
  mutate(observer = fct_recode(observer, "1" = "LP", "2" = "RC")) %>% 
  group_by(transid) %>% 
  filter(row_number() == 1) %>% 
  mutate(pentad = substr(transid,1,9))

# Create and populate arrays ----------------------------------------------
airtemp_arr <- array(NA,c(n_pen,n_reps))
dimnames(airtemp_arr) <- list(pen_ref, rep_ref)

for(p in 1:n_pen){
  ex_df1 <- filter(airtemp, pentad == pen_ref[p])
  airtemp_arr[p,1:nrow(ex_df1)] <- ex_df1$meantemp
}

observer_arr <- array(NA,c(n_pen,n_reps))
dimnames(observer_arr) <- list(pen_ref, rep_ref)

for(p in 1:n_pen){
  ex_df1 <- filter(observer, pentad == pen_ref[p])
  observer_arr[p,1:nrow(ex_df1)] <- ex_df1$observer
}

Y <- array(NA,c(n_pen, n_reps,n_spp))
dimnames(Y) <- list(pen_ref, rep_ref,spp_ref)
dimnames(Y)

for(i in c(1:n_spp)){
   ex_df1 <- scorp %>% filter(spp == spp_ref[i])
   for(p in 1:n_pen){
    ex_df2 <- filter(ex_df1, pentad == pen_ref[p])
    if(nrow(ex_df2) >= 1){
      trans <- ex_df2 %>% select(transect) %>% pull()
      Y[p,trans,i] <- 1
    }
  }
}


# Arrays check ------------------------------------------------------------

## Quick summary
for(i in c(1:n_spp)){
  ex_df1 <- scorp %>% filter(spp == spp_ref[i])
  print(nrow(ex_df1))
  print(spp_ref[i])
  
}

## Check  
Y[,,"Opistophthalmus sp."]
scorp %>% filter(spp == "Opistophthalmus sp.")

Y[,,"Opistophthalmus laticauda"]
scorp %>% filter(spp == "Opistophthalmus laticauda")

# Import, sort & scale environmental covariates ---------------------------
ndvi <- read_csv(here("data input","modis_ndvi.csv")) %>% 
  rename_all(funs(str_replace(., "July", "Jul"))) %>% 
  rename_all(funs(str_replace(.,"June","Jun")))

pencovs_all <- read_csv(here("data input","pentadCovs.csv")) %>% 
  rename_all(tolower) %>% 
  select(pentad,map,map_ctn,elev,elev_range) %>% 
  left_join(select(ndvi, pentad, mean_all),by = "pentad") %>% 
  rename(ndvi = mean_all)
pencovs_all

pencovs <- pencovs_all %>% filter(pentad %in% pen_ref) %>% arrange(pentad)
pencovs_sc <- pencovs %>% mutate_at(vars(map:ndvi),scale) %>% mutate_at(vars(map:ndvi),as.numeric)
airtemp_arr_sc <- scale(airtemp_arr)

rm(list = ls()[!ls() %in% c("data","scorp","Y","n_pen","n_reps","n_spp","airtemp_arr","airtemp_arr_sc",
                            "pen_ref", "rep_ref","spp_ref","rep_df","pencovs",
                            "pencovs_sc","pencovs_all","observer","observer_arr")])

save.image(here("data output","Arrays.RData"))
