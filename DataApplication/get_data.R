setwd("/Users/shushuz/Dropbox (University of Michigan)/ES_HighD/Empirical_studies/DataApplication")
library(haven)
library(quantreg)
library(xtable)
# For data wrangling
#dplyr::select, mutate, select, recode
library(dplyr)
library(survey)
library(conquer)
library(glmnet)
library(janitor)

yearletters = "P"
years = "2017-2018"
yearidx = 1

## 1. Demographic data
data<- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_DEMO.XPT",
                            years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data$wt = data$WTMECPRP
data$wt8yr= 1/4 * data$WTMECPRP
data$psu = data$SDMVPSU
data$strata = data$SDMVSTRA
data$sex = data$RIAGENDR
data$age = data$RIDAGEYR
race3=data$RIDRETH3
data$race = ifelse(race3 %in% c(7,NA), NA,
                   case_when(race3 %in% c(1,2)~"M",
                             race3%in% c(3)~ "W",
                             race3%in% c(4)~"B",
                             race3 %in% c(6)~ "A"))
#data$race= ifelse(race3c=="W",1, ifelse(race3c=="B", 2, ifelse(race3c=="M",3, ifelse(race3c== "A",4,NA))))
data$wave=rep(yearletters[yearidx], nrow(data) )     
data$pir = data$INDFMPIR
data$edu= data$DMDEDUC2
data$edu[data$edu %in% c(NA,7,9)] = NA

Demo <- data %>% dplyr::select(race,id, sex, age, pir, edu)
#Demo <- data %>% dplyr::select(race,id, wt, wt8yr, psu, strata, sex, age, pir, edu)

## 2. Hospital Utilization & Access to Care
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_HUQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data$access_care = ifelse(data$HUQ030 %in% c(NA,7,9),NA,
                          case_when(data$HUQ030==2 ~ 1,#no place to go for health care
                                    data$HUQ030==1 ~ 2, # one place
                                    data$HUQ030==3 ~ 3)) #more than one
data$HUQ051[data$HUQ051 %in% c(NA,77,99)] = NA
data$times_care = data$HUQ051
data$HUD062[data$HUD062 %in% c(NA,77,99)] = NA
data$time_since_last_visit = data$HUD062
data$mental = ifelse(data$HUQ090 %in% c(NA,7,9),NA,
                     case_when(data$HUQ090==1 ~ 1,
                               data$HUQ090==2 ~ 2))
data$overnight = ifelse(data$HUQ071 %in% c(NA,7,9),NA,
                        case_when(data$HUQ071==1 ~ 1,
                                  data$HUQ071==2 ~ 2))
HUQ <- data %>% dplyr::select(id,access_care,mental,overnight)


## 3. Insurance
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_HIQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data$insurance = ifelse(data$HIQ011 %in% c(NA,7,9),NA,
                        case_when(data$HIQ011==2 ~ 1,
                                  data$HIQ011==1 ~ 2))
data$ever_no_insurance = ifelse(data$HIQ210 %in% c(NA,7,9),NA,
                                case_when(data$HIQ210==2 ~ 2,
                                          data$HIQ210==1 ~ 1))
HIQ <- data %>% dplyr::select(id,insurance,ever_no_insurance)


#4. physical activity
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_PAQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
vigworknumdays = data$PAQ610
vigworknumdays[vigworknumdays %in% c(77,99,NA)] = NA
vigworkminday = data$PAD615
vigworkminday[vigworkminday %in% c(7777,9999,NA)] = NA
vigwork = data$PAQ605
vigwork[vigwork %in% c(7,9,NA)] = NA
vigwork[vigwork == 2] = 0
vigwork[vigwork == 1&!is.na(vigwork)] = (vigworknumdays*vigworkminday)[vigwork == 1&!is.na(vigwork)]
vigwork_wt = vigwork * 8

modworknumdays = data$PAQ625
modworknumdays[modworknumdays %in% c(77,99,NA)] = NA
modworkminday = data$PAD630
modworkminday[modworkminday %in% c(7777,9999,NA)] = NA
modwork = data$PAQ620
modwork[modwork %in% c(7,9,NA)] = NA
modwork[modwork == 2] = 0
modwork[modwork == 1&!is.na(modwork)] = (modworknumdays*modworkminday)[modwork == 1&!is.na(modwork)]
modwork_wt = modwork * 4

walknumdays = data$PAQ640
walknumdays[walknumdays %in% c(77,99,NA)] = NA
walkminday = data$PAD645
walkminday[walkminday %in% c(7777,9999,NA)] = NA
walk = data$PAQ635
walk[walk %in% c(7,9,NA)] = NA
walk[walk==2] = 0
walk[walk == 1&!is.na(walk)] = (walknumdays*walkminday)[walk == 1&!is.na(walk)]
walk_wt = walk * 4

vigrecnumdays = data$PAQ655
vigrecnumdays[vigrecnumdays %in% c(77,99,NA)] = NA
vigrecminday = data$PAD660
vigrecminday[vigrecminday %in% c(7777,9999,NA)] = NA
vigrec = data$PAQ650
vigrec[vigrec %in% c(7,9,NA)] = NA
vigrec[vigrec == 2] = 0
vigrec[vigrec == 1&!is.na(vigrec)] = (vigrecnumdays*vigrecminday)[vigrec==1&!is.na(vigrec)]
vigrec_wt = vigrec*8

modrecnumdays = data$PAQ670
modrecnumdays[modrecnumdays %in% c(77,99,NA)] = NA
modrecminday = data$PAD675
modrecminday[modrecminday %in% c(7777,9999,NA)] = NA
modrec = data$PAQ665
modrec[modrec %in% c(7,9,NA)] = NA
modrec[modrec == 2] = 0
modrec[modrec == 1&!is.na(modrec)] = (modrecnumdays*modrecminday)[modrec == 1&!is.na(modrec)]
modrec_wt = modrec * 4
data$vigwork_wt = vigwork_wt
data$modwork_wt = modwork_wt
data$walk_wt = walk_wt
data$vigrec_wt = vigrec_wt
data$modrec_wt = modrec_wt
#data$pa = vigwork_wt + modwork_wt + walk_wt + vigrec_wt + modrec_wt
PAQ = data %>% dplyr::select(id,vigwork_wt,modwork_wt, walk_wt, vigrec_wt, modrec_wt)


## 5. Medical Conditions
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_MCQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data$asthma = data$MCQ010
data$asthma[data$asthma %in% c(7,9,NA)] = NA
data$asthma_still = data$MCQ035
data$asthma_still[data$asthma_still %in% c(7,9,NA)] = NA
data$asmthma_past_year = data$MCQ040
data$asmthma_past_year[data$asmthma_past_year %in% c(7,9,NA)] = NA
data$asmthma_ER = data$MCQ050
data$asmthma_ER[data$asmthma_ER %in% c(7,9,NA)] = NA
data$fever = data$AGQ030
data$fever[data$fever %in% c(7,9,NA)] = NA
data$treat = data$MCQ053
data$treat[data$treat %in% c(7,9,NA)] = NA
data$overweight = data$MCQ080
data$overweight[data$overweight %in% c(7,9,NA)] = NA
data$blood_trans = data$MCQ092
data$blood_trans[data$blood_trans %in% c(7,9,NA)] = NA
data$trans_year = data$MCD093
data$trans_year[data$trans_year %in% c(7,9,NA)] = NA
data$period = data$MCQ149
data$period[data$period %in% c(7,9,NA)] = NA
data$age_period = data$MCQ151
data$age_period[data$age_period %in% c(77,99,NA)] = NA
data$age_menarche = data$RHD018
data$age_menarche[data$age_menarche %in% c(7777,9999,NA)] = NA
data$arthritis = data$MCQ160A
data$arthritis[data$arthritis %in% c(7,9,NA)] = NA
data$arthritis_type = data$MCQ195
data$arthritis_type[data$arthritis_type %in% c(7,9,NA)] = NA
data$heart = data$MCQ160B
data$heart[data$heart %in% c(7,9,NA)] = NA
data$age_heart = data$MCD180B
data$age_heart[data$age_heart %in% c(77777,99999,NA)] = NA
data$coronary = data$MCQ160C
data$coronary[data$coronary %in% c(7,9,NA)] = NA
data$age_coronary = data$MCD180C
data$age_coronary[data$age_coronary %in% c(77777,99999,NA)] = NA
data$angina = data$MCQ160D
data$angina[data$angina %in% c(7,9,NA)] = NA
data$age_angina = data$MCD180D
data$age_angina[data$age_angina %in% c(77777,99999,NA)] = NA
data$heart_attach = data$MCQ160E
data$heart_attach[data$heart_attach %in% c(7,9,NA)] = NA
data$age_heart_attach = data$MCD180E
data$age_heart_attach[data$age_heart_attach %in% c(77777,99999,NA)] = NA
data$stroke = data$MCQ160F
data$stroke[data$stroke %in% c(7,9,NA)] = NA
data$age_stroke = data$MCD180F
data$age_stroke[data$age_stroke %in% c(77777,99999,NA)] = NA
data$thyroid = data$MCQ160M
data$thyroid[data$thyroid %in% c(7,9,NA)] = NA
data$thyroid_still = data$MCQ170M
data$thyroid_still[data$thyroid_still %in% c(7,9,NA)] = NA
data$age_thyroid = data$MCD180M
data$age_thyroid[data$age_thyroid %in% c(77777,99999,NA)] = NA
data$COPD = data$MCQ160P
data$COPD[data$COPD %in% c(7,9,NA)] = NA
data$liver = data$MCQ160L
data$liver[data$liver %in% c(7,9,NA)] = NA
data$liver_still = data$MCQ170L
data$liver_still[data$liver_still %in% c(7,9,NA)] = NA
data$age_liver = data$MCD180L
data$age_liver[data$age_liver %in% c(77777,99999,NA)] = NA
data$liver_cond = data$MCQ500
data$liver_cond[data$liver_cond %in% c(7,9,NA)] = NA
data$fat_liver = data$MCQ510A
data$fat_liver[data$fat_liver %in% c(77,99,NA)] = NA
data$fibrosis_liver = data$MCQ510B
data$cirrhosis_liver = data$MCQ510C
data$Viral_liver = data$MCQ510D
data$other_liver = data$MCQ510F
data$pain = data$MCQ520
data$pain[data$pain %in% c(7,9,NA)] = NA
data$pain_unconf = data$MCQ530
data$pain_unconf[data$pain_unconf %in% c(7,9,NA)] = NA
data$pain_DR = data$MCQ540
data$pain_DR[data$pain_DR %in% c(7,9,NA)] = NA
data$gallstones = data$MCQ550
data$gallstones[data$gallstones %in% c(7,9,NA)] = NA
data$gallstones_surgery = data$MCQ560
data$gallstones_surgery[data$gallstones_surgery %in% c(7,9,NA)] = NA
data$age_gallbladder = data$MCQ570
data$age_gallbladder[data$age_gallbladder %in% c(77777,99999,NA)] = NA
data$cancer = data$MCQ220
data$cancer[data$cancer %in% c(77,99,NA)] = NA
data$cancer1 = data$MCQ230A
data$cancer1[data$cancer1 %in% c(77,99,NA)] = NA
data$cancer2 = data$MCQ230B
data$cancer2[data$cancer2 %in% c(77,99,NA)] = NA
data$cancer3 = data$MCQ230C
data$cancer3[data$cancer3 %in% c(77,99,NA)] = NA
data$cancer4 = data$MCQ230D
data$cancer4[data$cancer4 %in% c(77,99,NA)] = NA
data$asthma_relative = data$MCQ300B
data$asthma_relative[data$asthma_relative %in% c(7,9,NA)] = NA
data$diabete_relative = data$MCQ300C
data$diabete_relative[data$diabete_relative %in% c(7,9,NA)] = NA
data$heart_relative = data$MCQ300A
data$heart_relative[data$heart_relative %in% c(7,9,NA)] = NA
data$lose_weight = data$MCQ366A
data$lose_weight[data$lose_weight %in% c(7,9,NA)] = NA
data$excercise = data$MCQ366B
data$excercise[data$excercise %in% c(7,9,NA)] = NA
data$add_salt = data$MCQ366C
data$add_salt[data$add_salt %in% c(7,9,NA)] = NA
data$reduce_fat = data$MCQ366D
data$reduce_fat[data$reduce_fat %in% c(7,9,NA)] = NA
data$lose_weight_now = data$MCQ371A
data$lose_weight_now[data$lose_weight_now %in% c(7,9,NA)] = NA
data$excercise_now = data$MCQ371B
data$excercise_now[data$excercise_now %in% c(7,9,NA)] = NA
data$reduce_salt_now = data$MCQ371C
data$reduce_salt_now[data$reduce_salt_now %in% c(7,9,NA)] = NA
data$reduce_fat_now = data$MCQ371D
data$reduce_fat_now[data$reduce_fat_now %in% c(7,9,NA)] = NA
data$metal = data$OSQ230
data$metal[data$metal %in% c(7,9,NA)] = NA

#data[,-which(names(data)) %*% in ]
Medical = data[,64:ncol(data)]

## 6. Mental Health
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_DPQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data[data %in% c(7,9,NA)] = NA
Mental <- data[,-1]

## 7. Occupation
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_OCQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data$work_type = data$OCD150
data$work_type[data$work_type %in% c(7,9,NA)] = NA
data$work_hours = data$OCQ180
data$work_hours[data$work_hours %in% c(77777,99999,NA)] = NA
data$long_hours = data$OCQ210
data$long_hours[data$long_hours %in% c(7,9,NA)] = NA
data$work_schedule = data$OCQ670
data$work_schedule[data$work_schedule %in% c(7,9,NA)] = NA
data$not_work = data$OCD383
data$not_work[data$not_work %in% c(77,99,NA)] = NA
Work= data %>% dplyr::select(id,work_type,work_hours,long_hours,work_schedule,not_work)

## 8. Asprin Use
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_RXQASA.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
Asprin <- data[,-1]

# ## 9. Reproductive Health
# data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_RHQ.XPT",
#                              years[yearidx],yearletters[yearidx])))
# data$id = data$SEQN
# data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
# RHQ <- data[,-1]

## 9. Weight History
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_WHQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
WHQ <- data[,-1]

## 10. Volatile Toxicant (P_VTQ)
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_VTQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
VTQ <- data[,-1]

## 11. Smoking - Household Smokers (P_SMQFAM)
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_SMQFAM.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
SMQFAM <- data[,-1]

## 12. Smoking - Cigarette Use (P_SMQ)
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_SMQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
SMQ <- data[,-1]

## 13. Smoking - Recent Tobacco Use (P_SMQRTU)
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_SMQRTU.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
SMQRTU <- data[,-1]

## 14. Smoking - Secondhand Smoke Exposure (P_SMQSHS)
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_SMQSHS.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
SMQSHS <- data[,-1]

## 15. Osteoporosis (P_OSQ)
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_OSQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
OSQ <- data[,-1]

## 16. Pesticide Use (P_PUQMEC)
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_PUQMEC.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
PUQMEC <- data[,-1]

## 17. Food Security (FSQ)
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_FSQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
FSQ <- data[,-1]

## 18. Body Measures (BMX)
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_BMX.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = data %>% subset(select=-c(BMDSTATS))
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
BMX <- data[,-1]

## 19. Income (P_INQ)
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_INQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
INQ <- data[,-1]

## 20. Pesticide Use (P_PUQMEC)
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_PUQMEC.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
PUQMEC <- data[,-1]

## 21. Diet Behavior & Nutrition (P_DBQ)
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_DBQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
DBQ <- data[,-1]

## 22. Current Health Status (P_HSQ)
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_HSQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
HSQ <- data[,-1]

## 23. Blood Pressure & Cholesterol (P_BPQ)
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_BPQ.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
BPQ <- data[,-1]

## 24. Prescription Medications (P_RXQ_RX)
#data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_RXQ_RX.XPT",
#                             years[yearidx],yearletters[yearidx])))
# data$id = data$SEQN
# data = apply(data, c(1,2), FUN = function(x) ifelse(x %in% c(7,9,77,99,777,999,7777,9999,77777,99999,NA), NA, x))
# RXQ_RX <- data %>% dplyr::select(id, RXDUSE)

## Response: Cotinine and Hydroxycotinine
data <- read_xpt(url(sprintf("https://wwwn.cdc.gov/Nchs/Nhanes/%s/%s_COT.XPT",
                             years[yearidx],yearletters[yearidx])))
data$id = data$SEQN
data$cotinine = data$LBXCOT
LBXCOT = data %>% dplyr::select(id,cotinine)

###merge tables
tablesToMerge = c("Demo","HUQ", "HIQ", "PAQ","Medical","Work",
                  "WHQ","SMQFAM","SMQ","SMQRTU","SMQSHS",
                  "LBXCOT","FSQ","BMX","INQ","PUQMEC","DBQ","HSQ","BPQ",
                  "Asprin","VTQ","Mental")
firstTable = tablesToMerge[1]
merged = eval(parse(text = firstTable))
for (i in 2:length(tablesToMerge)){
  nextTable = tablesToMerge[i]
  merged = merge(merged, eval(parse(text = nextTable)), by.x = "id", by.y = "id")
}

#write.csv(merged,sprintf("NHANES_Cotinine_%s.csv",yearletters[yearidx]))
write.csv(merged,sprintf("NHANES_Cotinine_%s_new.csv",yearletters[yearidx]),row.names = FALSE)
# ES Regression
#y = total_w$cotinine
#x = total_w[c("race","sex","age","pir","edu","access_care","mental","overnight",
#              "insurance","ever_no_insurance","pa")]

setwd("/Users/shushuz/Dropbox (University of Michigan)/ES_HighD/Empirical_studies/DataApplication")
yearletters = "P"
years = "2017-2018"
yearidx = 1
merged = read.csv(sprintf("NHANES_Cotinine_%s_new.csv",yearletters[yearidx]))

#mapply(function(x) length(unique(x)), merged)
filter = Filter(function(x) length(unique(x))>1, merged) #remove constant columns
filter = filter[!is.na(filter$cotinine),]  # drop those cotinine NA
filter$race = relevel(factor(filter$race),"W")
filter = subset(filter,select=-c(id))
# remove those variables with a lot of missing values
##colnames(total_w)[sapply(total_w, class)=="numeric" & colMeans(is.na(total_w)) > 0.9]


#remove_cols = filter[colMeans(is.na(filter)) < 0.3]


factorize = function(x){
  if(length(unique(x))<=7){
    return(factor(x,exclude = NULL))
  } else{
    return(x)
  }
}
data_clean = data.frame(lapply(filter,factorize))
# remove numeric variables with >10% missingness
data_clean = data_clean[colMeans(is.na(data_clean)) < 0.1]
# remove categorical variables with >50% missingness
data_clean = data_clean[colMeans(data_clean==NA) < 0.1]
# remove numeric variables with >30% missingness
#num_col_remove = colnames(data_clean)[sapply(data_clean, class)=="numeric" & colMeans(is.na(data_clean)) > 0.1]
#data_clean = data_clean[,!(names(data_clean) %in% num_col_remove)]
#remove numeric rows with missing values
remove_rows = data_clean[rowMeans(is.na(data_clean[,sapply(data_clean, class)=="numeric"]))==0,]

design_matrix = model.matrix(~.,data=remove_rows)
design_matrix = remove_constant(design_matrix, na.rm = FALSE, quiet = TRUE)
dim(design_matrix)
write.csv(design_matrix,"design_matrix_new.csv",row.names = FALSE)

#############################
mapply(function(x) length(unique(x)), filter)
factorize = function(x){
  if(length(unique(x))<=7){
    return(factor(x,exclude = NULL))
  } else{
      return(x)
    }
}
total_w = data.frame(lapply(filter,factorize))
# remove those numeric variables with a lot of missing values
##colnames(total_w)[sapply(total_w, class)=="numeric" & colMeans(is.na(total_w)) > 0.9]
remove_cols = total_w[colMeans(is.na(total_w)) < 0.1]
#remove rows with missing values
remove_rows = remove_cols[rowMeans(is.na(remove_cols))==0,]
# remove sampling variables
data_clean = subset(remove_rows,select=-c(id))#,wt,wt8yr,psu,strata
dim(data_clean)
design_matrix = model.matrix(~.,data=data_clean)
design_matrix = remove_constant(design_matrix, na.rm = FALSE, quiet = TRUE)
dim(design_matrix)
write.csv(design_matrix,"design_matrix_new.csv",row.names = FALSE)
