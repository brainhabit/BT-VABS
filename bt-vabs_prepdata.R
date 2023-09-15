###########################################################################################################################################
#
# Script to prepare data set for twin analysis of the BT-VABS paper investigating
# etiological structure underlying different domains of the Vineland instrument.
#
# G. Bussu, 28 November 2022  [adapted from a script tutorial provided by M. Taylor (Feb 2022)]
#
####################################
# The script needs as input adaptive functioning and demographic scores from the entire BATSS sample in the long format.
#
# The script gives as output the selected data set based on inclusion criteria, in wide (twin_data.csv) and long format (phenotypic_data.csv).
#
###########################################################################################################################################


rm(list=ls())

### set a working directory (otherwise create a new R project):

setwd('C:/Users/giobu365/Documents/BT-VABS/data')
list.files() 

### import and check data:

data_vabs <- readxl::read_xlsx('Vineland_5m.xlsx')
names(data_vabs);dim(data_vabs)
names(data_vabs)[names(data_vabs)==" Namn/ID"]<-"id"
names(data_vabs)[names(data_vabs)=="Twin pair no."]<-"twinpair"
names(data_vabs)[names(data_vabs)=="Ålder"]<-"age"
empty_col<-colSums(is.na(data_vabs)) == nrow(data_vabs)
data_vabs<-data_vabs[,!empty_col]

data_demo <- readxl::read_xlsx('C:/Users/giobu365/Documents/BT-MSEL/data_batss_msel_vabs/Demographics_5m.xlsx')
names(data_demo);dim(data_demo)

data_exclusion<-readxl::read_xlsx('C:/Users/giobu365/Documents/BT-MSEL/data_batss_msel_vabs/BabyTwins background summary and exclusion 20201102.xlsx')
names(data_exclusion); dim(data_exclusion)

data_background<- readxl::read_xlsx('C:/Users/giobu365/Documents/BT-MSEL/data_batss_msel_vabs/Background_5m.xlsx')

data_pregnancy<- readxl::read_xlsx('C:/Users/giobu365/Documents/BT-MSEL/data_batss_msel_vabs/pregnancy_time_data.xlsx')

### matching participants across files
indx<-match(data_vabs$id,data_demo$`Participant EASE Code`)
data_demo<-data_demo[indx,]

indx_ex<-match(data_vabs$id,data_exclusion$Code)
data_exclusion<-data_exclusion[indx_ex,]

length(which(data_vabs$id!=data_demo$`Participant EASE Code`| data_vabs$id!=data_exclusion$Code | data_demo$`Participant EASE Code`!=data_exclusion$Code))

# add dataset background info

backmatch<-match(data_vabs$id,data_background$Kod)
data_background<-data_background[backmatch,]

length(which(data_vabs$id!=data_background$Kod))

# add data on pregnancy term
indx_preg<-match(data_vabs$id,data_pregnancy$ID)

data_pregnancy<-data_pregnancy[indx_preg,]

### exclusion based on general criteria

excl_indx<-which(data_exclusion$`Exclusion version A`==1)

# clean data files from excluded participants
data_demo_clean<-data_demo[-excl_indx,]
data_vabs_clean<-data_vabs[-excl_indx,]
data_background_clean<-data_background[-excl_indx,]
data_pregnancy_clean<-data_pregnancy[-excl_indx,]

# check pregnancy
preg_day<-data_pregnancy_clean$Gdays
preg_day[which(is.na(preg_day))]<-0

pregnancy_term<-data_pregnancy_clean$Gweeks*7+preg_day

# composite score Communication & Socialization
data_vabs_clean$SC<-data_vabs_clean$Soc_Su+data_vabs_clean$Kom_Su
### dataset included in our analyses

data<-data.frame(cbind(data_vabs_clean$id,data_vabs_clean$age,data_vabs_clean$Fys_Su,data_vabs_clean$SC,data_demo_clean$Gender,pregnancy_term,data_demo_clean$`Bio Mum Age`,data_demo_clean$`Bio Dad Age`,data_demo_clean$`A. Highest level of education`,data_demo_clean$`B. Highest level of education`,data_background_clean$`F16. Ungefär hur hög är familjens* gemensamma månadsinkomst innan skatt (lön och andra ersättningar/bidrag)?`,data_demo_clean$TWAB,data_demo_clean$`Twin pair no.`))
names(data)<-c('id','age','motor','social_comm','sex','term_age','Mum_age','Dad_age','A_edu','B_edu','family_income','TWAB','Twinpair')

# parental education level to MAX within family
data$A_edu_num<-as.numeric(as.factor(data$A_edu))
data$A_edu_num[which(data$A_edu_num==4)]<-6
data$A_edu_num[which(data$A_edu_num==5)]<-4
data$A_edu_num[which(data$A_edu_num==6)]<-5

data$B_edu_num<-as.numeric(as.factor(data$B_edu))
data$B_edu_num[which(data$B_edu_num==4)]<-6
data$B_edu_num[which(data$B_edu_num==5)]<-4
data$B_edu_num[which(data$B_edu_num==6)]<-5

#data$edu_max<-pmax(data$A_edu_num,data$B_edu_num)
data$edu_mean<-rowMeans(cbind(data$A_edu_num,data$B_edu_num),na.rm = T)


# parental age to mean across mum and dad (NA.RM=TRUE for 3 missing Dad age)
data$Mum_age<-as.numeric(data$Mum_age)
data$Dad_age<-as.numeric(data$Dad_age)
data$parental_age<-rowMeans(cbind(data$Mum_age,data$Dad_age),na.rm = T)

## add zygosity
data_exclusion_clean<-data_exclusion[-excl_indx,]
data$zygosity<-data_exclusion_clean$Zygosity

names(data);dim(data)

##### basic data tidying #####

data<-data[,-c(7:10,14,15)] # rm parent A/B info

table(data$zygosity)

# check how many pairs are incomplete by making a variable that counts 
#   the frequency of their pair ID:

data <- merge(data,data.frame(table(Twinpair=data$Twinpair)),by='Twinpair')
table(data$Freq) # check how many pair IDs appear only once 

# missing data

miss_motor<-length(which(is.na(data$motor)))
miss_soccomm<-length(which(is.na(data$social_comm)))
miss_income<-length(which(is.na(data$family_income)))

# carry out the exclusion:

indx<-which(data$Freq==1) # length 4

data<-rbind(data[1:indx[1],],data[indx[1]:indx[2],],data[indx[2]:indx[3],],data[indx[3]:indx[4],],data[indx[4]:length(data$id),])
data <- merge(data,data.frame(table(id=data$id)),by='id')
naindx<-matrix(which(data$Freq.y==2),nrow=2)[1,]

data$motor[naindx]<-NA
data$social_comm[naindx]<-NA

indx1<-which(data$TWAB[naindx]=="1")
indx2<-which(data$TWAB[naindx]=="2")
data$TWAB[naindx[indx1]]<-2
data$TWAB[naindx[indx2]]<-1

dim(data); names(data)

# remove the frequency variable:

data <- data[-c(13,14)]

##### transform data for twin analysis #####

# binary sex: 0=Females; 1=Males.
data$sex[which(data$sex=='Male')]<-1
data$sex[which(data$sex=='Female')]<-0
data$sex<-as.numeric(data$sex)

# ordinal discrete income
data$income<-as.numeric(as.factor(data$family_income))
data$income[which(data$income==2)]<-12
data$income[which(data$income==1)]<-2
data$income[which(data$income==11)]<-NA
data$income[which(data$income==12)]<-11

data<-data[,-c(8)]

# numeric variables
vars <- colnames(data)[c(2:10,12)]
data[vars]<-lapply(data[vars],as.numeric)

### check the distribution of the VABS scales:

vars <- colnames(data)[c(4:5)]
lapply(data[vars],psych::describe)

# reasonable skewness, try without transforming. IF twin assumptions not met, reduce skewness.


# scale variables
vars <- colnames(data)[c(3:5,7,9,10,12)]
data[vars]<-lapply(data[vars],scale)


### control for the effects of covariates: sex, age, parental age, income, education.

library(drgee)

summary(gee( formula = motor~age,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = social_comm~age,data = data, clusterid = "Twinpair", cond = F))
# significant age effect on both

summary(gee( formula = motor~sex,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = social_comm~sex,data = data, clusterid = "Twinpair", cond = F))
# no significant effect of sex

summary(gee( formula = motor~parental_age,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = social_comm~parental_age,data = data, clusterid = "Twinpair", cond = F))
# no significant effect of parental age

summary(gee( formula = motor~term_age,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = social_comm~term_age,data = data, clusterid = "Twinpair", cond = F))
# significant effect of term age on both

summary(gee( formula = motor~edu_mean,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = social_comm~edu_mean,data = data, clusterid = "Twinpair", cond = F))
# no significant effect of education

summary(gee( formula = motor~income,data = data, clusterid = "Twinpair", cond = F))
summary(gee( formula = social_comm~income,data = data, clusterid = "Twinpair", cond = F))
# significant effect of income on motor scale, trend level on social communication


# fit regression models that adjust for significant covariates and then save the residuals:

data$resMotor <- resid(lm(motor~age+term_age+income,data=data,na.action=na.exclude))
data$resSocComm <- resid(lm(social_comm~age+term_age+income,data=data,na.action=na.exclude))

################# data imputation ##############################################################################
library(mice)

data_original<-data

data<-data[-c(13:17)]

stand <- colnames(data)[c(3,4,5,7,9,10,12)]
data[stand]<-lapply(data[stand],as.numeric)


tempData <- mice(data,m=100,maxit=100,meth='cart',seed=2022)

# pooling 
complete_data<-complete(tempData)
data<-complete_data

################################################################################################################

# check that the standardization worked:

stand <- colnames(data)[c(13:14)]
lapply(data[stand],psych::describe)


############################ Randomized twin order ################################################

### give each person a random number (only use two possible numbers, and never 0 and 1):

set.seed(2022)
npairs <- nrow(data)/2
rand <- c(rep(4,npairs),rep(5,npairs))
data$rand <- sample(rand) # randomly choose 1 number for each person

### divide the twins into two subsets: twin 1 and twin 2

twin1 <- subset(data,TWAB==1)
twin2 <- subset(data,TWAB==2)

### now we need to flip twin 2's random number so that it is the opposite of twin 1:

# basically we are going to add twin 1's random number to the twin 2 dataframe and 
#  then recode it so it's the opposite of what it currently is. We have to start by
#   making a dataframe containing just the pair ID number and twin 1's random number:

randVar <- data.frame(cbind(twin1$Twinpair,twin1$rand))
colnames(randVar) <- c('Twinpair','rand')

# now take away twin 2's random number and replace it with the twin 1 random number:

twin2 <- twin2[-c(15)]
twin2 <- merge(twin2,randVar,by='Twinpair')

# now we check that the frequency of each number is the same in both files (i.e.
#   twin 1 and twin 2 should both have the same random number):

table(twin1$rand);table(twin2$rand)

# for twin 1, convert the random numbers to 0s and 1s:

twin1$rand[twin1$rand==4] <- 0
twin1$rand[twin1$rand==5] <- 1

# and then recode twin 2, but in the opposite direction to twin 1:

twin2$rand[twin2$rand==4] <- 1
twin2$rand[twin2$rand==5] <- 0

# check that it worked:

table(twin1$rand);table(twin2$rand) # twin 1 should have as many 0s as twin 2 has 1s

# now join twin 1 and twin 2 back together:

data <- data.frame(rbind(twin1,twin2))
names(data);dim(data)

# sort the data by pairnr and check the first few rows:

data <- data[order(data[,1]),]
head(data)

# split the data by random number:

twin1 <- subset(data,rand==1)
twin2 <- subset(data,rand==0)

# check that the dimensions are the same:

dim(twin1);dim(twin2)

# add '1' to the end of the variables for twin 1 and '2' for twin 2:

twin1labs <- paste(colnames(twin1),1,sep='')
twin2labs <- paste(colnames(twin2),2,sep='')

names(twin1) <- twin1labs
names(twin2) <- twin2labs

##########################################################################################################################

# combine data so that 1 pair each line
dataD <- data.frame(cbind(twin1,twin2))

# remove unused variables:

dataD <- dataD[-c(1,3:10,12,15:27,30)]

# relabel a few variables:

names(dataD)[names(dataD)=='Twinpair1'] <- 'Twinpair'
names(dataD)[names(dataD)=='zygosity1'] <- 'zygosity'

#### done. Save two versions of the data: 

write.csv(data,'phenotypic_resid_data_VABS_NAincome.csv')
write.csv(dataD,'twin_data_VABS_NAincome.csv')
