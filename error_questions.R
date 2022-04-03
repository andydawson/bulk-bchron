library(dplyr)
library(reshape2)
data = read.csv("chroncontrol_summary_pollen_full.csv")


### Find NA's and future dates for limityounger and age ###

# limityounger
sum(is.na(data$limityounger))

na_young = data[,c("datasetid", "limityounger")]
filter(na_young, is.na(limityounger))

future_young = data[,c("datasetid", "limityounger")]
filter(future_young, limityounger <= -72)

# age
sum(is.na(data$age))

na_age = data[,c("datasetid", "age")]
filter(na_age, is.na(age))

future_age = data[,c("datasetid", "age")]
filter(future_age, age <= -72)


### find dates outside of intcal20 range (95 to 50193) ###
young_range = data[,c("datasetid", "age", "type")]
young_range = filter(young_range, age <95 & type == "Radiocarbon")

old_range = data[,c("datasetid", "age")]
old_range = filter(old_range, age >50193)
