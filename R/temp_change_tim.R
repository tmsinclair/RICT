list.of.packages <- c("devtools",
                      "tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(devtools)
library(tidyverse)
# install experimental rict package
list.of.packages <- c("rict")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install_github("aquaMetrics/rict", dependencies = T)

install("U:/R/R-3.4.4/library/rict-master")

library(rict) # experimental
library(rgdal)
temperature_prediction <- rict_predict(demo_observed_values)
data <- cbind(demo_observed_values,
              "MEAN.AIR.TEMP" = temperature_prediction$MEAN.AIR.TEMP,
              "AIR.TEMP.RANGE" = temperature_prediction$AIR.TEMP.RANGE)
temperature_increases <- seq(0, 2, 0.2)
modelled_temperatures <- map_df(temperature_increases, function(increment){
    data$MEAN.AIR.TEMP <- data$MEAN.AIR.TEMP + increment
    data$AIR.TEMP.RANGE <- data$AIR.TEMP.RANGE + increment
    predict_modelled <- rict_predict(data)
    classify_modelled <- rict_classify(predict_modelled, year_type = "single")
    modelled <- cbind(predict_modelled, classify_modelled)
    modelled$TEMP.INCREASE <- increment
    modelled$SITE  <- as.character(modelled$SITE)
    return(modelled)
})


# Filter on one year
modelled_data <- modelled_temperatures %>%
    filter(YEAR == 2016)
ggplot(modelled_data, aes(x = TEMP.INCREASE, y  = ASPT_aver_spr_aut, color = SITE)) +
    geom_line()
select(modelled_data, SITE, TEMP.INCREASE, SuitCode, belongs_to_end_grp) %>%
    arrange(TEMP.INCREASE)
