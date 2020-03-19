#RICT model processing results from model

#load necessary packages
library(dplyr)
library(ggplot2)

#read in the data
results <- read.csv("./input/rict_compare_results.csv")

#select just probability comparisons for the difference in water quality
comp_probs <- results %>%
    select(starts_with("Probability.Result.A.in."))

#The water quality categories
categories <-  c("High", "Good", "Moderate", "Poor", "Bad")

#Split each quality and site combination into a list
comp_probs_list <- split(comp_probs, seq(nrow(comp_probs)))

#Give them sensible names
names(comp_probs_list) <- paste0(results$EQR.metric.compared,"_Result.A.",results$Result.A,"_Result.B.",results$Result.B)

#####
#create a function to turn the results into matrix format (although data frame in R)
quality_matrix <- function(x,y){
 z <- NULL

   for(i in y){
       w <- x %>%
           select(starts_with(paste0("Probability.Result.A.in.", i)))

        colnames(w) <- paste0("Probability.Result.A.in.",y)

       z <- rbind(z,w)

   }
   rownames(z) <- paste0("Probability.Result.B.in.",y)
    print(z)
}

#repeat the function for all the EQR and site combinations
quality_matrices <- lapply(comp_probs_list, quality_matrix, y = categories)

#####
#Visualising the matrices

for(j in 1:length(comp_probs_list)){
    comp_probs_list[j]
}

plot_quality_comparison <- function(probability, categories, names){
    #Expand the combinations of quality categorisations into rows
    expanded_categories <- expand.grid(Result_A = categories, Result_B = categories)

    expanded_categories$Probability <- as.numeric(probability)

    probability_plot <- ggplot(expanded_categories, aes(x = Result_A, y = Result_B, fill= Probability)) +
        geom_tile()+
        scale_fill_gradient(limits = c(0, 100))

    ggsave(paste0("./figures/quality_probability_matrix_",names,".jpg"), device = "jpg", plot = probability_plot, width = 10, height = 8)
}

lapply(comp_probs_list, plot_quality_comparison, categories = categories, names = names(comp_probs_list))

