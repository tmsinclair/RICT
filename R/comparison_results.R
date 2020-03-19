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
    return(z)
}


quality_matrix(comp_probs_list$`1`, y = categories)
#repeat the function for all the EQR and site combinations
quality_matrices <- lapply(comp_probs_list, quality_matrix, y = categories)

#####
#Visualising the matrices

#Expand the combinations of quality categorisations into rows
expanded_categories <- expand.grid(Result_A = categories, Result_B = categories)
j <- 1
for(j in 1:length(comp_probs_list)){
    #attach the probabilities into the site quality combinations
    expanded_categories$Probability <- as.numeric(unlist(comp_probs_list[j]))

    #visualise the data in a matrix
    probability_plot <- ggplot(expanded_categories, aes(x = Result_A, y = Result_B, fill= Probability)) +
        geom_tile()+
        scale_fill_gradientn(colours = c("#340044","#2C4279","#1E8644","#64CA44","#FFE51E"), limits = c(0, 100))

    #save it as a jpg file
    ggsave(paste0("./figures/quality_probability_matrix_",as.character(names(comp_probs_list)[j]),".jpg"), device = "jpg", plot = probability_plot, width = 10, height = 8)
}

