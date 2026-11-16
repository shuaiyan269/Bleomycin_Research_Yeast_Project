#other

library(readxl)

##PMA1_PTK2_mutation_in_clones
#statistical test
data <- read_excel("E:/Bleomycin_Research_Yeast_Project/07_other/input/PMA1_PTK2_mutation_in_clones.xlsx", sheet = 1)
str(data)
tab <- table(data$PMA1_mutation, data$PTK2_mutation)
fisher.test(tab)
