library(tidyverse)
library(rhdf5)
library(edgeR)

plot_data <- read_tsv("C:/Users/Xinyi/Desktop/omics/human kidney gluconeogenesis/dotplot expression.txt")


diabetic_dt <- plot_data[,1:3]
colnames(diabetic_dt) <- c("Gene", "Expression Level", "Percentage of PCT cell")
diabetic_dt <- diabetic_dt %>%
  mutate("Group" = "Diabetic") 


control_dt <- cbind(plot_data[1],plot_data[,4:5])
colnames(control_dt) <- c("Gene", "Expression Level", "Percentage of PCT cell")
control_dt <- control_dt %>%
  mutate("Group" = "non-Diabetic") 

dt <- rbind(control_dt,diabetic_dt)
Expression_Level <- as.numeric(dt$'Expression Level')
Percentage_of_PCT_cell <- as.numeric(dt$'Percentage of PCT cell')

library(RColorBrewer)
dt$Group <- factor(dt$Group,levels = c("non-Diabetic", "Diabetic"))

(dotplot <- ggplot(dt, (aes(x=Group, y=Gene,
                           color = Expression_Level, 
                           size=Percentage_of_PCT_cell))) + 
  geom_point() +
    theme(axis.text = element_text(face="bold"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
