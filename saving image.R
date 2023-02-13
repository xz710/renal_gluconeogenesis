png(file="C:/Users/Xinyi/Desktop/omics/human kidney gluconeogenesis/gluconeogenesis expression plot.png",
    width=300, height=600)
ggplot(dt, (aes(x=Group, y=Gene,
                color = Percentage_of_PCT_cell, 
                size=Expression_Level))) + 
  geom_point() +
  theme(axis.text = element_text(face="bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

png(file="C:/Users/Xinyi/Desktop/omics/human kidney gluconeogenesis/gluconeogenesis expression plot2.png",
    width=300, height=500)
ggplot(dt, (aes(x=Group, y=Gene,
                color = Expression_Level, 
                size=Percentage_of_PCT_cell))) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()