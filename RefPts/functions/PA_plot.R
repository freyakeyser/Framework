library(tidyverse)
library(ggthemes)
library(cowplot)

theme_set(theme_tufte(base_size = 20))

d=data.frame(x1=c(0,0.2,0.4), x2=c(0.2,0.4,1), y1=c(0,0,0), y2=c(1,1,1), 
             t=factor(c('Critical Zone','Cautious Zone','Healthy Zone'),
                      levels=c('Critical Zone','Cautious Zone','Healthy Zone')), 
             r=c("firebrick2","goldenrod1",'chartreuse2'))
ln <- data.frame(x=c(0.0,0.2,0.2,0.4,0.4,1),y=c(0.01,0.01,0.01,0.2,0.2,0.2),grp = c(rep("Critical Zone",2),
                                                                      rep("Cautious Zone",2),
                                                                      rep("Healthy Zone",2)))
labs <- data.frame(x=c(0.2,0.4,0.55),y=c(0.5,0.5,0.2),label = c(1,2,"3a"))

rp.plt <- ggplot() + 
                  geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t,color=t)) +
                  geom_text(data=d, aes(x=x1+(x2-x1)/2, y=1.015, label=t),size=9) +
                  geom_vline(xintercept = c(0.2,0.4),size=1.5,color='darkgrey') + 
                  geom_line(data=ln,aes(x=x,y=y,group=grp),size=1.5,linetype  ='dashed') + 
                  geom_label(data = labs,aes(x=x,y=y,label=label),size=10) +
                  theme(axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")),
                        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")),
                        #axis.title.y = element_text(angle = 0),
                        legend.position = 'none') +
                  scale_x_continuous(limits=c(0, 1.05), expand = c(0, 0),breaks = NULL,name="Stock Status") +
                  scale_y_continuous(limits=c(0, 1.05), expand = c(0, 0),breaks=NULL,name="Removal Rate") +
                  scale_fill_manual(values = d$r) + scale_color_manual(values = d$r) 

save_plot("D:/Framework/SFA_25_26_2024/RefPts/Figures/Generic_RP_plot.png",base_height = 10,base_width = 15,plot = rp.plt)




d=data.frame(x1=c(0,0.2,0.4,0.7), x2=c(0.2,0.4,0.7,1), y1=c(0,0,0,0), y2=c(1,1,1,1), 
             t=factor(c('Critical Zone','Cautious Zone','Healthy Zone below TRP','Healthy Zone above TRP'),
                      levels=c('Critical Zone','Cautious Zone','Healthy Zone below TRP','Healthy Zone above TRP')), 
             r=c("firebrick2","goldenrod1",'chartreuse2','chartreuse2'))
ln <- data.frame(x=c(0.0,0.2,0.2,0.4,0.4,0.7,0.7,1),y=c(0.01,0.01,0.01,0.2,0.2,0.2,0.3,0.3),grp = c(rep("Critical Zone",2),
                                                                                                    rep("Cautious Zone",2),
                                                                                                    rep("Healthy Zone",2),
                                                                                                    rep("Healthy Zone 2",2)))
labs <- data.frame(x=c(0.2,0.4,0.55,0.7,0.85),y=c(0.5,0.5,0.2,0.5,0.3),label = c(1,2,"3a",4,"3b"))

rp.plt2 <- ggplot() + 
  geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t,color=t)) +
  geom_text(data=d, aes(x=x1+(x2-x1)/2, y=1.015, label=t),size=9) +
  geom_vline(xintercept = c(0.2,0.4,0.7),size=1.5,color='darkgrey') + 
  geom_line(data=ln,aes(x=x,y=y,group=grp),size=1.5,linetype  ='dashed') + 
  geom_label(data = labs,aes(x=x,y=y,label=label),size=10) +
  theme(axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")),
        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")),
        #axis.title.y = element_text(angle = 0),
        legend.position = 'none') +
  scale_x_continuous(limits=c(0, 1.05), expand = c(0, 0),breaks = NULL,name="Stock Status") +
  scale_y_continuous(limits=c(0, 1.05), expand = c(0, 0),breaks=NULL,name="Removal Rate") +
  scale_fill_manual(values = d$r) + scale_color_manual(values = d$r) 

save_plot("D:/Framework/SFA_25_26_2024/RefPts/Figures/TRP_RP_plot.png",base_height = 10,base_width = 15,plot = rp.plt2)
