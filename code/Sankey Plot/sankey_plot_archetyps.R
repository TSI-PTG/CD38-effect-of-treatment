### Sankey Plot of Archetyp Analysis #### 

#Make a sankey using ggsankey 

#devtools::install_github("davidsjoberg/ggsankey")

library(ggsankey)
library(dplyr)
library(ggplot2)

#load(here::here("data_processed", "felz_final.rda"))


#d<- d |> filter(STUDY_EVALUATION_ID != c(15,18))


# Felz Patients

d_sankey<- d  %>%  filter(Felzartamab == 1) %>% 
  select(IndexBx_Max_Archetype, FU1Bx_Max_Archetype, FU2Bx_Max_Archetype) %>% 
  mutate("Index" = case_when(IndexBx_Max_Archetype == "Early_stage_ABMR" ~ "EABMR", 
                             IndexBx_Max_Archetype == "Fully_developed_ABMR" ~ "FABMR", 
                             IndexBx_Max_Archetype == "Late_stage_ABMR"~ "LABMR", 
                             IndexBx_Max_Archetype == "Non-rejecting"~ "NR"), 
         "Week 24" = case_when(FU1Bx_Max_Archetype == "Early_stage_ABMR" ~ "EABMR", 
                               FU1Bx_Max_Archetype == "Fully_developed_ABMR" ~ "FABMR", 
                               FU1Bx_Max_Archetype == "Late_stage_ABMR"~ "LABMR", 
                               FU1Bx_Max_Archetype == "Non-rejecting"~ "NR"),
         "Week 52" = case_when(FU2Bx_Max_Archetype == "Early_stage_ABMR" ~ "EABMR", 
                               FU2Bx_Max_Archetype == "Fully_developed_ABMR" ~ "FABMR", 
                               FU2Bx_Max_Archetype == "Late_stage_ABMR"~ "LABMR", 
                               FU2Bx_Max_Archetype == "Non-rejecting"~ "NR"))


df <- d_sankey %>%
  make_long("Index", "Week 24", "Week 52")



abmr_arche_felz<- 
  ggplot(df, aes(x = x, 
                 next_x = next_x, 
                 node = node, 
                 next_node = next_node,
                 fill = factor(node),
                 label= node)) +
  geom_sankey(flow.alpha = 0.4, flow.fill = "grey")  + 
  geom_sankey_label(size = 2, fill="white") + 
  theme_sankey() +
  labs(x=c("Biopsy")) +
  scale_fill_manual(values= alpha(c("#bc3c29", "#20854e", "#bc3c29", "#20854e", "#bc3c29","#20854e"),0.8)) + 
  theme(legend.position = "none") + 
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(size = 8.3, color = "black", family = "Arial"), 
        axis.title.x = element_text(size = 10.97, color = "black", vjust = -1.5), 
        legend.text=element_text(size=8, family = "Arial"), 
        plot.title = element_text(size=11, face = "bold", family = "Arial", hjust = 0.2, vjust = -3)) + 
  ggtitle("Felzartamab") 




# Placebo

d_sankey<- d  %>%  filter(Felzartamab == 0) %>% 
  select(IndexBx_Max_Archetype, FU1Bx_Max_Archetype, FU2Bx_Max_Archetype) %>% 
  mutate("Index" = case_when(IndexBx_Max_Archetype == "Early_stage_ABMR" ~ "EABMR", 
                             IndexBx_Max_Archetype == "Fully_developed_ABMR" ~ "FABMR", 
                             IndexBx_Max_Archetype == "Late_stage_ABMR"~ "LABMR", 
                             IndexBx_Max_Archetype == "Non-rejecting"~ "NR"), 
         "Week 24" = case_when(FU1Bx_Max_Archetype == "Early_stage_ABMR" ~ "EABMR", 
                               FU1Bx_Max_Archetype == "Fully_developed_ABMR" ~ "FABMR", 
                               FU1Bx_Max_Archetype == "Late_stage_ABMR"~ "LABMR", 
                               FU1Bx_Max_Archetype == "Non-rejecting"~ "NR"),
         "Week 52" = case_when(FU2Bx_Max_Archetype == "Early_stage_ABMR" ~ "EABMR", 
                               FU2Bx_Max_Archetype == "Fully_developed_ABMR" ~ "FABMR", 
                               FU2Bx_Max_Archetype == "Late_stage_ABMR"~ "LABMR", 
                               FU2Bx_Max_Archetype == "Non-rejecting"~ "NR", 
                               FU2Bx_Max_Archetype == "TCMR"~ "TCMR"))


df <- d_sankey %>%
  make_long("Index", "Week 24", "Week 52")

#df<- df |> drop_na(node)


### Plot 

abmr_arche_placebo <-  ggplot(df, aes(x = x, 
                                      next_x = next_x, 
                                      node = node, 
                                      next_node = next_node,
                                      fill = factor(node),
                                      label= node)) +
  geom_sankey(flow.alpha = 0.4, flow.fill= "grey")  + 
  geom_sankey_label(size = 2, fill="white") + 
  theme_sankey() +
  labs(x="Biopsy") +
  scale_fill_manual(values= alpha(c("#bc3c29", "#20854e", "#bc3c29", "#20854e", "#bc3c29","#20854e"),0.8)) + 
  theme(legend.position = "none") + 
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(size = 8.3, color = "black", family = "Arial"), 
        axis.title.x = element_text(size = 10.97, color = "black", vjust = -1.5), 
        legend.text=element_text(size=8, family = "Arial"), 
        plot.title = element_text(size=11, face = "bold", family = "Arial", hjust = 0.2, vjust = -3)) + 
  ggtitle("Placebo")  
#coord_cartesian(ylim = c(-7.5,6.4)) 



### Combine in one graphic

sankey_abmr_arche<- ggpubr::ggarrange(abmr_arche_placebo, abmr_arche_felz,  ncol = 2, nrow = 1)

ggsave(here::here("graphics", "sankey_abmr_activity.tiff"),
       width = 20, height = 10, 
       units = "cm",  dpi = 300, compression = "lzw", type="cairo" )
