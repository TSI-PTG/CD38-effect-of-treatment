# Make individual scores trajectories for "non-responders"  ======

m<- data_K1208 |> filter(STUDY_EVALUATION_ID == 7 | STUDY_EVALUATION_ID == 13)

# Plot ABMR Scores ==== 
# Preparing a delta variable 
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(ABMRpm_diff = c(0, diff(ABMRpm))) %>%
  mutate(ABMRpm_diff = if_else(Group == "Index", lead(ABMRpm_diff), ABMRpm_diff)) %>%
  ungroup()


# Shift ABMRpm_diff values to align with the line they should color
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(LineColor = lead(ABMRpm_diff)) %>%
  ungroup()

# Plot
p1<- 
  ggplot(m, aes(x = Group, y = ABMRpm, group = STUDY_EVALUATION_ID, col = LineColor)) +
  geom_line(alpha = 0.5) + 
  geom_point(aes(color = ABMRpm_diff)) +
  theme_bw() +
  labs(x = "", y = bquote("ABMR classifier" ~ ABMR[prob])) +
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15, unit = "pt"), 
        legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_color_gradient2(
    low = "#00ff00bc",  
    mid = "grey60",  
    high = "red",  
    breaks = c(min(m$LineColor, na.rm = TRUE), max(m$LineColor, na.rm = TRUE)))


# Plot ggt0 score  =======
# Preparing a delta variable 
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(ggt0_diff = c(0, diff(ggt0))) %>%
  mutate(ggt0_diff = if_else(Group == "Index", lead(ggt0_diff), ggt0_diff)) %>%
  ungroup()


# Shift ggt0_diff values to align with the line they should color
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(LineColor = lead(ggt0_diff)) %>%
  ungroup()

# Plot
p2<- ggplot(m, aes(x = Group, y = ggt0, group = STUDY_EVALUATION_ID, col = LineColor)) +
  geom_line(alpha = 0.5) + 
  geom_point(aes(color = ggt0_diff)) +
  theme_bw() +
  labs(x = "", y=bquote("Glomerulitis classifier"~ g>0[prob])) +
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15, unit = "pt"), 
        legend.position = "none",          
        panel.grid.major = element_blank(),          
        panel.grid.minor = element_blank()) +
  scale_color_gradient2(
    low = "#00ff00bc",  
    mid = "grey60",  
    high = "red",  
    breaks = c(min(m$LineColor, na.rm = TRUE), max(m$LineColor, na.rm = TRUE)))


# Plot ptcgt0 score  =======
# Preparing a delta variable 
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(ptcgt0_diff = c(0, diff(ptcgt0))) %>%
  mutate(ptcgt0_diff = if_else(Group == "Index", lead(ptcgt0_diff), ptcgt0_diff)) %>%
  ungroup()


# Shift ptcgt0_diff values to align with the line they should color
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(LineColor = lead(ptcgt0_diff)) %>%
  ungroup()

# Plot
p3<- ggplot(m, aes(x = Group, y = ptcgt0, group = STUDY_EVALUATION_ID, col = LineColor)) +
  geom_line(alpha = 0.5) + 
  geom_point(aes(color = ptcgt0_diff)) +
  theme_bw() +
  labs(x = "", y = bquote("Peritubular capillaritis classifier" ~ ptc>0[prob])) +
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15, unit = "pt"), 
        legend.position = "none",          
        panel.grid.major = element_blank(),          
        panel.grid.minor = element_blank()) +
  scale_color_gradient2(
    low = "#00ff00bc",  
    mid = "grey60",  
    high = "red",  
    breaks = c(min(m$LineColor, na.rm = TRUE), max(m$LineColor, na.rm = TRUE)))


# Plot NKB score  =======
# Preparing a delta variable 
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(NKB_diff = c(0, diff(NKB))) %>%
  mutate(NKB_diff = if_else(Group == "Index", lead(NKB_diff), NKB_diff)) %>%
  ungroup()


# Shift NKB_diff values to align with the line they should color
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(LineColor = lead(NKB_diff)) %>%
  ungroup()

# Plot
p4<- ggplot(m, aes(x = Group, y = NKB, group = STUDY_EVALUATION_ID, col = LineColor)) +
  geom_line(alpha = 0.5) + 
  geom_point(aes(color = NKB_diff)) +
  theme_bw() +
  labs(x = "", y = bquote("NK cell burden NKB")) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15, unit = "pt"), 
        legend.position = "none",          
        panel.grid.major = element_blank(),          
        panel.grid.minor = element_blank()) +
  scale_color_gradient2(
    low = "#00ff00bc",  
    mid = "grey60",  
    high = "red",  
    breaks = c(min(m$LineColor, na.rm = TRUE), max(m$LineColor, na.rm = TRUE)))


# Plot DSAST score  =======
# Preparing a delta variable 
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(DSAST_diff = c(0, diff(DSAST))) %>%
  mutate(DSAST_diff = if_else(Group == "Index", lead(DSAST_diff), DSAST_diff)) %>%
  ungroup()


# Shift DSAST_diff values to align with the line they should color
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(LineColor = lead(DSAST_diff)) %>%
  ungroup()

# Plot
p5<- ggplot(m, aes(x = Group, y = DSAST, group = STUDY_EVALUATION_ID, col = LineColor)) +
  geom_line(alpha = 0.5) + 
  geom_point(aes(color = DSAST_diff)) +
  theme_bw() +
  labs(x = "", y = bquote("DSA-selective transcripts (DSAST)")) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15, unit = "pt"), 
        legend.position = "none",          
        panel.grid.major = element_blank(),          
        panel.grid.minor = element_blank()) +
  scale_color_gradient2(
    low = "#00ff00bc",  
    mid = "grey60",  
    high = "red",  
    breaks = c(min(m$LineColor, na.rm = TRUE), max(m$LineColor, na.rm = TRUE)))


# Plot TCMRt score  =======
# Preparing a delta variable 
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(TCMRt_diff = c(0, diff(TCMRt))) %>%
  mutate(TCMRt_diff = if_else(Group == "Index", lead(TCMRt_diff), TCMRt_diff)) %>%
  ungroup()


# Shift TCMRt_diff values to align with the line they should color
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(LineColor = lead(TCMRt_diff)) %>%
  ungroup()

# Plot
p6<- ggplot(m, aes(x = Group, y = TCMRt, group = STUDY_EVALUATION_ID, col = LineColor)) +
  geom_line(alpha = 0.5) + 
  geom_point(aes(color = TCMRt_diff)) +
  theme_bw() +
  labs(x = "", y = bquote("TCMR classifier" ~ TCMR[prob])) +
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15, unit = "pt"), 
        legend.position = "none",          
        panel.grid.major = element_blank(),          
        panel.grid.minor = element_blank()) +
  scale_color_gradient2(
    low = "#00ff00bc",  
    mid = "grey60",  
    high = "red",  
    breaks = c(min(m$LineColor, na.rm = TRUE), max(m$LineColor, na.rm = TRUE)))

# Plot tgt1 score  =======
# Preparing a delta variable 
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(tgt1_diff = c(0, diff(tgt1))) %>%
  mutate(tgt1_diff = if_else(Group == "Index", lead(tgt1_diff), tgt1_diff)) %>%
  ungroup()


# Shift tgt1_diff values to align with the line they should color
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(LineColor = lead(tgt1_diff)) %>%
  ungroup()

# Plot
p7<- ggplot(m, aes(x = Group, y = tgt1, group = STUDY_EVALUATION_ID, col = LineColor)) +
  geom_line(alpha = 0.5) + 
  geom_point(aes(color = tgt1_diff)) +
  theme_bw() +
  labs(x = "", y = bquote("Tubulitis classifier" ~ t>1[prob])) +
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15, unit = "pt"), 
        legend.position = "none",          
        panel.grid.major = element_blank(),          
        panel.grid.minor = element_blank()) +
  scale_color_gradient2(
    low = "#00ff00bc",  
    mid = "grey60",  
    high = "red",  
    breaks = c(min(m$LineColor, na.rm = TRUE), max(m$LineColor, na.rm = TRUE)))

# Plot igt1 score  =======
# Preparing a delta variable 
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(igt1_diff = c(0, diff(igt1))) %>%
  mutate(igt1_diff = if_else(Group == "Index", lead(igt1_diff), igt1_diff)) %>%
  ungroup()


# Shift igt1_diff values to align with the line they should color
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(LineColor = lead(igt1_diff)) %>%
  ungroup()

# Plot
p8<- ggplot(m, aes(x = Group, y = igt1, group = STUDY_EVALUATION_ID, col = LineColor)) +
  geom_line(alpha = 0.5) + 
  geom_point(aes(color = igt1_diff)) +
  theme_bw() +
  labs(x = "", y = bquote("Interstitial infiltrate classifier" ~ i>1[prob])) +
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15, unit = "pt"), 
        legend.position = "none",          
        panel.grid.major = element_blank(),          
        panel.grid.minor = element_blank()) +
  scale_color_gradient2(
    low = "#00ff00bc",  
    mid = "grey60",  
    high = "red",  
    breaks = c(min(m$LineColor, na.rm = TRUE), max(m$LineColor, na.rm = TRUE)))



# Plot TCB score  =======
# Preparing a delta variable 
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(TCB_diff = c(0, diff(TCB))) %>%
  mutate(TCB_diff = if_else(Group == "Index", lead(TCB_diff), TCB_diff)) %>%
  ungroup()


# Shift TCB_diff values to align with the line they should color
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(LineColor = lead(TCB_diff)) %>%
  ungroup()

# Plot
p9<- ggplot(m, aes(x = Group, y = TCB, group = STUDY_EVALUATION_ID, col = LineColor)) +
  geom_line(alpha = 0.5) + 
  geom_point(aes(color = TCB_diff)) +
  theme_bw() +
  labs(x = "", y = bquote("T-cell burden TCB")) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15, unit = "pt"), 
        legend.position = "none",          
        panel.grid.major = element_blank(),          
        panel.grid.minor = element_blank()) +
  scale_color_gradient2(
    low = "#00ff00bc",  
    mid = "grey60",  
    high = "red",  
    breaks = c(min(m$LineColor, na.rm = TRUE), max(m$LineColor, na.rm = TRUE)))


# Plot QCAT score  =======
# Preparing a delta variable 
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(QCAT_diff = c(0, diff(QCAT))) %>%
  mutate(QCAT_diff = if_else(Group == "Index", lead(QCAT_diff), QCAT_diff)) %>%
  ungroup()


# Shift QCAT_diff values to align with the line they should color
m <- m %>%
  mutate(STUDY_EVALUATION_ID = factor(STUDY_EVALUATION_ID)) %>%
  group_by(STUDY_EVALUATION_ID) %>%
  mutate(LineColor = lead(QCAT_diff)) %>%
  ungroup()

# Plot
p10<- ggplot(m, aes(x = Group, y = QCAT, group = STUDY_EVALUATION_ID, col = LineColor)) +
  geom_line(alpha = 0.5) + 
  geom_point(aes(color = QCAT_diff)) +
  theme_bw() +
  labs(x = "", y = bquote("Cytotoxic T-cell-associated QCAT")) +
  #scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15, unit = "pt"), 
        legend.position = "none",          
        panel.grid.major = element_blank(),          
        panel.grid.minor = element_blank()) +
  scale_color_gradient2(
    low = "#00ff00bc",  
    mid = "grey60",  
    high = "red",  
    breaks = c(min(m$LineColor, na.rm = TRUE), max(m$LineColor, na.rm = TRUE)))


# Arrange ====
ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 2, ncol = 5)

# Save ====
ggsave(here::here("graphics", "non_responders_score.tiff"),
       width = 30, height = 20, 
       units = "cm",  dpi = 300, compression = "lzw", type="cairo" )

