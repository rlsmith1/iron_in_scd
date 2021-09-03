

# Our main question is whether diet had an effect on each endpoint; a secondary question is whether there were sex differences.




# libraries ---------------------------------------------------------------

    library(tidyverse)
    library(tidymodels)
    library(ggforce)
    library(ggrepel)
    library(MASS)
    library(FactoMineR)
    library(factoextra)



# data --------------------------------------------------------------------


    # import
    df_mice <- read.csv("data/Majed's mice.csv") %>% 
      as_tibble() %>% 
      dplyr::select(-c(X, X.1, X.2)) %>% 
      filter(!is.na(Mouse.ID))
    
    # format
    df_mice <- df_mice %>% 
      mutate_if(is.character, factor) %>% 
      mutate_if(is.integer, factor)


    
#### Models   

    
# overall plot 
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      dplyr::select(-c(hepatic.necrosis..yes.no., iron.deposition.score..0.4...spleen., iron.deposition.score..0.4...liver., iron.deposition.score..0.4..kidney.)) %>% 
      pivot_longer(cols = 5:ncol(.), names_to = "variable", values_to = "value") %>% 
      
      ggplot(aes(x = treatment, y = value)) +
      geom_violin(aes(color = treatment)) +
      facet_grid(variable ~ timepoint..month., scales = "free_y") +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    

# ferritin ----------------------------------------------------------------


    # identify outliers?
    df_mice %>% 
      group_by(treatment) %>% 
      mutate(Q1 = quantile(Ferritin..ng.ml., 0.25, na.rm = TRUE), 
             Q3 = quantile(Ferritin..ng.ml., 0.75, na.rm = TRUE),
             IQR = Q3 - Q1,
             high = Q3 + 1.5*IQR,
             low = Q1 - 1.5*IQR,
             outlier = ifelse(Ferritin..ng.ml. > high | Ferritin..ng.ml. < low, 1, 0)) %>% 
      ungroup() %>% 
      dplyr::select(c(Mouse.ID, sex, treatment, timepoint..month., Ferritin..ng.ml., outlier))
    
    # model
    aov_ferritin <- aov(Ferritin..ng.ml. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_ferritin %>% tidy()
    aov_ferritin %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05) %>% print(n = nrow(.))
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = Ferritin..ng.ml.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   treatment effect
    #   female sufficient different than deficient, esp after 6 months

    
# HAMP1 -------------------------------------------------------------------

    # model
    aov_hamp1 <- aov(HAMP1..cDNA.ng. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_hamp1 %>% tidy()
    aov_hamp1 %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05) %>% print(n = nrow(.))
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = HAMP1..cDNA.ng.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   treatment effect (across all three timepoints, sufficient different than deficient)
    #   timepoint-sex interaction (M different from F at 6 mos)
    #   timepoint-sex-treatment interaction
    #       2: F sufficient different than F deficient
    #       4: M sufficient different than M deficient
    #       6: F sufficient different than F deficient
    

# Hepcidin ----------------------------------------------------------------

    
    aov_hepcidin <- aov(Hepcidin..ng.ml. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_hepcidin %>% tidy()
    aov_hepcidin %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05) %>% print(n = nrow(.))
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = Hepcidin..ng.ml.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   treatment effect
    #   timepoint-treatment interaction (timepoints 4 & 6, sufficient different from deficient)
    #   timepoint-sex interaction
    #       F: 4 different from 6
    #       M: 2 different from 4
    #   timepoint-sex-treatment interaction
    #       6: F sufficient different than F deficient
    #       F sufficient: 4 different from 6
    #       4 sufficient: M different from F
    #       M: 4 sufficient different than 4 deficient
    #       M sufficient: 2 different than 4; 4 different than 6
    
    

# Plasma heme -------------------------------------------------------------


    aov_plasma_heme <- aov(Plasma.Heme..µM. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_plasma_heme %>% tidy()
    aov_plasma_heme %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05) %>% print(n = nrow(.))
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = Plasma.Heme..µM.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   timepoint effect (4 diff from 2, 6 diff from 2; esp in deficient)
    #   treatment effect
    #   sex effect
    #   timepoint-treatment interaction (timepoints 4, sufficient different from deficient)
    #   timepoint-sex interaction
    #       F: 4 & 6 different from 2
    #       2: M different from F


# Haptoglobin ------------------------------------------------------------


    aov_haptoglobin <- aov(Haptoglobin..ng.ml. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_haptoglobin %>% tidy()
    aov_haptoglobin %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05)
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = Haptoglobin..ng.ml.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   sex effect
    #   timepoint-sex interaction
    #       M: 4 & 6 different from 2
    #       2: M different from F
    
    

# Hemopexin ---------------------------------------------------------------


    aov_hemopexin <- aov(Hemopexin..ng.ml. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_hemopexin %>% tidy()
    aov_hemopexin %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05) %>% print(n = nrow(.))
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = Hemopexin..ng.ml.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   timepoint effect (4 different from 2)
    #   sex effect
    #   timepoint-sex interaction
    #       M: 4 different from 2
    #       2: M different from F
    

    
# P50 ---------------------------------------------------------------------

# ***not enough data for timepoint comparison (only timepoint 6)
    aov_p50 <- aov(P50..mm.Hg. ~ treatment*sex, data = df_mice)
    aov_p50 %>% tidy()
    aov_p50 %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05)
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = P50..mm.Hg.)) +
      geom_violin(aes(color = treatment)) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   treatment effect (4 different from 2)
    #   sex effect


# Round cell half-life ----------------------------------------------------

  
    aov_round_cell <- aov(Round.cell.half.life..minutes. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_round_cell %>% tidy()
    aov_round_cell %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05) %>% print(n = nrow(.))
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = Round.cell.half.life..minutes.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    
    # conclusions:
    #   timepoint effect (2 & 6 different from 4)
    #   treatment effect (4 different from 2)
    #   treatment-sex interaction
    #       F & M: 2 & 6 different from 4
    

    
# RBC count ---------------------------------------------------------------

    
    aov_rbc <- aov(RBC..10.6.uL. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_rbc %>% tidy()
    aov_rbc %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05)
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = RBC..10.6.uL.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   timepoint effect (6 different from 4)
    #   treatment effect

    

# Hemoglobin --------------------------------------------------------------

    
    aov_hb <- aov(Hb..g.dL. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_hb %>% tidy() %>% filter(p.value < 0.05)
    aov_hb %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05)
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = Hb..g.dL.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   timepoint effect (6 different from 4)
    #   treatment effect
    #   sex effect
    #   timepoint-treatment interaction
    #       deficient: 2 different from 4
    #       4: suff different from def
    #   timepoint-sex interaction
    #       deficient: M different from F
    #       M: sufficient different from deficient
    

    
# HCT ---------------------------------------------------------------------


    aov_hct <- aov(HCT.... ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_hct %>% tidy() %>% filter(p.value < 0.05)
    aov_hct %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05)
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = HCT....)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   timepoint effect (2 different from 4)
    #   treatment effect
    #   sex effect

    

# MCV ---------------------------------------------------------------------


    aov_mcv <- aov(MCV..fL. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_mcv %>% tidy() %>% filter(p.value < 0.05)
    aov_mcv %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05)
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = MCV..fL.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   timepoint effect (2 different from 4, 2 different from 6, 4 different from 6)

    

# MCH ---------------------------------------------------------------------

    
    aov_mch <- aov(MCH..pg. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_mch %>% tidy() %>% filter(p.value < 0.05)
    aov_mch %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05)
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = MCH..pg.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   timepoint effect (2 different from 4, 2 different from 6)
    #   sex effect
    #   treatment-sex interaction
    #       deficient: M different from F
    #       M: sufficient different from deficient
    
    

    
# MCHC --------------------------------------------------------------------


    aov_mchc <- aov(MCHC..g.dL. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_mchc %>% tidy() %>% filter(p.value < 0.05)
    aov_mchc %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05)
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = MCHC..g.dL.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   timepoint effect (6 different from 2, 6 different from 4)
    #   timepoint-sex interaction effect
    #       2: M different from F
    #       M: 6 different from 2 & 4
    #   treatment-sex interaction
    #       deficient: M different from F
    #       M: sufficient different from deficient
    

    
# MCHC-O ------------------------------------------------------------------


    aov_mchc_o <- aov(MCHC.O..g.dL. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_mchc_o %>% tidy() %>% filter(p.value < 0.05)
    aov_mchc_o %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05) %>% print(n = nrow(.))
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = MCHC.O..g.dL.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   timepoint effect (6 different from 2, 6 different from 4)
    #   treatment effect
    #   timepoint-sex interaction 
    #       F: 6 different from 4
    #       M: 6 different from 2 & 4
    #   treatment-sex interaction
    #       deficient: M different from F
    #       M: sufficient different from deficient
    
    

# RDW ---------------------------------------------------------------------


    aov_rdw <- aov(RDW.... ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_rdw %>% tidy() %>% filter(p.value < 0.05)
    aov_rdw %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05)
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = RDW....)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   treatment-sex interaction
    #       deficient: M different from F
    #       M: sufficient different from deficient
    
    

# RET ---------------------------------------------------------------------


    aov_ret <- aov(Ret.He..pg. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_ret %>% tidy() %>% filter(p.value < 0.05)
    aov_ret %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05)
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = Ret.He..pg.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   sex effect
    #   timepoint-sex interaction (4: M different from F)
    #   treatment-sex interaction
    #       deficient: M different from F
    #       M: sufficient different from deficient
    
    

# ARC ---------------------------------------------------------------------


    aov_arc <- aov(ARC..10.6.uL. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_arc %>% tidy() %>% filter(p.value < 0.05)
    aov_arc %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05)
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = ARC..10.6.uL.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   timepoint effect (6 different from 2)
    #   treatment effect 
    
    

# Ret-He ------------------------------------------------------------------


    aov_ret_he <- aov(Ret.He..pg. ~ timepoint..month.*treatment*sex, data = df_mice)
    aov_ret_he %>% tidy() %>% filter(p.value < 0.05)
    aov_ret_he %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05)
    
    # plot
    cols <- c("2" = "#828282", "4" = "#828282", "6" = "#828282", "M" = "#0072B2", "F" = "#CC79A7")
    
    df_mice %>% 
      ggplot(aes(x = treatment, y = Ret.He..pg.)) +
      geom_violin(aes(color = treatment)) +
      facet_wrap(~timepoint..month.) +
      # geom_boxplot(aes(color = treatment), width = .3, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.1) +
      stat_summary(fun = "median", geom = "crossbar", aes(color = timepoint..month.), position = position_dodge(), size = .3, lty = 2) +
      
      geom_sina(aes(fill = timepoint..month., color = sex), size = 2.5, shape = 1) +
      stat_summary(fun = "median", geom = "crossbar", aes(fill = timepoint..month., color = sex), position = position_dodge(), size = .2) +
      
      scale_color_manual(values = cols) +
      theme_bw() +
      theme(legend.position = "none")
    
    # conclusions:
    #   sex effect
    #   timepoint-sex interaction effect
    #       4: M different from F
    #   treatment-sex interaction effect 
    #       deficient: M different from F
    #       M: sufficient different from deficient
    
    

# hepatic necrosis --------------------------------------------------------


    # select variables of interest
    df_hep_nec <- df_mice %>% 
      dplyr::select(c(timepoint..month., treatment, sex, hepatic.necrosis..yes.no.)) %>% 
      dplyr::filter(!is.na(hepatic.necrosis..yes.no.))
    
    # regression
    glm(hepatic.necrosis..yes.no. ~ timepoint..month. + treatment + sex, data = df_hep_nec, family = "binomial") %>% 
      summary() %>% 
      .$coefficients %>% 
      as_tibble(rownames = "term") %>% 
      mutate(odds_ratio = exp(Estimate))
    
    # plot
    df_hep_nec %>% 
      # count(treatment, timepoint..month., sex, hepatic.necrosis..yes.no.) %>% 
      filter(!is.na(hepatic.necrosis..yes.no.)) %>% 
      
      ggplot(aes(x = treatment, y = hepatic.necrosis..yes.no., color = hepatic.necrosis..yes.no.)) +
      geom_point(position = position_jitter(0.2)) +
      facet_grid(sex ~ timepoint..month.) +
      
      theme_bw()
    
    # conclusions:
    # treatment is the only significant variable, with sufficient treatment associated with ~4x increased likelihood of hepatic necrosis
  
      

# iron deposition score ---------------------------------------------------

    
    # spleen
    
        # select variables of interest
        df_spleen <- df_mice %>% 
          dplyr::select(c(timepoint..month., treatment, sex, iron.deposition.score..0.4...spleen.)) %>% 
          dplyr::filter(!is.na(iron.deposition.score..0.4...spleen.))
    
        # regression (select model with lowest AIC)
        model_spleen <- polr(iron.deposition.score..0.4...spleen. ~ timepoint..month.*sex + treatment, data = df_spleen, Hess = TRUE, method = "probit")
        coef_spleen <- model_spleen %>% summary() %>% coef() %>% as_tibble(rownames = "term")
        coef_spleen %>% mutate(p.val = pnorm(abs(coef_spleen$`t value`), lower.tail = FALSE)*2, OR = exp(Value))

        # plot
        df_spleen %>% ggplot(aes(x = treatment, y = iron.deposition.score..0.4...spleen., color = iron.deposition.score..0.4...spleen.)) +
          geom_point(position = position_jitter(0.2)) +
          facet_grid(sex ~ timepoint..month., margins = TRUE) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))          
    
    # liver
        
        # select variables of interest
        df_liver <- df_mice %>% 
          dplyr::select(c(timepoint..month., treatment, sex, iron.deposition.score..0.4...liver.)) %>% 
          dplyr::filter(!is.na(iron.deposition.score..0.4...liver.))
        
        # regression (select model with lowest AIC)
        model_liver <- polr(iron.deposition.score..0.4...liver. ~ timepoint..month.*sex + treatment, data = df_liver, Hess = TRUE, method = "probit")
        coef_liver <- model_liver %>% summary() %>% coef() %>% as_tibble(rownames = "term")
        coef_liver %>% mutate(p.val = pnorm(abs(coef_liver$`t value`), lower.tail = FALSE)*2, OR = exp(Value))
        
        polr(iron.deposition.score..0.4...liver. ~ timepoint..month.*sex*treatment, data = df_liver, Hess = TRUE, method = "probit")
        polr(iron.deposition.score..0.4...liver. ~ timepoint..month. + sex + treatment, data = df_liver, Hess = TRUE, method = "probit")
        polr(iron.deposition.score..0.4...liver. ~ timepoint..month.*sex + treatment, data = df_liver, Hess = TRUE, method = "probit")
        polr(iron.deposition.score..0.4...liver. ~ timepoint..month.*treatment + sex, data = df_liver, Hess = TRUE, method = "probit")
        
        
        # plot
        df_liver %>% ggplot(aes(x = treatment, y = iron.deposition.score..0.4...liver., color = iron.deposition.score..0.4...liver.)) +
          geom_point(position = position_jitter(0.2)) +
          facet_grid(sex ~ timepoint..month., margins = TRUE) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))          
        
    
    # kidney
        
        # select variables of interest
        df_kidney <- df_mice %>% 
          dplyr::select(c(timepoint..month., treatment, sex, iron.deposition.score..0.4..kidney.)) %>% 
          dplyr::filter(!is.na(iron.deposition.score..0.4..kidney.))
        
        # regression (select model with lowest AIC)
        model_kidney <- polr(iron.deposition.score..0.4..kidney. ~ timepoint..month. + treatment + sex, data = df_kidney, Hess = TRUE, method = "probit")
        coef_kidney <- model_kidney %>% summary() %>% coef() %>% as_tibble(rownames = "term")
        coef_kidney %>% mutate(p.val = pnorm(abs(coef_kidney$`t value`), lower.tail = FALSE)*2, OR = exp(Value))
        
        # plot
        df_kidney %>% ggplot(aes(x = treatment, y = iron.deposition.score..0.4..kidney., color = iron.deposition.score..0.4..kidney.)) +
          geom_point(position = position_jitter(0.2)) +
          facet_grid(sex ~ timepoint..month., margins = TRUE) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))          
        
    
        polr(iron.deposition.score..0.4..kidney. ~ timepoint..month.*sex*treatment, data = df_kidney, Hess = TRUE, method = "probit")
        polr(iron.deposition.score..0.4..kidney. ~ timepoint..month.*sex + treatment, data = df_kidney, Hess = TRUE, method = "probit")
        polr(iron.deposition.score..0.4..kidney. ~ timepoint..month. + sex*treatment, data = df_kidney, Hess = TRUE, method = "probit")
        polr(iron.deposition.score..0.4..kidney. ~ timepoint..month.*treatment + sex, data = df_kidney, Hess = TRUE, method = "probit")
        polr(iron.deposition.score..0.4..kidney. ~ timepoint..month. + treatment + sex, data = df_kidney, Hess = TRUE, method = "probit")
        polr(iron.deposition.score..0.4..kidney. ~ timepoint..month. + treatment, data = df_kidney, Hess = TRUE, method = "probit")
        polr(iron.deposition.score..0.4..kidney. ~ timepoint..month. + sex, data = df_kidney, Hess = TRUE, method = "probit")
        polr(iron.deposition.score..0.4..kidney. ~ treatment + sex, data = df_kidney, Hess = TRUE, method = "probit")
        
    
                                                       
    
    
    
    
    
    
    
        
            
    
    
    