

# Our main question is whether diet had an effect on each endpoint; a secondary question is whether there were sex differences.




# libraries ---------------------------------------------------------------

    library(tidyverse)
    library(tidymodels)
    library(corrr)
    library(ggforce)
    library(ggrepel)
    library(ggpubr)
    library(MASS)
    library(FactoMineR)
    library(factoextra)
    library(ggfortify)



# data --------------------------------------------------------------------


    # import
    df_mice <- read.csv("data/Majed's mice.csv") %>% 
      as_tibble() %>% 
      dplyr::select(-c(X, X.1, X.2)) %>% 
      filter(!is.na(Mouse.ID))
    
    # format
    df_mice <- df_mice %>% 
      mutate_if(is.character, factor) %>% 
      mutate_if(is.integer, as.numeric) %>% 
      mutate(timepoint..month. = as.factor(timepoint..month.)) %>% 
      mutate(hepatic.necrosis..yes.no. = factor(ifelse(hepatic.necrosis..yes.no. == "Y", 1, 0)))



# Exploratory analyses ----------------------------------------------------


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
    
    
    # correlations
   
        df_correlate <- df_mice %>% dplyr::select(-c(1:4, 24)) 
        mice_cor <- df_correlate %>% correlate(method = "spearman", diagonal = NA)
        
        # plot
        rplot(mice_cor, colors = c("deepskyblue1", "white", "coral"), print_cor = TRUE) +
          theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
        
        # test for significance
        f_calc_ttest_p_value <- function(vec_a, vec_b){
          
          t.test(vec_a, vec_b)$p.value
          
        }
        
        mice_cor_long <- mice_cor %>% stretch()
        
        colpair_map(df_correlate, f_calc_ttest_p_value) %>% 
          pivot_longer(2:ncol(.), names_to = "y", values_to = "p.val") %>% 
          .$p.val %>% 
          p.adjust(method = "bonferroni") %>% 
          as_tibble() %>% 
          
          bind_cols(mice_cor_long) %>% 
          dplyr::rename("FDR" = "value") %>% 
          dplyr::select(c(x, y, r, FDR))
          
        
    # PCA
        
        # format data & impute
        df_mice_long <- df_mice %>% 
          dplyr::select(-c(P50..mm.Hg., iron.deposition.score..0.4...spleen., iron.deposition.score..0.4...liver., iron.deposition.score..0.4..kidney.)) %>%  # too many missing values        
          pivot_longer(5:ncol(.), values_to = "value", names_to = "variable")
        
        df_pca_long <- df_mice_long %>% 
          group_by(sex, timepoint..month., treatment, variable) %>% 
          summarise(median(value, na.rm = TRUE)) %>% # calculate median of each group
          left_join(df_mice_long) %>% 
          mutate(value = ifelse(is.na(value), `median(value, na.rm = TRUE)`, value)) # impute group median
        
        df_pca <- df_pca_long %>% 
          ungroup() %>% 
          dplyr::select(-5) %>% 
          pivot_wider(id_cols = c(1:3,5), names_from = variable, values_from = value) %>% 
          dplyr::select(-c(1:4)) %>% 
          as.data.frame()
        
        rownames(df_pca) <- df_mice$Mouse.ID
          
        # run PCA
        pca_mice <- prcomp(df_pca, na.rm = TRUE)
        fviz_eig(pca_mice)

        df_pca_individuals <- pca_mice$x %>% 
          as_tibble(rownames = "Mouse.ID") %>% 
          right_join(df_mice_long %>% dplyr::select(-c(variable, value)) %>% mutate(Mouse.ID = as.character(Mouse.ID))) %>% 
          unique() %>% 
          mutate(treatment_timepoint = paste0(treatment, sep = "_", timepoint..month.),
                 sex_timepoint = paste0(sex, sep = "_", timepoint..month.),
                 sex_treatment = paste0(sex, sep = "_", treatment)) 
        
        df_pca_variables <- pca_mice$rotation %>% as_tibble(rownames = "variable")
          
        # plot
        ggplot() +
          geom_point(data = df_pca_individuals, aes(x = PC1, y = PC2, color = sex_treatment)) +
          stat_conf_ellipse(data = df_pca_individuals, aes(x = PC1, y = PC2, color = sex_treatment), alpha = 0.7) +

          geom_point(data = df_pca_variables, aes(x = PC1*2^14*1.5, y = PC2*2^14*1.5)) +
          geom_text_repel(data = df_pca_variables, aes(x = PC1*2^14*1.5, y = PC2*2^14*1.5, label = variable)) +
          geom_segment(data = df_pca_variables %>% filter(variable %in% c("Ferritin..ng.ml.", "Haptoglobin..ng.ml.", "Hemopexin..ng.ml.")), 
                       aes(x = 0, y = 0, xend = PC1*2^14*1.5, yend = PC2*2^14*1.5), arrow = arrow()) +

          theme_bw() 
        

    
#### Models   

    
    
    

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
    
    lm(Hepcidin..ng.ml. ~ timepoint..month.*treatment*sex, data = df_mice) %>% autoplot()
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
    aov_p50 <- aov(P50..mm.Hg. ~ treatment + sex, data = df_mice)
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

  
    aov_round_cell <- aov(Round.cell.half.life..minutes. ~ timepoint..month. + treatment*sex, data = df_mice)
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

    
    aov_rbc <- aov(RBC..10.6.uL. ~ timepoint..month.*treatment, data = df_mice)
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


    aov_hct <- aov(HCT.... ~ timepoint..month.*treatment + sex, data = df_mice)
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


    aov_mcv <- aov(MCV..fL. ~ timepoint..month.*sex, data = df_mice)
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

    
    aov_mch <- aov(MCH..pg. ~ timepoint..month. + treatment*sex, data = df_mice)
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


    aov_ret <- aov(RET. ~ timepoint..month.*treatment*sex, data = df_mice)
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


    aov_arc <- aov(Ret.He..pg. ~ timepoint..month. + treatment*sex, data = df_mice)
    aov_arc %>% tidy() %>% filter(p.value < 0.05)
    aov_arc %>% TukeyHSD() %>% tidy() %>% filter(adj.p.value < 0.05)
    
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
    glm(hepatic.necrosis..yes.no. ~ treatment + sex, data = df_hep_nec, family = "binomial") %>% 
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
    
        # ordinal regression (select model with lowest AIC)
        ord_spleen <- polr(iron.deposition.score..0.4...spleen. ~ timepoint..month.*sex + treatment, data = df_spleen, Hess = TRUE, method = "probit")
        coef_spleen <- ord_spleen %>% summary() %>% coef() %>% as_tibble(rownames = "term")
        coef_spleen %>% mutate(p.val = pnorm(abs(coef_spleen$`t value`), lower.tail = FALSE)*2, OR = exp(Value))
        
        df_spleen %>% count(treatment, iron.deposition.score..0.4...spleen.)
        
        # linear regression
        lm_spleen <- lm(iron.deposition.score..0.4...spleen. ~ timepoint..month. + treatment + sex, data = df_spleen)
        lm_spleen %>% summary()
        
        aov_spleen <- aov(iron.deposition.score..0.4...spleen. ~ timepoint..month. + treatment + sex, data = df_spleen)
        aov_spleen %>% tidy()
        
          # plot residuals
          augment(lm_spleen) %>% 
            ggplot(aes(x = .fitted, y = .resid)) +
            geom_point()
          
          autoplot(lm_spleen)
          
          

        # plot
        df_spleen %>% mutate(timepoint..month. = as.numeric(as.character(timepoint..month.)),
                             iron.deposition.score..0.4...spleen. = as.numeric(as.character(iron.deposition.score..0.4...spleen.)),
                             sex_treatment = factor(paste0(sex, sep = "_", treatment))) %>% 
          group_by(treatment, timepoint..month., sex, sex_treatment) %>% 
          summarise(mean(iron.deposition.score..0.4...spleen.)) %>% 
          
          ggplot(aes(x = timepoint..month., y = `mean(iron.deposition.score..0.4...spleen.)`)) +
          geom_point(aes(color = sex_treatment), size = 3, position = position_dodge(0.5)) +
          geom_line(aes(color = sex_treatment)) +
          theme_bw()
          
    
    # liver
        
        # select variables of interest
        df_liver <- df_mice %>% 
          dplyr::select(c(timepoint..month., treatment, sex, iron.deposition.score..0.4...liver.)) %>% 
          dplyr::filter(!is.na(iron.deposition.score..0.4...liver.))
        
        # ordinal regression (select model with lowest AIC)
        ord_liver <- polr(iron.deposition.score..0.4...liver. ~ timepoint..month.*sex + treatment, data = df_liver, Hess = TRUE, method = "probit")
        coef_liver <- ord_liver %>% summary() %>% coef() %>% as_tibble(rownames = "term")
        coef_liver %>% mutate(p.val = pnorm(abs(coef_liver$`t value`), lower.tail = FALSE)*2, OR = exp(Value))
        
        # linear regression
        lm_liver <- lm(iron.deposition.score..0.4...liver. ~ timepoint..month.*treatment, data = df_liver)
        lm_liver %>% summary()
        
        aov_liver <- aov(iron.deposition.score..0.4...liver. ~ timepoint..month.*treatment, data = df_liver)
        aov_liver %>% tidy()
        
          # plot residuals
          augment(lm_liver) %>% 
            ggplot(aes(x = .fitted, y = .resid)) +
            geom_point()
          
          autoplot(lm_liver)
        
        # plot
          df_liver %>% mutate(timepoint..month. = as.numeric(as.character(timepoint..month.)),
                              iron.deposition.score..0.4...liver. = as.numeric(as.character(iron.deposition.score..0.4...liver.)),
                               sex_treatment = factor(paste0(sex, sep = "_", treatment))) %>% 
            group_by(treatment, timepoint..month., sex, sex_treatment) %>% 
            summarise(mean(iron.deposition.score..0.4...liver.)) %>% 
            
            ggplot(aes(x = timepoint..month., y = `mean(iron.deposition.score..0.4...liver.)`)) +
            geom_point(aes(color = sex_treatment), size = 3, position = position_dodge(0.5)) +
            geom_line(aes(color = sex_treatment)) +
            theme_bw()
          
    
    # kidney
        
        # select variables of interest
        df_kidney <- df_mice %>% 
          dplyr::select(c(timepoint..month., treatment, sex, iron.deposition.score..0.4..kidney.)) %>% 
          dplyr::filter(!is.na(iron.deposition.score..0.4..kidney.))
        
        # ordinal regression (select model with lowest AIC)
        ord_kidney <- polr(iron.deposition.score..0.4..kidney. ~ timepoint..month. + treatment + sex, data = df_kidney, Hess = TRUE, method = "probit")
        coef_kidney <- ord_kidney %>% summary() %>% coef() %>% as_tibble(rownames = "term")
        coef_kidney %>% mutate(p.val = pnorm(abs(coef_kidney$`t value`), lower.tail = FALSE)*2, OR = exp(Value))
        
        # linear regression
        lm_kidney <- lm(iron.deposition.score..0.4..kidney. ~ timepoint..month.*treatment, data = df_kidney)
        lm_kidney %>% summary()
        
        aov_kidney <- aov(iron.deposition.score..0.4..kidney. ~ timepoint..month.*treatment, data = df_kidney)
        aov_kidney %>% tidy()
        
          # plot residuals
          augment(lm_kidney) %>% 
            ggplot(aes(x = .fitted, y = .resid)) +
            geom_point()
          
          autoplot(lm_kidney)
          
        # plot
          df_kidney %>% mutate(timepoint..month. = as.numeric(as.character(timepoint..month.)),
                              iron.deposition.score..0.4..kidney. = as.numeric(as.character(iron.deposition.score..0.4..kidney.)),
                              sex_treatment = factor(paste0(sex, sep = "_", treatment))) %>% 
            group_by(treatment, timepoint..month., sex, sex_treatment) %>% 
            summarise(mean(iron.deposition.score..0.4..kidney.)) %>% 
            
            ggplot(aes(x = timepoint..month., y = `mean(iron.deposition.score..0.4..kidney.)`)) +
            geom_point(aes(color = sex_treatment), size = 3, position = position_dodge(0.5)) +
            geom_line(aes(color = sex_treatment)) +
            theme_bw()
          

    
    
# adjust for multiple comparisons ---------------------------------------------------
          
    
          # combine all model results
          df_all_model_rs <- bind_rows(
            
            aov_ferritin %>% tidy() %>% mutate(variable = "ferritin"),
            aov_hamp1 %>% tidy() %>% mutate(variable = "hamp1"),
            aov_hepcidin %>% tidy() %>% mutate(variable = "hepcidin"),
            aov_plasma_heme %>% tidy() %>% mutate(variable = "plasma_heme"),
            aov_haptoglobin %>% tidy() %>% mutate(variable = "haptoglobin"),
            aov_hemopexin %>% tidy() %>% mutate(variable = "hemopexin"),
            aov_p50 %>% tidy() %>% mutate(variable = "p50"),
            aov_round_cell %>% tidy() %>% mutate(variable = "round_cell"),
            aov_rbc %>% tidy() %>% mutate(variable = "rbc"),
            aov_hb %>% tidy() %>% mutate(variable = "Hb"),
            aov_hct %>% tidy() %>% mutate(variable = "hct"),
            aov_mcv %>% tidy() %>% mutate(variable = "mcv"),
            aov_mch %>% tidy() %>% mutate(variable = "mch"),
            aov_mchc %>% tidy() %>% mutate(variable = "mchc"),
            aov_mchc_o %>% tidy() %>% mutate(variable = "mchc_o"),
            aov_rdw %>% tidy() %>% mutate(variable = "rdw"),
            aov_ret %>% tidy() %>% mutate(variable = "ret"),
            aov_arc %>% tidy() %>% mutate(variable = "arc"),
            aov_ret_he %>% tidy() %>% mutate(variable = "ret_he"),
            aov_spleen %>% tidy() %>% mutate(variable = "spleen"),
            aov_liver %>% tidy() %>% mutate(variable = "liver"),
            aov_kidney %>% tidy() %>% mutate(variable = "kidney")
            
          )  
          
          df_all_models_sig_rs <- df_all_model_rs %>% 
            dplyr::filter(!is.na(p.value)) %>% 
            
            # adjust p-values using BH method
            mutate(p.adj = p.adjust(p.value, method = "BH")) %>% 
            
            # find significant comparisons
            filter(p.adj < 0.05)
          
        

# run Tukey HSD on significant results ------------------------------------

          
          f_tukey_rs <- function(aov, var) {
            
            aov %>% 
              TukeyHSD() %>% 
              tidy() %>% 
              filter(term %in% filter(df_all_models_sig_rs, variable == var)$term & adj.p.value < 0.05) %>% 
              print(n = nrow(.))
            
          }
          
          # ferritin
          df_all_models_sig_rs %>% filter(variable == "ferritin")
          f_tukey_rs(aov_ferritin, "ferritin")
          
          # hamp1
          df_all_models_sig_rs %>% filter(variable == "hamp1")
          f_tukey_rs(aov_hamp1, "hamp1")
          
          # hepcidin
          df_all_models_sig_rs %>% filter(variable == "hepcidin")
          f_tukey_rs(aov_hepcidin, "hepcidin")
          
          # plasma heme
          df_all_models_sig_rs %>% filter(variable == "plasma_heme")
          f_tukey_rs(aov_plasma_heme, "plasma_heme")
          
          # haptoglobin
          df_all_models_sig_rs %>% filter(variable == "haptoglobin")
          f_tukey_rs(aov_haptoglobin, "haptoglobin")
          
          # hemopexin
          df_all_models_sig_rs %>% filter(variable == "hemopexin")
          f_tukey_rs(aov_hemopexin, "hemopexin")
          
          # p50
          df_all_models_sig_rs %>% filter(variable == "p50")
          f_tukey_rs(aov_p50, "p50")
          
          # round_cell
          df_all_models_sig_rs %>% filter(variable == "round_cell")
          f_tukey_rs(aov_round_cell, "round_cell")

          # RBC count
          df_all_models_sig_rs %>% filter(variable == "rbc")
          f_tukey_rs(aov_rbc, "rbc")

          # Hb
          df_all_models_sig_rs %>% filter(variable == "Hb")
          f_tukey_rs(aov_hb, "Hb")
          
          # hct          
          df_all_models_sig_rs %>% filter(variable == "hct")
          f_tukey_rs(aov_hct, "hct")
          
          # mcv
          df_all_models_sig_rs %>% filter(variable == "mcv")
          f_tukey_rs(aov_mcv, "mcv")
          
          # mch
          df_all_models_sig_rs %>% filter(variable == "mch")
          f_tukey_rs(aov_mch, "mch")
          
          # mchc
          df_all_models_sig_rs %>% filter(variable == "mchc")
          f_tukey_rs(aov_mchc, "mchc")
          
          # mchc-o
          df_all_models_sig_rs %>% filter(variable == "mchc_o")
          f_tukey_rs(aov_mchc_o, "mchc_o")
          
          # rdw
          df_all_models_sig_rs %>% filter(variable == "rdw")
          f_tukey_rs(aov_rdw, "rdw")
          
          # ret
          df_all_models_sig_rs %>% filter(variable == "ret")
          f_tukey_rs(aov_ret, "ret")
          
          # arc
          df_all_models_sig_rs %>% filter(variable == "arc")
          f_tukey_rs(aov_arc, "arc")
          
          # ret-he
          df_all_models_sig_rs %>% filter(variable == "ret_he")
          f_tukey_rs(aov_ret_he, "ret_he")
          
          # spleen
          df_all_models_sig_rs %>% filter(variable == "spleen")
          f_tukey_rs(aov_spleen, "spleen")
          
          # liver
          df_all_models_sig_rs %>% filter(variable == "liver")
          f_tukey_rs(aov_liver, "liver")
          
          # kidney
          df_all_models_sig_rs %>% filter(variable == "kidney")
          f_tukey_rs(aov_kidney, "kidney")
          
          
          
          
          
          
          
          
          
          
          
          
          
          