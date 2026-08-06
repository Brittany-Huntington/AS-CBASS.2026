#Code for merging processed site level ED data

#load packages
library(tidyverse)
library(readxl)
library(rstudioapi)
library(RColorBrewer)
library(CBASSED50) #Voolstra R package
library(car) #Leveneʻs test
library(broom)
library(FSA)
library(lme4)
library(broom.mixed)
library(lmerTest)
library(emmeans)
library(purrr)
library(performance)
library(here)
library(multcomp)

rm(list = ls())
custom_shapes <- c(23, 22, 15, 16, 17, 18, 21, 24)
species_colors <- c(
  "AABR" = "#046E8F", 
  "AGLO" = "#83BA75FF", 
  "AHYA" = "#D44D5C",
  "ICRA" = "#462255"  
)

species_labels <- c(
  "AGLO" = "A. globiceps",
  "AABR" = "A. abrotanoides",
  "AHYA" = "A. hyacinthus",
  "ICRA" = "I. crateriformis"
)

#####Load and Prep Data-----------------------------

EDdf <- read.csv(here("Data", "EDsdf_all_ramps.csv")) %>%
  filter(Dataset == "All") %>%
  mutate(across(where(is.character), as.factor))

#add decline width
EDdf <- EDdf %>% mutate(DW = ED95 - ED5)

#supplemental table 
write.csv(
  EDdf,
  paste0("ED_Supplemental_table.csv"),
  row.names = FALSE
)


#long format 
ED_long <- EDdf %>%
  pivot_longer(
    cols = c(ED5, ED50, ED95, DW),
    names_to = "ED_level",
    values_to = "ED_value"
  ) %>%
  mutate(Species = factor(Species, levels = c("AABR", "AHYA", "AGLO", "ICRA")))


#stats: does ED-value differ among species; mixed model for each ED level-------------

#subset by ED
ED5_df  <- filter(ED_long, ED_level == "ED5")
ED50_df <- filter(ED_long, ED_level == "ED50")
ED95_df <- filter(ED_long, ED_level == "ED95")
ED_DW_df <-filter(ED_long, ED_level == "DW")

#evaluate parametric assumptions
#Shapiro-Wilk Normality
ED5_df %>%
  group_by(Species) %>%
  summarise(p_value = shapiro.test(ED_value)$p.value)

ED50_df %>%
  group_by(Species) %>%
  summarise(p_value = shapiro.test(ED_value)$p.value)

ED95_df %>%
  group_by(Species) %>%
  summarise(p_value = shapiro.test(ED_value)$p.value)
#not very normal --> log transform didnʻt help, go non parametric

ED_DW_df%>%
  group_by(Species) %>%
  summarise(p_value = shapiro.test(ED_value)$p.value)

#Homogeneity of variance
leveneTest(ED_value ~ Species, data = ED5_df)
leveneTest(ED_value ~ Species, data = ED50_df)
leveneTest(ED_value ~ Species, data = ED95_df)
leveneTest(ED_value ~ Species, data = ED_DW_df)
#all fine here


models <- ED_long %>%
  split(.$ED_level) %>%
  map(~ lmer(ED_value ~ Species + (1|Site), data = .))

lapply(models, anova)
lapply(models, function(m) emmeans(m, pairwise ~ Species))

# extract ANOVA tables and add ED_level column
anova_table <- lapply(names(models), function(ed) {
  m <- models[[ed]]
  an <- anova(m)
  an_df <- as.data.frame(an) %>%
    rownames_to_column("term") %>%
    mutate(ED_level = ed)
  an_df
}) %>%
  bind_rows() 
# %>%
#   select(-term)

anova_table

#create variance table
var_table <- imap_dfr(models, function(m, ed) {
  vc <- as.data.frame(VarCorr(m))
  site_var <- vc$vcov[vc$grp == "Site"]
  resid_var <- attr(VarCorr(m), "sc")^2
  total <- site_var + resid_var
  tibble(
    ED_level = ed,
    Site_variance = site_var,
    Residual_variance = resid_var,
    Site_prop = site_var / total,
    Residual_prop = resid_var / total
  )
})

var_table

# extract emmeans pairwise comparisons
emmeans_table <- purrr::imap_dfr(models, function(m, ed) {
  em <- emmeans(m, ~ Species)
  tidy(pairs(em)) %>%
    mutate(ED_level = ed)
})

emmeans_table

#check residual diagnostics for ED50 only
m_ED50 <- lmer(ED_value ~ Species + (1|Site), data = ED50_df)
summary(m_ED50)
performance(m_ED50)
anova(m_ED50)
species_emms <- emmeans(m_ED50, ~ Species)
species_pairs <- pairs(species_emms)
print(species_pairs)
species_letters <- cld(species_emms, Letters = letters, adjust = "tukey")
print(species_letters)
# residual diagnostics
par(mfrow = c(1,2))
plot(resid(m_ED50) ~ fitted(m_ED50))
qqnorm(resid(m_ED50))
qqline(resid(m_ED50))



####plotting----------------------------------

#summary stats; mean + SE
ED_summary <- ED_long %>%
  group_by(Species, ED_level) %>%
  summarise(
    n = n(),
    Mean = mean(ED_value),
    SE = sd(ED_value) / sqrt(n()),
    CI_lower = Mean - qt(0.975, df = n - 1) * SE,
    CI_upper = Mean + qt(0.975, df = n - 1) * SE,
    .groups = "drop"
  )

site_summary <- ED_long %>%
  group_by(Site, Species, ED_level) %>%
  summarise(mean_ED = mean(ED_value), .groups="drop")


#boxplot
g1 <- ggplot(ED_long, aes(x = Species, y = ED_value, fill = Species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(aes(color = Species), width = 0.15, alpha = 0.5, size = 1, show.legend = FALSE) +
  facet_wrap(~ ED_level, scales = "free_y") +
  theme_classic(base_size = 14) +
  labs(
    x = "Species",
    y = expression(paste("Temperature (", degree, "C)")),
   # title = "CBASS-Derived Thermal Thresholds",
   # subtitle = "American Samoa"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  theme(legend.position = "none")+
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors)+ 
  scale_x_discrete(labels = species_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))

g1

#Mean 
g2 <- ggplot(ED_long, aes(x = Species, y = ED_value)) +
geom_jitter(aes(color = Species),width = 0.15, alpha = 0.4, size = 1, show.legend = FALSE) +
  
  
  # geom_errorbar(data = ED_summary,
  #               aes(x = Species,
  #                   ymin = Mean - SE,
  #                   ymax = Mean + SE),
  #               width = 0.2,
  #               size = 0.8,
  #               inherit.aes = FALSE) +
  geom_errorbar(
    data = ED_summary,
    aes(
      x = Species,
      ymin = CI_lower,
      ymax = CI_upper
    ),
    width = 0.2,
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  geom_point(data = ED_summary,
             aes(x = Species, y = Mean, fill = Species),
             shape = 21,
             color = "black",
             size = 2,
             stroke = 0.8,
             inherit.aes = FALSE) +
  
  facet_wrap(~ ED_level, scales = "free_y") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Species",
    y = expression(paste("Temperature (", degree, "C)")),
    #title = "CBASS-Derived Thermal Thresholds",
    #subtitle = "American Samoa"
  )+
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors)+ 
  scale_x_discrete(labels = species_labels) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "italic"))

g2


#ggsave(filename = here ("Plots", "All_sites_boxplot.pdf"), plot = g1,  width = 16, height = 9, device = "pdf")
#ggsave(filename = here ("Plots", "All_sites_points&means.pdf"), plot = g2,  width = 16, height = 9, device = "pdf")


#ED50 only mean +/- SE plot ------

ED50_summary <- ED_summary %>% filter(ED_level == "ED50")

cld_df <- data.frame(
  Species = c("AABR", "AHYA", "AGLO", "ICRA"),
  .group  = c("a", "b", "ab", "b")
)

ED50_summary <- ED50_summary %>%
  left_join(cld_df, by = "Species")

g3 <- ggplot(ED50_df, aes(x = Species, y = ED_value)) +
  geom_jitter(aes(color = Species),width = 0.1, alpha = 0.3, size = 2, show.legend = FALSE) +
  geom_point(data = ED50_summary,
             aes(x = Species, y = Mean, color = Species),
             size = 7,
             inherit.aes = FALSE) +
  
  # geom_errorbar(data = ED50_summary,
  #               aes(x = Species,
  #                   ymin = Mean - SE,
  #                   ymax = Mean + SE),
  #               width = 0.2,
  #               size = 1,
  #               inherit.aes = FALSE) +
  geom_errorbar(
    data = ED50_summary,
    aes(
      x = Species,
      ymin = CI_lower,
      ymax = CI_upper
    ),
    width = 0.2,
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  geom_text(data = ED50_summary,
            aes(x = Species,
                y = Mean + SE + 0.1,
                label = .group),
            nudge_x = 0.1,
            hjust = 0,
            size = 7,
            fontface = "bold",
            inherit.aes = FALSE)+
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Species",
    y = expression(paste("ED50 (", degree, "C)")),
   # title = "CBASS-Derived Thermal Thresholds",
  )+
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors)+ 
  scale_x_discrete(labels = species_labels) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "italic"))


g3


#ggsave(filename = here ("Plots", "ED50_post-hocs.pdf"), plot = g3,  width = 16, height = 9, device = "pdf")


#Q1B.  How important are site level differences to our ESA taxa?---------------------------

common_sites <-  c(1,3,4,9,10,11)

dat_esa <- EDdf %>%
  filter(
    Site %in% common_sites,
    Species %in% c("AGLO", "ICRA")
  )

mod_esa <- lmer(ED50 ~ Species + (1 | Site), data = dat_esa)
anova(mod_esa)


vc <- as.data.frame(VarCorr(mod_esa))

site_var <- vc$vcov[vc$grp == "Site"]
resid_var <- attr(VarCorr(mod_esa), "sc")^2

site_prop <- site_var / (site_var + resid_var)
site_prop

library(performance)
r2_nakagawa(mod_esa)
