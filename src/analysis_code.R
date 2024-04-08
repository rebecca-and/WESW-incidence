library(tidyverse)
library(metafor)
library(INLA)
library(countrycode)

study_data <- read_csv("inputs/study_data.csv")
region <- moz.utils::region()

######################################################################
### Contents:
### Line 12 -Formatting and Meta-analyis (presented in forest plots in figures 1, S1, S2)
### Line 353 - Median IRRs
### Line 371 - Log-Linear Model (SSA, Kenya, ZWE)
### Line 583- Figure 2 Plots (Gives reported correlation between WESW and total population incidence)
######################################################################


######################################################################
####  Calculating standard error and variance for meta-analysis   ####
######################################################################

study_data <- study_data %>%
  mutate(se_log_irr = sqrt(1/new_infections + 1/infections),
         log_irr = log(IRR),
         variance = se_log_irr^2) %>%
  left_join(region)


####################################################
####    Meta-analysis Forest plot formatting    ####
####################################################

# study_data <- study_data %>%
#   arrange(region, iso3, year) %>%
#   mutate(idx = as.numeric(factor(ref))) %>%
#   mutate(slab = paste(year, "|", iso3, "|", ref)) %>%
#   arrange(region, iso3, ref, slab) %>%
#   mutate(region = factor(region),
#          idx = row_number())

custom_match = moz.utils::cc_plot()

names_and_iso3 <- data.frame(name = moz.utils::ssa_names()) %>%
  mutate(iso3 = countrycode(name, "country.name", "iso3c")) %>%
  left_join(data.frame(custom_match) %>% rownames_to_column() %>% rename(iso3 = rowname)) %>%
  mutate(plot_name = ifelse(is.na(custom_match), name, custom_match),
         plot_name = ifelse(plot_name == "United Republic of Tanzania", "Tanzania", plot_name)) %>%
  select(iso3, plot_name)

study_data <- study_data %>%
  left_join(names_and_iso3) %>%
  mutate(slab2 = paste0(year, " | ", plot_name, " | ", ref)) %>%
  arrange(region, plot_name, ref, slab2)  %>%
  mutate(region = factor(region),
         idx = row_number())



#######################################################
####    Meta-analysis + Forest Plot in Figure 1    ####
#######################################################

#### SSA ####

fsw.model <- rma.mv(yi = log_irr,
                    V = variance,
                    slab = slab2,
                    data = study_data,
                    random = ~ 1 | ref/area,
                    test = "t",
                    method = "REML")

##  I^2 calculated using the method described by Viechtbauer in the `metafor` guidance:
##  https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

W <- diag(1/fsw.model$vi)
X <- model.matrix(fsw.model)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

fsw.modelI2 <- round((100 * sum(fsw.model$sigma2) / (sum(fsw.model$sigma2) + (fsw.model$k-fsw.model$p)/sum(diag(P)))), 1)


#### WCA ####

wca.fswmodel <- rma.mv(yi = log_irr,
                       V = variance,
                       slab = slab2,
                       data = study_data %>% filter(region == "WCA"),
                       random = ~ 1 | ref/area,
                       test = "t",
                       method = "REML")

# I^2 #

W <- diag(1/wca.fswmodel$vi)
X <- model.matrix(wca.fswmodel)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

wca.fswmodelI2 <- round((100 * sum(wca.fswmodel$sigma2) / (sum(wca.fswmodel$sigma2) + (wca.fswmodel$k-wca.fswmodel$p)/sum(diag(P)))), 1)


#### ESA ####

esa.fswmodel <- rma.mv(yi = log_irr,
                       V = variance,
                       slab = slab2,
                       data = study_data %>% filter(region == "ESA"),
                       random = ~ 1 | ref/area,
                       test = "t",
                       method = "REML")

# I^2 #

W <- diag(1/esa.fswmodel$vi)
X <- model.matrix(esa.fswmodel)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

esa.fswmodelI2 <- round((100 * sum(esa.fswmodel$sigma2) / (sum(esa.fswmodel$sigma2) + (esa.fswmodel$k-esa.fswmodel$p)/sum(diag(P)))), 1)


### Generate plot

png("figures/Fig 1 Forest.png", width = 600, height = 1000)

forest.rma(fsw.model,
           showweights = F,
           xlim = c(-15,15),
           ylim = c(-3,92),
           atransf = exp,
           header = "Year | Country | Ref",
           cex = 0.75,
           rows=c(3:20, 24:88),
           mlab = paste0("SSA regional estimate   I^2 =", fsw.modelI2, "%"),
           plotwidth=unit(30,"cm"),
           xlab="Incidence Rate Ratio",
           order= -order(study_data$idx))

op <- par(cex=0.8, font=2)
par(font=4)
text(-15, c(89,21), pos=4, c("Eastern and Southern Africa (ESA)",
                             "Western and Central Africa (WCA)"
))
par(op)


addpoly(esa.fswmodel, row=22.5, atransf = exp, mlab = paste0("ESA regional estimate   I^2 =", esa.fswmodelI2, "%"))

addpoly(wca.fswmodel, row=1.5, atransf = exp, mlab = paste0("WCA regional estimate   I^2 =", wca.fswmodelI2, "%"))

dev.off()


########################################################################################
####                  Sensitivity:  Only high quality studies                       ####
####                    Forest Plot in Supplement Figure S2                         ####
########################################################################################

quality <- study_data %>%
  filter(!ref %in% c("Priddy et al", "Van Damme et al", "Laga et al", "Kilburn et al", "Fowke et al", "Kerrigan et al"))



#### SSA ####

qual_fsw.model <- rma.mv(yi = log_irr,
                    V = variance,
                    slab = slab2,
                    data = quality,
                    random = ~ 1 | ref/area,
                    test = "t",
                    method = "REML")

# I^2 #

W <- diag(1/qual_fsw.model$vi)
X <- model.matrix(qual_fsw.model)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

qual_fsw.modelI2 <- round((100 * sum(qual_fsw.model$sigma2) / (sum(qual_fsw.model$sigma2) + (qual_fsw.model$k-qual_fsw.model$p)/sum(diag(P)))), 1)


#### WCA ####

qual_wca.fswmodel <- rma.mv(yi = log_irr,
                       V = variance,
                       slab = slab2,
                       data = quality %>% filter(region == "WCA"),
                       random = ~ 1 | ref/area,
                       test = "t",
                       method = "REML")

# I^2 #

W <- diag(1/qual_wca.fswmodel$vi)
X <- model.matrix(qual_wca.fswmodel)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

qual_wca.fswmodelI2 <- round((100 * sum(qual_wca.fswmodel$sigma2) / (sum(qual_wca.fswmodel$sigma2) + (qual_wca.fswmodel$k-qual_wca.fswmodel$p)/sum(diag(P)))), 1)


#### ESA ####

qual_esa.fswmodel <- rma.mv(yi = log_irr,
                       V = variance,
                       slab = slab2,
                       data = quality %>% filter(region == "ESA"),
                       random = ~ 1 | ref/area,
                       test = "t",
                       method = "REML")

# I^2 #

W <- diag(1/qual_esa.fswmodel$vi)
X <- model.matrix(qual_esa.fswmodel)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

qual_esa.fswmodelI2 <- round((100 * sum(qual_esa.fswmodel$sigma2) / (sum(qual_esa.fswmodel$sigma2) + (qual_esa.fswmodel$k-qual_esa.fswmodel$p)/sum(diag(P)))), 1)


### Generate plot

png("figures/Supp Fig S1 Forest.png", width = 600, height = 1000)

forest.rma(qual_fsw.model,
           showweights = F,
           xlim = c(-15,15),
           ylim = c(-3,83.5),
           atransf = exp,
           header = "Year | Country | Ref",
           cex = 0.75,
           rows=c(3:18, 22:80),
           mlab = paste0("SSA regional estimate   I^2 =", qual_fsw.modelI2, "%"),
           plotwidth=unit(30,"cm"),
           xlab="Incidence Rate Ratio",
           order= -order(quality$idx))

op <- par(cex=0.8, font=2)
par(font=4)
text(-15, c(81,19), pos=4, c("Eastern and Southern Africa (ESA)",
                             "Western and Central Africa (WCA)"
))
par(op)


addpoly(qual_esa.fswmodel, row=20.5, atransf = exp, mlab = paste0("ESA regional estimate   I^2 =", qual_esa.fswmodelI2, "%"))

addpoly(qual_wca.fswmodel, row=1.5, atransf = exp, mlab = paste0("WCA regional estimate   I^2 =", qual_wca.fswmodelI2, "%"))

dev.off()


########################################################################################
####    Sensitivity: Matching WESW to National total female population incidence    ####
####                    Forest Plot in Supplement Figure S2                         ####
########################################################################################


######################################################################
####  Calculating standard error and variance for meta-analysis   ####
######################################################################

nat_study_data <- study_data %>%
  mutate(nat_se_log_irr = sqrt(1/new_infections + 1/nat_infections),
         nat_log_irr = log(nat_IRR),
         nat_variance = nat_se_log_irr^2) %>%
  left_join(region)


#### SSA Model ####

nat_fsw.model <- rma.mv(yi = nat_log_irr,
                    V = variance,
                    slab = slab2,
                    data = nat_study_data,
                    random = ~ 1 | ref/area,
                    test = "t",
                    method = "REML")
# I^2 #

W <- diag(1/nat_fsw.model$vi)
X <- model.matrix(nat_fsw.model)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

nat_fsw.modelI2 <- round((100 * sum(nat_fsw.model$sigma2) / (sum(nat_fsw.model$sigma2) + (nat_fsw.model$k-nat_fsw.model$p)/sum(diag(P)))), 1)


#### WCA ####

nat_wca.fswmodel <- rma.mv(yi = nat_log_irr,
                       V = nat_variance,
                       slab = slab2,
                       data = nat_study_data %>% filter(region == "WCA"),
                       random = ~ 1 | ref/area,
                       test = "t",
                       method = "REML")

# I^2 #

W <- diag(1/nat_wca.fswmodel$vi)
X <- model.matrix(nat_wca.fswmodel)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

nat_wca.fswmodelI2 <- round((100 * sum(nat_wca.fswmodel$sigma2) / (sum(nat_wca.fswmodel$sigma2) + (nat_wca.fswmodel$k-nat_wca.fswmodel$p)/sum(diag(P)))), 1)



#### ESA ####

nat_esa.fswmodel <- rma.mv(yi = nat_log_irr,
                       V = nat_variance,
                       slab = slab2,
                       data = nat_study_data %>% filter(region == "ESA"),
                       random = ~ 1 | ref/area,
                       test = "t",
                       method = "REML")

# I^2 #

W <- diag(1/nat_esa.fswmodel$vi)
X <- model.matrix(nat_esa.fswmodel)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

nat_esa.fswmodelI2 <- round((100 * sum(nat_esa.fswmodel$sigma2) / (sum(nat_esa.fswmodel$sigma2) + (nat_esa.fswmodel$k-nat_esa.fswmodel$p)/sum(diag(P)))), 1)

png("figures/Supp Fig S2 Forest.png", width = 600, height = 1000)

forest.rma(nat_fsw.model,
           showweights = F,
           xlim = c(-15,15),
           ylim = c(-3,92),
           atransf = exp,
           header = "Year | Country | Ref",
           cex = 0.75,
           rows=c(3:20, 24:88),
           mlab = paste0("SSA regional estimate   I^2 =", fsw.modelI2, "%"),
           plotwidth=unit(30,"cm"),
           xlab="Incidence Rate Ratio",
           order= -order(nat_study_data$idx))

op <- par(cex=0.8, font=2)
par(font=4)
text(-15, c(89,21), pos=4, c("Eastern and Southern Africa (ESA)",
                             "Western and Central Africa (WCA)"
))
par(op)


addpoly(nat_esa.fswmodel, row=22.5, atransf = exp, mlab = paste0("ESA regional estimate   I^2 =", nat_esa.fswmodelI2, "%"))

addpoly(nat_wca.fswmodel, row=1.5, atransf = exp, mlab = paste0("WCA regional estimate   I^2 =", nat_wca.fswmodelI2, "%"))

dev.off()

########################################################################################
####                                  Median IRRs                                   ####
####                   Figure 4: IRRs over time & Supp. Table S2                    ####
########################################################################################


esa <- study_data %>%
  filter(region == "ESA")

wca <- study_data %>%
  filter(region == "WCA")

quantile(study_data$IRR)
quantile(esa$IRR)
quantile(wca$IRR)



########################################################################################
####                               Log-linear model                                 ####
####                   Figure 4: IRRs over time & Supp. Table S3                    ####
########################################################################################

# `inla_data` is a dataframe used to run the INLA models.
# type == "timeseries" is used to produce log_irr's over time for all of SSA . It contains the pooled HIV incidence for SSA from 1985-2021. This was calculated using the total number of new infections and population size for women aged 15-39 across the region as predicted by Spectrum.
# type == "country dat" is used for the individual country models (Kenya and Zimbabwe), and contains country specific incidence trajectories for women aged 15-39.
# type == NA is study data (the same as `study_data` above, with some additional indexing columns (id.iso3))
inla_data <- read_csv("inputs/inla_data.csv")


# Full model SSA
ssa_dat <- inla_data %>%
  filter(type == "timeseries" | is.na(type))


ref.iid.prec.prior <- list(prec = list(prior = "normal", param = c(1.6,2)))

formula <- new_infections ~ f(id.ref, model = "iid", hyper = ref.iid.prec.prior) + centre_year + f(id.ref2, centre_year, model = "iid", hyper = ref.iid.prec.prior)

inla_mod <- INLA::inla(formula,
                        data = ssa_dat,
                        family = "Xpoisson",
                        E = ssa_dat$incidence,
                        offset = log(ssa_dat$pys),
                        control.compute = list(config = TRUE,
                                               dic = TRUE),
                        control.predictor=list(compute=TRUE,
                                               link = 1),
                        verbose = FALSE)


summarya <-  summary(inla_mod)

modelsummary2 <- data.frame(summarya$fixed) %>%
  bind_rows(data.frame(summarya$hyperpar))

modelsummary2 <- round(modelsummary2, digits = 3) %>%
  rownames_to_column() %>%
  mutate(`95% CI` = paste0(X0.025quant, ", ",X0.975quant),
         Mean = mean,
         Covariate = rowname) %>%
  select(Covariate, Mean, `95% CI`) %>%
  mutate(Covariate = case_when(Covariate == "(Intercept)" ~ "Intercept (Year = 2003)",
                               Covariate == "centre_year" ~ "Year",
                               Covariate == "Precision for id.ref" ~ "Study Random Intercepts",
                               Covariate == "Precision for id.ref2" ~ "Study Random Slopes with respect to year"
  ))

modelsummary2

filt_df <- ssa_dat %>%
  filter(!is.na(logit_kp_incidence))

samples <- inla.posterior.sample(1000, inla_mod)

contents = inla_mod$misc$configs$contents

effect = "Predictor"
id.effect = which(contents$tag==effect)
ind.effect = contents$start[id.effect]-1 + (1:contents$length[id.effect])

ind.effect <- 1:(nrow(ssa_dat) - nrow(filt_df))

samples.effect = lapply(samples, function(x) x$latent[ind.effect])

incidence_samples <- matrix(sapply(samples.effect, cbind), ncol=1000)

ident <- ssa_dat[ind.effect, ]

qtls <- apply(incidence_samples, 1, quantile, c(0.025, 0.5, 0.975))

dat <- ident %>%
  ungroup() %>%
  mutate(
    lower = qtls[1,],
    median = qtls[2,],
    upper = qtls[3,]
  )

png("figures/Fig 4 IRR over time.png", width = 800, height = 400)

dat %>%
  moz.utils::name_region() %>%
  filter(type == "timeseries") %>%
  ggplot(aes(x = year, y = exp(median))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = ((exp(lower))), ymax = (exp(upper))), alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = 2, colour = "darkred", linewidth = 0.75) +
  geom_point(data = (study_data  %>% moz.utils::name_region() %>% mutate("Person Years" = total_py, Region = region)), aes(y = (IRR), x = year, size = `Person Years`, color = Region), show.legend = TRUE) +
  moz.utils::standard_theme() +
  theme(plot.title = element_text(size = 8)) +
  xlab("Year") +
  theme(axis.title.x = element_blank()) +
  ylab("Incidence Rate Ratio \n(Log Scale)") +
  scale_alpha(guide = "none") +
  scale_y_continuous(trans = "log", breaks = c(0, 1, 10, 100)) +
  scale_colour_manual(values = c("#3B9AB2", "#E1AF00")) +
  theme(aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size=5))) +  # Increase dot size in legend
  theme(legend.text = element_text(size = 11))+
  theme(legend.title = element_text(size = 12, face = "bold"))

dev.off()


# `kenya` and`zwe` are the incidence estimates extracted from studies
kenya <- inla_data %>%
  filter(id.iso3 == 18,
         is.na(type))

zwe <- inla_data %>%
  filter(id.iso3 == 38,
         is.na(type))


datasets2 <- list(

  kenya_dat <- inla_data %>%
    filter(type == "country dat",
           iso3 == "KEN") %>%
    bind_rows(kenya),


  zwe_dat <- inla_data %>%
    filter(type == "country dat",
           iso3 == "ZWE") %>%
    bind_rows(zwe)

)

country_case_study_summaries <- list()
country_case_study_dat <- list()

for (i in 1:length(datasets2)){

  fulldat2a <- datasets2[[i]]


  formula10a <-  new_infections ~ centre_year + f(id.ref, model = "iid", hyper = ref.iid.prec.prior)

  inla_moda <- INLA::inla(formula10a,
                          data = fulldat2a,
                          family = "Xpoisson",
                          E = fulldat2a$incidence,
                          offset = log(fulldat2a$pys),
                          control.compute = list(config = TRUE,
                                                 dic = TRUE),
                          control.predictor=list(compute=TRUE,
                                                 link = 1),
                          verbose = FALSE)


  summarya <-  summary(inla_moda)

  modelsummary2 <- data.frame(summarya$fixed) %>%
    bind_rows(data.frame(summarya$hyperpar))

  modelsummary2 <- round(modelsummary2, digits = 3) %>%
    rownames_to_column() %>%
    mutate(`95% CI` = paste0(X0.025quant, ", ",X0.975quant),
           Mean = mean,
           Covariate = rowname) %>%
    select(Covariate, Mean, `95% CI`) %>%
    mutate(Covariate = case_when(Covariate == "(Intercept)" ~ "Intercept (Year = 2003)",
                                 Covariate == "centre_year" ~ "Year",
                                 Covariate == "Precision for id.ref" ~ "Study Random Intercepts"
    ))

  country_case_study_summaries[[i]] <- modelsummary2

  filt_df <- fulldat2a %>%
    filter(!is.na(logit_kp_incidence))

  samples <- inla.posterior.sample(1000, inla_mod)

  contents = inla_mod$misc$configs$contents

  effect = "Predictor"
  id.effect = which(contents$tag==effect)
  ind.effect = contents$start[id.effect]-1 + (1:contents$length[id.effect])

  ind.effect <- 1:(nrow(fulldat2a) - nrow(filt_df))

  samples.effect = lapply(samples, function(x) x$latent[ind.effect])

  incidence_samples <- matrix(sapply(samples.effect, cbind), ncol=1000)

  ident <- fulldat2a[ind.effect, ]

  qtls <- apply(incidence_samples, 1, quantile, c(0.025, 0.5, 0.975))

  dat <- ident %>%
    ungroup() %>%
    mutate(
      lower = qtls[1,],
      median = qtls[2,],
      upper = qtls[3,]
    )
  country_case_study_dat[[i]] <- dat

}

write.csv(country_case_study_summaries[[1]], "figures/kenya_irr_mod_summary.csv", row.names = F)
write.csv(country_case_study_summaries[[2]], "figures/zwe_irr_mod_summary.csv", row.names = F)






########################################################################################
####           Correlation between WESW and total population incidence              ####
####                          Figure 2: Incidence plots                             ####
########################################################################################


genpop_fsw_incid <- study_data %>%
  moz.utils::name_region() %>%
  mutate(Region = region,
         `Person Years` = total_py,
         incidence = incidence*100,
         incidence_kp = incidence_kp*100) %>%
  ggplot() +
  geom_point(aes(x = exp(log(incidence)), y = exp(log(incidence_kp)), color = Region, size = `Person Years`)) +
  moz.utils::standard_theme() +
  scale_color_manual(values = c("#3B9AB2", "#E1AF00")) +
  lims(x = c(0,4)) +
  geom_abline(linetype = 2, linewidth = 0.75) +
  scale_y_continuous(trans = "log", labels = scales::percent_format(accuracy = 0.01, scale = 1, suffix = ""), breaks = c(0, 0.2, 2, 20)) +
  scale_x_continuous(trans = "log", labels = scales::percent_format(accuracy = 0.01, scale = 1, suffix = ""), breaks = c(0.02, 0.2, 2)) +
  labs( y = "WESW Incidence/100py (Log Scale)", x = "Total Female Population Incidence/100py (Log Scale)", tag = "A")+
  coord_cartesian(y = c(0.1,25)) +
  geom_smooth(aes(x = exp(log(incidence)), y =exp(log(incidence_kp))), method = "lm", formula = y ~ x, linetype = 0, se = FALSE) +
  ggpmisc::stat_poly_eq(aes(x = exp(log(incidence)), y =exp(log(incidence_kp)))) +
  theme(aspect.ratio = 1)+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(legend.text = element_text(size = 11))+
  theme(legend.title = element_text(size = 13, face = "bold"))


natural_incid_time <- study_data %>%
  # group_by(ref) %>%
  moz.utils::name_region() %>%
  mutate(Region = region,
         `Person Years` = total_py,
         incidence = incidence*100,
         incidence_kp = incidence_kp*100) %>%
  # ungroup() %>%
  # left_join(year_seq) %>%
  # group_by(year_group) %>%
  # mutate(median = median(incidence_kp)) %>%
  # ungroup() %>%
  ggplot() +
  geom_point(aes(x = year, y = exp(log(incidence_kp)), color = Region, size = `Person Years`)) +
  lims(x = c(1985, 2020)) +
  moz.utils::standard_theme() +
  scale_color_manual(values = c("#3B9AB2", "#E1AF00")) +
  labs( y = "WESW Incidence/100py", x = "Year", tag = "B") +
  coord_cartesian(y = c(0.1,20))+
  theme(aspect.ratio = 1)+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(legend.text = element_text(size = 11))+
  theme(legend.title = element_text(size = 13, face = "bold"))


incidplots <-  ggpubr::ggarrange(genpop_fsw_incid, natural_incid_time, nrow = 1, common.legend = T, legend = "bottom")

png("figures/Fig 3 IRR plots 2703.png", width = 800, height = 400)

incidplots

dev.off()
