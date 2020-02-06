# Load packages ---- 
library(tidyverse)
library(broom)


## function to compute student and Welch contrasts

t.contrast <- function(dv, groups, contrast, true_means) {
  means <- by(dv, groups, mean)
  vars <- by(dv, groups, var)
  Ns <- by(dv, groups, length)
  ihat <- contrast %*% means
  
  # student's t test
  df_student <- sum(Ns) - length(Ns)
  mse <- sum(vars * (Ns - 1)) / df_student
  se_student <- sqrt(contrast^2 %*% (mse / Ns))
  t_student <- ihat / se_student
  p_student <- 2 * (1 - pt(abs(t_student), df_student))
  
  # welch's t test
  df_welch <- (sum(vars / Ns))^2 / sum((vars^2 / (Ns^2 * (Ns - 1))))
  se_welch <- sqrt(contrast^2 %*% (vars / Ns))
  t_welch <- ihat / se_welch
  p_welch <- 2 * (1 - pt(abs(t_welch), df_welch))
  
  # coverage probability
  true_ihat <- contrast %*% true_means
  
  t_crit_student <- qt(.025, df = df_student)        
  lb_student <- ihat - t_crit_student * se_student
  ub_student <- ihat + t_crit_student * se_student
  coverage_student <- (ub_student - true_ihat) * (true_ihat - lb_student) > 0
  
  t_crit_welch <- qt(.025, df = df_welch)        
  lb_welch <- ihat - t_crit_welch * se_welch
  ub_welch <- ihat + t_crit_welch * se_welch
  coverage_welch <- (ub_welch - true_ihat) * (true_ihat - lb_welch) > 0
  
  result <- list(p_student = p_student, p_welch = p_welch, coverage_student = coverage_student, coverage_welch = coverage_welch)
  return(result)
}


# Two groups sims ----
nsims <- 10000
set.seed(2184)

two_groups_sims <- tibble(
  
  # specify all the conditions you want 
  # these need to be the same length so it's okay if there are some duplicates within a variable as long as there aren't duplicate combinations of all variables
  # don't worry about matching conditions for now
  sample_ratio = c(1, 3/2, 2, 2, 2),
  min_sample = c(20, 50, 100, 100, 100),
  effect = c(0, .2, .5, .8, .8),
  var_ratio = c(.2, .5, 1, 2, 5)
) %>% 
  
  # complete() will create all possible combinations of the conditions, which is why you don't have to worry about matching above
  complete(sample_ratio, min_sample, effect, var_ratio) %>% 
  
  # now we'll compute the sample sizes, means, and variances for each combination of conditions
  rowwise() %>% 
  mutate(
    ns = map2(min_sample, sample_ratio, ~ c(min_sample, min_sample * sample_ratio)),
    means = map(effect, ~ case_when(
      effect == 0 ~ c(6, 6),
      effect == .2 ~ c(6, 6.28),
      effect == .5 ~ c(6, 6.71),
      effect == .8 ~ c(6, 7.13)
      )
    ),
    vars = map(var_ratio, ~ case_when(
      var_ratio == .2 ~ c(2, 10),
      var_ratio == .5 ~ c(2, 4),
      var_ratio == 1 ~ c(1, 1),
      var_ratio == 2 ~ c(4, 2),
      var_ratio == 5 ~ c(10, 2)
    ))
  ) %>% 
  unnest() %>% 
  
  # group by the conditions and then nest the ns, means and vars
  # this makes it convenient to run the simulations
  group_by(sample_ratio, min_sample, effect, var_ratio) %>% 
  nest() %>% 
  
  # run the simulations
  mutate(
    sim_results = map(data, ~ 
      replicate(nsims, # use replicate to repeat the process for each simulation
        t.contrast(
          dv = c(rnorm(n = .$ns[1], mean = .$means[1], sd = sqrt(.$vars[1])), rnorm(n = .$ns[2], mean = .$means[2], sd = sqrt(.$vars[2]))),
          groups = rep(c(1, 2), times = .$ns),
          contrast = c(-1, 1),
          true_means = .$means
        )
      ) %>% 
        t() %>% # replicate produces a matrix where each row is one output, so transpose to make columns the different outputs
        as.tibble() %>% # coerce to a tibble
        mutate_all(unlist) # unlist the columns of the tibble
    )
  ) %>% 
  
  select(-data) %>% # remove the data so we can unnest
  unnest() %>% 
  
  group_by(sample_ratio, min_sample, effect, var_ratio) %>% # group by conditions and summarize the results
  summarize(
    reject_student = mean(p_student < .05),
    reject_welch = mean(p_welch < .05),
    coverage_student = mean(coverage_student),
    coverage_welch = mean(coverage_welch)
    )







     

# Four groups sims ----
t.multicontrast <- function(dv, groups, true_means, contrast_names, joint_contrasts, ...) {
  
  contrast <- matrix(c(...), nrow = length(list(...)), byrow = TRUE)

  means <- by(dv, groups, mean)
  vars <- by(dv, groups, var)
  Ns <- by(dv, groups, length)
  ihat <- contrast %*% means
  
  # student's t test
  df_student <- sum(Ns) - length(Ns)
  mse <- sum(vars * (Ns - 1)) / df_student
  se_student <- sqrt(mse * (contrast^2 %*% (1 / Ns)))
  t_student <- ihat / se_student
  p_student <- 2 * (1 - pt(abs(t_student), df_student))
  
  # welch's t test
  df_welch <- (sum(vars / Ns))^2 / sum((vars^2 / (Ns^2 * (Ns - 1))))
  se_welch <- sqrt(contrast^2 %*% (vars / Ns))
  t_welch <- ihat / se_welch
  p_welch <- 2 * (1 - pt(abs(t_welch), df_welch))
  
  # coverage probability
  true_ihat <- contrast %*% true_means
  
  t_crit_student <- qt(.025, df = df_student)      
  lb_student <- ihat - t_crit_student * se_student
  ub_student <- ihat + t_crit_student * se_student
  coverage_student <- (ub_student - true_ihat) * (true_ihat - lb_student) > 0
  
  t_crit_welch <- qt(.025, df = df_welch)        
  lb_welch <- ihat - t_crit_welch * se_welch
  ub_welch <- ihat + t_crit_welch * se_welch
  coverage_welch <- (ub_welch - true_ihat) * (true_ihat - lb_welch) > 0
  
  joint_rejects_student = sum(p_student[which(contrast_names %in% joint_contrasts)] < .05) %>% rep(length(list(...))) # rep() just makes other coding convenient by making the the same length as the other outputs but it isn't necessary
  joint_rejects_welch = sum(p_welch[which(contrast_names %in% joint_contrasts)] < .05) %>% rep(length(list(...))) # rep() just makes other coding convenient by making the the same length as the other outputs but it isn't necessary
  
  result <- list(
    contrast_names = contrast_names, 
    p_student = p_student, p_welch = p_welch, 
    coverage_student = coverage_student, coverage_welch = coverage_welch, 
    joint_rejects_student = joint_rejects_student, joint_rejects_welch = joint_rejects_welch
    # ihat = ihat, true_ihat = true_ihat, 
    # se_student = se_student, 
    # df_student = df_student
    )

  return(result)
}




# No effects ====
nsims <- 200
set.seed(2184)

no_effects <- tibble(
  
  # specify all the conditions you want 
  # these need to be the same length so it's okay if there are some duplicates within a variable as long as there aren't duplicate combinations of all variables
  # don't worry about matching conditions for now
  sample_ratio = c(1, 3/2, 2, 2, 2),
  min_sample = c(20, 30, 50, 100, 100),
  effect = c(0),
  var_ratio = c(.2, .5, 1, 2, 5)
) %>% 
  
  # complete() will create all possible combinations of the conditions, which is why you don't have to worry about matching above
  complete(sample_ratio, min_sample, effect, var_ratio) %>% 
  
  # now we'll compute the sample sizes, means, and variances for each combination of conditions
  rowwise() %>% 
  mutate(
    ns = map2(min_sample, sample_ratio, ~ c(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio)),
    means = map(effect, ~ case_when(
      effect == 0 ~ c(6, 6, 6, 6),
      effect == .2 ~ c(6, 6.28, 6.28, 6),
      effect == .5 ~ c(6, 6.71, 6.71, 6),
      effect == .8 ~ c(6, 7.13, 7.13, 6)
    )
    ),
    vars = map(var_ratio, ~ case_when(
      var_ratio == .2 ~ c(2, 10, 10, 10),
      var_ratio == .5 ~ c(2, 4, 4, 4),
      var_ratio == 1 ~ c(2, 2, 2, 2),
      var_ratio == 2 ~ c(4, 2, 2, 2),
      var_ratio == 5 ~ c(10, 2, 2, 2)
    ))
  ) %>% 
  unnest() %>% 
  
  # group by the conditions and then nest the ns, means and vars
  # this makes it convenient to run the simulations
  group_by(sample_ratio, min_sample, effect, var_ratio) %>% 
  nest() %>% 
  
  # run the simulations
  mutate(
    sim_results = map(data, ~ 
                        replicate(nsims, 
                                  t.multicontrast(
                                    dv = c(
                                      rnorm(n = .$ns[1], mean = .$means[1], sd = sqrt(.$vars[1])), 
                                      rnorm(n = .$ns[2], mean = .$means[2], sd = sqrt(.$vars[2])),
                                      rnorm(n = .$ns[3], mean = .$means[3], sd = sqrt(.$vars[3])),
                                      rnorm(n = .$ns[4], mean = .$means[4], sd = sqrt(.$vars[4]))
                                    ),
                                    groups = rep(c(1, 2, 3, 4), times = .$ns),
                                    true_means = .$means,
                                    contrast_names = c(
                                      'ME (1 & 2 vs 3 & 4)', 
                                      'ME (1 & 3 vs 2 & 4)', 
                                      'Interaction', 
                                      'SE (1 vs 2)', 
                                      'SE (3 vs 4)', 
                                      'SE (1 vs 3)', 
                                      'SE (2 vs 4)'
                                      ),
                                    joint_contrasts = c('ME (1 & 2 vs 3 & 4)', 'ME (1 & 3 vs 2 & 4)', 'Interaction'),
                                    c(-1, -1, 1, 1), # main effect 1
                                    c(-1, 1, -1, 1), # main effect 2
                                    c(-1, 1, 1, -1), # interaction effect
                                    c(-1, 1, 0, 0), # simple effect 1 vs 2
                                    c(0, 0, -1, 1), # simple effect 3 vs 4
                                    c(-1, 0, 1, 0), # simple effect 1 vs 3
                                    c(0, -1, 0, 1) # simple effect 2 vs 4
                                  )) %>% 
                        t() %>% 
                        as.tibble() %>% 
                        unnest()
    )) %>% 
  
  select(-data) %>% # remove the data so we can unnest
  unnest() %>% 
  
  group_by(sample_ratio, min_sample, effect, var_ratio) %>% # group by conditions and summarize the results
  group_by(sample_ratio, min_sample, effect, var_ratio, contrast_names) %>% 
  summarize(
    joint_reject_student = mean(joint_rejects_student > 0),
    joint_reject_welch = mean(joint_rejects_welch > 0),
    reject_student = mean(p_student < .05),
    reject_welch = mean(p_welch < .05),
    coverage_student = mean(coverage_student),
    coverage_welch = mean(coverage_welch)
  )

# save simulations
save(no_effects, file = '~/R-projects/Welch rule-ANOVA/no_effects_v2.R')



# false positives (per test)
no_effects %>% 
  ungroup() %>% 
  select(-coverage_student, -coverage_welch, -joint_reject_student, -joint_reject_welch) %>% 
  gather(key = test, value = reject_rate, reject_student, reject_welch) %>% 
  mutate(
    test = case_when(
      test == 'reject_student' ~ "Student's t test",
      test == 'reject_welch' ~ "Welch's t test"
    ),
    contrast_names = factor(contrast_names, levels = c('ME (1 & 2 vs 3 & 4)', 'ME (1 & 3 vs 2 & 4)', 'Interaction')),
    ns = paste(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio, sep = ', ') %>% factor(levels = c('50, 50, 50, 50', '50, 75, 75, 75', '50, 100, 100, 100')),
    var_ratio = factor(var_ratio),
    sample_ratio = factor(sample_ratio)
  ) %>% 
  filter(
    min_sample == 50,
    str_detect(contrast_names, 'ME|Interaction')
  ) %>% 
  ggplot(aes(y = reject_rate, x = var_ratio, group = ns, linetype = ns, shape = ns)) +
  geom_line(size = .5) +
  geom_point(size = 4) + 
  facet_grid(contrast_names ~ test) +
  scale_x_discrete(labels = factor(c('1/5', '1/2', '1', '2', '5'))) +
  scale_y_continuous(breaks = c(0, .05, .10), minor_breaks = seq(0, .13, .01)) +
  labs(
    x = bquote(frac(sigma[1]^2,{sigma[2]^2 == sigma[3]^2} == sigma[4]^2)),
    y = 'False Positive Rate'
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = 'gray80')
  ) +
  scale_shape_manual(name = 'Sample sizes \nGroup 1, 2, 3, 4',
                     labels = c(
                       '50, 50, 50, 50',
                       '50, 75, 75, 75',
                       '50, 100, 100, 100'
                     ),
                     values = rep(c(15, 16, 17), times = 3)
  ) +
  scale_linetype_manual(name = 'Sample sizes \nGroup 1, 2, 3, 4',
                        labels = c(
                          '50, 50, 50, 50',
                          '50, 75, 75, 75',
                          '50, 100, 100, 100'
                        ),
                        values = rep('solid', each = 9)
  ) 


# false positives (per test - including contrasts)
no_effects %>% 
  ungroup() %>% 
  select(-coverage_student, -coverage_welch, -joint_reject_student, -joint_reject_welch) %>% 
  gather(key = test, value = reject_rate, reject_student, reject_welch) %>% 
  mutate(
    test = case_when(
      test == 'reject_student' ~ "Student's t test",
      test == 'reject_welch' ~ "Welch's t test"
    ),
    contrast_names = factor(contrast_names, levels = c('ME (1 & 2 vs 3 & 4)', 'ME (1 & 3 vs 2 & 4)', 'Interaction', 'SE (1 vs 2)', 'SE (3 vs 4)')),
    ns = paste(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio, sep = ', ') %>% factor(levels = c('50, 50, 50, 50', '50, 75, 75, 75', '50, 100, 100, 100')),
    var_ratio = factor(var_ratio),
    sample_ratio = factor(sample_ratio)
  ) %>% 
  filter(
    min_sample == 50,
    #str_detect(contrast_names, 'ME|Interaction')
  ) %>% 
  ggplot(aes(y = reject_rate, x = var_ratio, group = ns, linetype = ns, shape = ns)) +
  geom_line(size = .5) +
  geom_point(size = 4) + 
  facet_grid(contrast_names ~ test) +
  scale_x_discrete(labels = factor(c('1/5', '1/2', '1', '2', '5'))) +
  scale_y_continuous(breaks = c(0, .05, .10), minor_breaks = seq(0, .13, .01)) +
  labs(
    x = bquote(frac(sigma[1]^2,{sigma[2]^2 == sigma[3]^2} == sigma[4]^2)),
    y = 'False Positive Rate'
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = 'gray80')
  ) +
  scale_shape_manual(name = 'Sample sizes \nGroup 1, 2, 3, 4',
                     labels = c(
                       '50, 50, 50, 50',
                       '50, 75, 75, 75',
                       '50, 100, 100, 100'
                     ),
                     values = rep(c(15, 16, 17), times = 3)
  ) +
  scale_linetype_manual(name = 'Sample sizes \nGroup 1, 2, 3, 4',
                        labels = c(
                          '50, 50, 50, 50',
                          '50, 75, 75, 75',
                          '50, 100, 100, 100'
                        ),
                        values = rep('solid', each = 9)
  ) 




# joint false positives (main effects and interaction)
no_effects %>% 
  ungroup() %>% 
  select(-coverage_student, -coverage_welch, -reject_student, -reject_welch) %>% 
  gather(key = test, value = reject_rate, joint_reject_student, joint_reject_welch) %>% 
  mutate(
    test = case_when(
      test == 'joint_reject_student' ~ "Student's t test",
      test == 'joint_reject_welch' ~ "Welch's t test"
    ),
    contrast_names = factor(contrast_names, levels = c('ME (1 & 2 vs 3 & 4)', 'ME (1 & 3 vs 2 & 4)', 'Interaction')),
    ns = paste(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio, sep = ', ') %>% factor(levels = c('50, 50, 50, 50', '50, 75, 75, 75', '50, 100, 100, 100')),
    var_ratio = factor(var_ratio),
    sample_ratio = factor(sample_ratio)
  ) %>% 
  filter(
    min_sample == 50,
    str_detect(contrast_names, 'Interaction')
  ) %>% 
  ggplot(aes(y = reject_rate, x = var_ratio, group = ns, linetype = ns, shape = ns)) +
  geom_line(size = .5) +
  geom_point(size = 4) + 
  facet_grid(contrast_names ~ test) +
  scale_x_discrete(labels = factor(c('1/5', '1/2', '1', '2', '5'))) +
  scale_y_continuous(breaks = seq(0, .30, .05)) +
  labs(
    x = bquote(frac(sigma[1]^2,{sigma[2]^2 == sigma[3]^2} == sigma[4]^2)),
    y = 'False Positive Rate'
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = 'gray80')
  ) +
  scale_shape_manual(name = 'Sample sizes \nGroup 1, 2, 3, 4',
                     labels = c(
                       '50, 50, 50, 50',
                       '50, 75, 75, 75',
                       '50, 100, 100, 100'
                     ),
                     values = rep(c(15, 16, 17), times = 3)
  ) +
  scale_linetype_manual(name = 'Sample sizes \nGroup 1, 2, 3, 4',
                        labels = c(
                          '50, 50, 50, 50',
                          '50, 75, 75, 75',
                          '50, 100, 100, 100'
                        ),
                        values = rep('solid', each = 9)
  ) 









# Crossover interaction, no main effect ====
nsims <- 10000
set.seed(2184)

crossover_sims <- tibble(
  
  # specify all the conditions you want 
  # these need to be the same length so it's okay if there are some duplicates within a variable as long as there aren't duplicate combinations of all variables
  # don't worry about matching conditions for now
  sample_ratio = c(1, 3/2, 2, 2, 2),
  min_sample = c(20, 30, 50, 100, 100),
  effect = c(.5),
  var_ratio = c(.2, .5, 1, 2, 5)
) %>% 
  
  # complete() will create all possible combinations of the conditions, which is why you don't have to worry about matching above
  complete(sample_ratio, min_sample, effect, var_ratio) %>% 
  
  # now we'll compute the sample sizes, means, and variances for each combination of conditions
  rowwise() %>% 
  mutate(
    ns = map2(min_sample, sample_ratio, ~ c(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio)),
    means = map(effect, ~ case_when(
      effect == 0 ~ c(6, 6, 6, 6),
      effect == .2 ~ c(6.28, 6, 6, 6.28),
      effect == .5 ~ c(6.71, 6, 6, 6.71),
      effect == .8 ~ c(7.13, 6, 6, 7.13)
    )
    ),
    vars = map(var_ratio, ~ case_when(
      var_ratio == .2 ~ c(2, 10, 10, 10),
      var_ratio == .5 ~ c(2, 4, 4, 4),
      var_ratio == 1 ~ c(2, 2, 2, 2),
      var_ratio == 2 ~ c(4, 2, 2, 2),
      var_ratio == 5 ~ c(10, 2, 2, 2)
    ))
  ) %>% 
  unnest() %>% 
  
  # group by the conditions and then nest the ns, means and vars
  # this makes it convenient to run the simulations
  group_by(sample_ratio, min_sample, effect, var_ratio) %>% 
  nest() %>% 
  
  # run the simulations
  mutate(
    sim_results = map(data, ~ 
      replicate(nsims, 
        t.multicontrast(
          dv = c(
            rnorm(n = .$ns[1], mean = .$means[1], sd = sqrt(.$vars[1])), 
            rnorm(n = .$ns[2], mean = .$means[2], sd = sqrt(.$vars[2])),
            rnorm(n = .$ns[3], mean = .$means[3], sd = sqrt(.$vars[3])),
            rnorm(n = .$ns[4], mean = .$means[4], sd = sqrt(.$vars[4]))
            ),
          groups = rep(c(1, 2, 3, 4), times = .$ns),
          true_means = .$means,
          contrast_names = c('ME (1 & 2 vs 3 & 4)', 'ME (1 & 3 vs 2 & 4)', 'Interaction', 'SE (1 vs 2)', 'SE (3 vs 4)', 'SE (1 vs 3)', 'SE (2 vs 4)'),
          joint_contrasts = c('ME (1 & 2 vs 3 & 4)', 'ME (1 & 3 vs 2 & 4)'),
          c(-1, -1, 1, 1), # main effect 1
          c(-1, 1, -1, 1), # main effect 2
          c(-1, 1, 1, -1), # interaction effect
          c(-1, 1, 0, 0), # simple effect 1 vs 2
          c(0, 0, -1, 1), # simple effect 3 vs 4
          c(-1, 0, 1, 0), # simple effect 1 vs 3
          c(0, -1, 0, 1) # simple effect 2 vs 4
        )) %>%
        t() %>% 
        as.tibble() %>% 
        unnest()
      )) %>% 
  
  select(-data) %>% # remove the data so we can unnest
  unnest() %>% 
  
  group_by(sample_ratio, min_sample, effect, var_ratio, contrast_names) %>% # group by conditions and summarize the results
  summarize(
    joint_reject_student = mean(joint_rejects_student > 1),
    joint_reject_welch = mean(joint_rejects_welch > 1),
    reject_student = mean(p_student < .05),
    reject_welch = mean(p_welch < .05),
    coverage_student = mean(coverage_student),
    coverage_welch = mean(coverage_welch)
  )


# save simulations
save(crossover_sims, file = '~/R-projects/Welch rule-ANOVA/crossover_sims.R')



# false positives
# joint false positives (main effects and interaction)
crossover_sims %>% 
  ungroup() %>% 
  select(-coverage_student, -coverage_welch) %>% 
  gather(key = test, value = reject_rate, reject_student, reject_welch) %>% 
  mutate(
    test = case_when(
      test == 'reject_student' ~ "Student's t test",
      test == 'reject_welch' ~ "Welch's t test"
    ),
    contrast_names = factor(contrast_names, levels = c('ME (1 & 2 vs 3 & 4)', 'ME (1 & 3 vs 2 & 4)')),
    ns = paste(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio, sep = ', ') %>% factor(levels = c('50, 50, 50, 50', '50, 75, 75, 75', '50, 100, 100, 100')),
    var_ratio = factor(var_ratio),
    sample_ratio = factor(sample_ratio)
  ) %>% 
  filter(
    min_sample == 50,
    str_detect(contrast_names, 'ME')
  ) %>% 
  ggplot(aes(y = reject_rate, x = var_ratio, group = ns, linetype = ns, shape = ns)) +
  geom_line(size = .5) +
  geom_point(size = 4) + 
  facet_grid(contrast_names ~ test) +
  scale_x_discrete(labels = factor(c('1/5', '1/2', '1', '2', '5'))) +
  scale_y_continuous(breaks = seq(0, .30, .05)) +
  labs(
    x = bquote(frac(sigma[1]^2,{sigma[2]^2 == sigma[3]^2} == sigma[4]^2)),
    y = 'False Positive Rate'
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = 'gray80')
  ) +
  scale_shape_manual(name = 'Sample sizes \nGroup 1, 2, 3, 4',
                     labels = c(
                       '50, 50, 50, 50',
                       '50, 75, 75, 75',
                       '50, 100, 100, 100'
                     ),
                     values = rep(c(15, 16, 17), times = 3)
  ) +
  scale_linetype_manual(name = 'Sample sizes \nGroup 1, 2, 3, 4',
                        labels = c(
                          '50, 50, 50, 50',
                          '50, 75, 75, 75',
                          '50, 100, 100, 100'
                        ),
                        values = rep('solid', each = 9)
  ) 


# power
crossover_sims %>% 
  ungroup() %>% 
  select(-coverage_student, -coverage_welch) %>% 
  gather(key = test, value = reject_rate, reject_student, reject_welch) %>% 
  mutate(
    test = case_when(
      test == 'reject_student' ~ "Student's t test",
      test == 'reject_welch' ~ "Welch's t test"
    ),
    contrast_names = factor(contrast_names, levels = c('Interaction', 'SE (1 vs 2)', 'SE (3 vs 4)')),
    ns = paste(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio, sep = ', ') %>% factor(levels = c('50, 50, 50, 50', '50, 75, 75, 75', '50, 100, 100, 100')),
    var_ratio = factor(var_ratio),
    sample_ratio = factor(sample_ratio)
  ) %>% 
  filter(
    min_sample == 50,
    str_detect(contrast_names, 'Int|SE')
  ) %>% 
  ggplot(aes(y = reject_rate, x = var_ratio, group = ns, linetype = ns, shape = ns)) +
  geom_line(size = .5) +
  geom_point(size = 4) + 
  facet_grid(contrast_names ~ test) +
  scale_x_discrete(labels = factor(c('1/5', '1/2', '1', '2', '5'))) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  labs(
    x = bquote(frac(sigma[1]^2,{sigma[2]^2 == sigma[3]^2} == sigma[4]^2)),
    y = 'Power'
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = 'gray80')
  ) +
  scale_shape_manual(name = 'Sample sizes \nGroup 1, 2, 3, 4',
                     labels = c(
                       '50, 50, 50, 50',
                       '50, 75, 75, 75',
                       '50, 100, 100, 100'
                     ),
                     values = rep(c(15, 16, 17), times = 3)
  ) +
  scale_linetype_manual(name = 'Sample sizes \nGroup 1, 2, 3, 4',
                        labels = c(
                          '50, 50, 50, 50',
                          '50, 75, 75, 75',
                          '50, 100, 100, 100'
                        ),
                        values = rep('solid', each = 9)
  ) + 
  coord_cartesian(ylim = c(0, 1))




# difference in power
crossover_sims %>% 
  ungroup() %>% 
  select(-coverage_student, -coverage_welch) %>% 
  mutate(diff = reject_welch - reject_student) %>% 
  gather(key = test, value = reject_rate, reject_student, reject_welch, diff) %>% 
  mutate(
    test = case_when(
      test == 'reject_student' ~ "Student's t test",
      test == 'reject_welch' ~ "Welch's t test",
      test == 'diff' ~ 'Welch - Student'
    ),
    ns = paste(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio, sep = ', ') %>% factor(),
    var_ratio = factor(var_ratio),
    sample_ratio = factor(sample_ratio)
  ) %>% 
  filter(
    test == 'Welch - Student',
    min_sample == 50,
    !str_detect(contrast_names, 'ME')
  ) %>% 
  ggplot(aes(y = reject_rate, x = var_ratio, group = ns)) +
  geom_line() +
  facet_grid(contrast_names ~ test) +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'gray80')

# coverage
crossover_sims %>% 
  ungroup() %>% 
  select(-reject_student, -reject_welch) %>% 
  mutate(diff = coverage_welch - coverage_student) %>% 
  gather(key = test, value = coverage, coverage_welch, coverage_student, diff) %>% 
  mutate(
    test = case_when(
      test == 'coverage_student' ~ "Student's t test",
      test == 'coverage_welch' ~ "Welch's t test",
      test == 'diff' ~ 'Welch - Student'
    ),
    ns = paste(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio, sep = ', ') %>% factor(),
    var_ratio = factor(var_ratio),
    sample_ratio = factor(sample_ratio)
  ) %>% 
  filter(
    test != 'Welch - Student',
    min_sample == 50
  ) %>% 
  ggplot(aes(y = coverage, x = var_ratio, group = ns)) +
  geom_hline(aes(yintercept = .95), color = 'gray80') +
  geom_line(aes(group = ns)) +
  facet_grid(contrast_names ~ test) +
  theme_classic()




# Main effect + interaction ====
nsims <- 10000
set.seed(2184)

me_int_sims <- tibble(
  
  # specify all the conditions you want 
  # these need to be the same length so it's okay if there are some duplicates within a variable as long as there aren't duplicate combinations of all variables
  # don't worry about matching conditions for now
  sample_ratio = c(1, 3/2, 2, 2, 2),
  min_sample = c(20, 30, 50, 100, 100),
  effect = c(.5),
  var_ratio = c(.2, .5, 1, 2, 5)
) %>% 
  
  # complete() will create all possible combinations of the conditions, which is why you don't have to worry about matching above
  complete(sample_ratio, min_sample, effect, var_ratio) %>% 
  
  # now we'll compute the sample sizes, means, and variances for each combination of conditions
  rowwise() %>% 
  mutate(
    ns = map2(min_sample, sample_ratio, ~ c(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio)),
    means = map(effect, ~ case_when(
      effect == .5 ~ c(6.71, 6, 6.28, 6)
    )
    ),
    vars = map(var_ratio, ~ case_when(
      var_ratio == .2 ~ c(2, 10, 10, 10),
      var_ratio == .5 ~ c(2, 4, 4, 4),
      var_ratio == 1 ~ c(2, 2, 2, 2),
      var_ratio == 2 ~ c(4, 2, 2, 2),
      var_ratio == 5 ~ c(10, 2, 2, 2)
    ))
  ) %>% 
  unnest() %>% 
  
  # group by the conditions and then nest the ns, means and vars
  # this makes it convenient to run the simulations
  group_by(sample_ratio, min_sample, effect, var_ratio) %>% 
  nest() %>% 
  
  # run the simulations
  mutate(
    sim_results = map(data, ~ 
                        replicate(nsims, 
                                  t.multicontrast(
                                    dv = c(
                                      rnorm(n = .$ns[1], mean = .$means[1], sd = sqrt(.$vars[1])), 
                                      rnorm(n = .$ns[2], mean = .$means[2], sd = sqrt(.$vars[2])),
                                      rnorm(n = .$ns[3], mean = .$means[3], sd = sqrt(.$vars[3])),
                                      rnorm(n = .$ns[4], mean = .$means[4], sd = sqrt(.$vars[4]))
                                    ),
                                    groups = rep(c(1, 2, 3, 4), times = .$ns),
                                    true_means = .$means,
                                    contrast_names = c('ME (1 & 2 vs 3 & 4)', 'ME (1 & 3 vs 2 & 4)', 'Interaction', 'SE (1 vs 2)', 'SE (3 vs 4)', 'SE (1 vs 3)', 'SE (2 vs 4)'),
                                    joint_contrasts = c('ME (1 & 2 vs 3 & 4)', 'ME (1 & 3 vs 2 & 4)', 'Interaction'),
                                    c(-1, -1, 1, 1), # main effect 1
                                    c(-1, 1, -1, 1), # main effect 2
                                    c(-1, 1, 1, -1), # interaction effect
                                    c(-1, 1, 0, 0), # simple effect 1 vs 2
                                    c(0, 0, -1, 1), # simple effect 3 vs 4
                                    c(-1, 0, 1, 0), # simple effect 1 vs 3
                                    c(0, -1, 0, 1) # simple effect 2 vs 4
                                  )) %>%
                        t() %>% 
                        as.tibble() %>% 
                        unnest()
    )) %>% 
  
  select(-data) %>% # remove the data so we can unnest
  unnest() %>% 
  
  group_by(sample_ratio, min_sample, effect, var_ratio, contrast_names) %>% # group by conditions and summarize the results
  summarize(
    joint_reject_student = mean(joint_rejects_student > 1),
    joint_reject_welch = mean(joint_rejects_welch > 1),
    reject_student = mean(p_student < .05),
    reject_welch = mean(p_welch < .05),
    coverage_student = mean(coverage_student),
    coverage_welch = mean(coverage_welch)
  )

save(me_int_sims, file = '~/R-projects/Welch rule-ANOVA/me_int_sims.R')


# power of tests
me_int_sims %>% 
  ungroup() %>% 
  select(-coverage_student, -coverage_welch) %>% 
  mutate(diff = reject_welch - reject_student) %>% 
  gather(key = test, value = reject_rate, reject_student, reject_welch, diff) %>% 
  mutate(
    test = case_when(
      test == 'reject_student' ~ "Student's t test",
      test == 'reject_welch' ~ "Welch's t test",
      test == 'diff' ~ 'Welch - Student'
    ),
    ns = paste(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio, sep = ', ') %>% factor(),
    var_ratio = factor(var_ratio),
    sample_ratio = factor(sample_ratio)
  ) %>% 
  filter(
    test != 'Welch - Student',
    min_sample == 50
  ) %>% 
  ggplot(aes(y = reject_rate, x = var_ratio, group = ns, color = reject_rate*4)) +
  geom_line() +
  facet_grid(contrast_names ~ test) +
  theme_classic() +
  scale_color_gradient2(low = 'darkred', high = 'darkblue', mid = 'gray80') +
  geom_hline(aes(yintercept = 0), color = 'gray80')


# difference in power
me_int_sims %>% 
  ungroup() %>% 
  select(-coverage_student, -coverage_welch) %>% 
  mutate(diff = reject_welch - reject_student) %>% 
  gather(key = test, value = reject_rate, reject_student, reject_welch, diff) %>% 
  mutate(
    test = case_when(
      test == 'reject_student' ~ "Student's t test",
      test == 'reject_welch' ~ "Welch's t test",
      test == 'diff' ~ 'Welch - Student'
    ),
    ns = paste(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio, sep = ', ') %>% factor(),
    var_ratio = factor(var_ratio),
    sample_ratio = factor(sample_ratio)
  ) %>% 
  filter(
    test == 'Welch - Student',
    min_sample == 50
    ) %>% 
  ggplot(aes(y = reject_rate, x = var_ratio, group = ns)) +
  geom_line() +
  facet_grid(contrast_names ~ test) +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'gray80')

# coverage
me_int_sims %>% 
  ungroup() %>% 
  select(-reject_student, -reject_welch) %>% 
  mutate(diff = coverage_welch - coverage_student) %>% 
  gather(key = test, value = coverage, coverage_welch, coverage_student, diff) %>% 
  mutate(
    test = case_when(
      test == 'coverage_student' ~ "Student's t test",
      test == 'coverage_welch' ~ "Welch's t test",
      test == 'diff' ~ 'Welch - Student'
    ),
    ns = paste(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio, sep = ', ') %>% factor(),
    var_ratio = factor(var_ratio),
    sample_ratio = factor(sample_ratio)
  ) %>% 
  filter(
    test != 'Welch - Student',
    min_sample == 50
  ) %>% 
  ggplot(aes(y = coverage, x = var_ratio, group = ns)) +
  geom_hline(aes(yintercept = .95), color = 'gray80') +
  geom_line() +
  facet_grid(contrast_names ~ test) +
  theme_classic()






# show joint probability of a type 1 error



# Main effect + interaction - new variance ====
nsims <- 10000
set.seed(2184)

me_int_sims_test <- tibble(
  
  # specify all the conditions you want 
  # these need to be the same length so it's okay if there are some duplicates within a variable as long as there aren't duplicate combinations of all variables
  # don't worry about matching conditions for now
  sample_ratio = c(1, 3/2, 2, 2, 2),
  min_sample = c(20, 30, 50, 100, 100),
  effect = c(.5),
  var_ratio = c(.2, .5, 1, 2, 5)
) %>% 
  
  # complete() will create all possible combinations of the conditions, which is why you don't have to worry about matching above
  complete(sample_ratio, min_sample, effect, var_ratio) %>% 
  
  # now we'll compute the sample sizes, means, and variances for each combination of conditions
  rowwise() %>% 
  mutate(
    ns = map2(min_sample, sample_ratio, ~ c(min_sample * sample_ratio, min_sample, min_sample, min_sample)),
    means = map(effect, ~ case_when(
      effect == .5 ~ c(6.71, 6, 6.28, 6)
    )
    ),
    vars = map(var_ratio, ~ case_when(
      var_ratio == .2 ~ c(.4, 2, 2, 2),
      var_ratio == .5 ~ c(1, 2, 2, 2),
      var_ratio == 1 ~ c(2, 2, 2, 2),
      var_ratio == 2 ~ c(2, 1, 1, 1),
      var_ratio == 5 ~ c(2, .4, .4, .4)
    ))
  ) %>% 
  unnest() %>% 
  
  # group by the conditions and then nest the ns, means and vars
  # this makes it convenient to run the simulations
  group_by(sample_ratio, min_sample, effect, var_ratio) %>% 
  nest() %>% 
  
  # run the simulations
  mutate(
    sim_results = map(data, ~ 
                        replicate(nsims, 
                                  t.multicontrast(
                                    dv = c(
                                      rnorm(n = .$ns[1], mean = .$means[1], sd = sqrt(.$vars[1])), 
                                      rnorm(n = .$ns[2], mean = .$means[2], sd = sqrt(.$vars[2])),
                                      rnorm(n = .$ns[3], mean = .$means[3], sd = sqrt(.$vars[3])),
                                      rnorm(n = .$ns[4], mean = .$means[4], sd = sqrt(.$vars[4]))
                                    ),
                                    groups = rep(c(1, 2, 3, 4), times = .$ns),
                                    true_means = .$means,
                                    contrast_names = c('ME (1 & 2 vs 3 & 4)', 'ME (1 & 3 vs 2 & 4)', 'Interaction', 'SE (1 vs 2)', 'SE (3 vs 4)', 'SE (1 vs 3)', 'SE (2 vs 4)'),
                                    joint_contrasts = c('ME (1 & 2 vs 3 & 4)', 'ME (1 & 3 vs 2 & 4)', 'Interaction'),
                                    c(-1, -1, 1, 1), # main effect 1
                                    c(-1, 1, -1, 1), # main effect 2
                                    c(-1, 1, 1, -1), # interaction effect
                                    c(-1, 1, 0, 0), # simple effect 1 vs 2
                                    c(0, 0, -1, 1), # simple effect 3 vs 4
                                    c(-1, 0, 1, 0), # simple effect 1 vs 3
                                    c(0, -1, 0, 1) # simple effect 2 vs 4
                                  )) %>%
                        t() %>% 
                        as.tibble() %>% 
                        unnest()
    )) %>% 
  
  select(-data) %>% # remove the data so we can unnest
  unnest() %>% 
  group_by(sample_ratio, min_sample, effect, var_ratio, contrast_names) %>% 
  summarize(
    joint_rejects_student = mean(joint_rejects_student >= 1),
    joint_rejects_welch = mean(joint_rejects_welch >= 1),
    reject_student = mean(p_student < .05),
    reject_welch = mean(p_welch < .05),
    coverage_student = mean(coverage_student),
    coverage_welch = mean(coverage_welch)
  )

save(me_int_sims_test, file = '~/R-projects/Welch rule-ANOVA/me_int_sims_test.R')


# power of tests
me_int_sims_test %>% 
  ungroup() %>% 
  select(-coverage_student, -coverage_welch) %>% 
  mutate(diff = reject_welch - reject_student) %>% 
  gather(key = test, value = reject_rate, reject_student, reject_welch, diff) %>% 
  mutate(
    test = case_when(
      test == 'reject_student' ~ "Student's t test",
      test == 'reject_welch' ~ "Welch's t test",
      test == 'diff' ~ 'Welch - Student'
    ),
    ns = paste(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio, sep = ', ') %>% factor(),
    var_ratio = factor(var_ratio),
    sample_ratio = factor(sample_ratio)
  ) %>% 
  filter(
    test != 'Welch - Student',
    min_sample == 50
  ) %>% 
  ggplot(aes(y = reject_rate, x = var_ratio, group = ns, color = reject_rate*4)) +
  geom_line() +
  facet_grid(contrast_names ~ test) +
  theme_classic() +
  scale_color_gradient2(low = 'darkred', high = 'darkblue', mid = 'gray80') +
  geom_hline(aes(yintercept = 0), color = 'gray80')


# difference in power
me_int_sims_test %>% 
  ungroup() %>% 
  select(-coverage_student, -coverage_welch) %>% 
  mutate(diff = reject_welch - reject_student) %>% 
  gather(key = test, value = reject_rate, reject_student, reject_welch, diff) %>% 
  mutate(
    test = case_when(
      test == 'reject_student' ~ "Student's t test",
      test == 'reject_welch' ~ "Welch's t test",
      test == 'diff' ~ 'Welch - Student'
    ),
    ns = paste(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio, sep = ', ') %>% factor(),
    var_ratio = factor(var_ratio),
    sample_ratio = factor(sample_ratio)
  ) %>% 
  filter(
    test == 'Welch - Student',
    min_sample == 50
  ) %>% 
  mutate(ns = factor(ns)) %>% 
  ggplot(aes(y = reject_rate, x = var_ratio, group = ns, color = ns)) +
  geom_line() +
  facet_grid(contrast_names ~ test) +
  theme_classic() +
  geom_hline(aes(yintercept = 0), color = 'gray80')

# coverage
me_int_sims_test %>% 
  ungroup() %>% 
  select(-reject_student, -reject_welch) %>% 
  mutate(diff = coverage_welch - coverage_student) %>% 
  gather(key = test, value = coverage, coverage_welch, coverage_student, diff) %>% 
  mutate(
    test = case_when(
      test == 'coverage_student' ~ "Student's t test",
      test == 'coverage_welch' ~ "Welch's t test",
      test == 'diff' ~ 'Welch - Student'
    ),
    ns = paste(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio, sep = ', ') %>% factor(),
    var_ratio = factor(var_ratio),
    sample_ratio = factor(sample_ratio)
  ) %>% 
  filter(
    test != 'Welch - Student',
    min_sample == 50
  ) %>% 
  ggplot(aes(y = coverage, x = var_ratio, group = ns)) +
  geom_hline(aes(yintercept = .95), color = 'gray80') +
  geom_line() +
  facet_grid(contrast_names ~ test) +
  theme_classic()


# Understand contrasts ----
vars1 <- c(2, 10, 10, 10)
vars2 <- c(.4, 2, 2, 2)
ns <- rep(50, 4)

pooled_var1 <- sum(vars1 * (ns - 1))/(sum(ns - 1))

sqrt(pooled_var1/50 + pooled_var1/50)

sqrt(10/50 + 10/50)


pooled_var2 <- sum(vars2 * (ns - 1))/(sum(ns - 1))

sqrt(pooled_var2/50 + pooled_var2/50)

sqrt(2/50 + 2/50)


vars1 <- c(10, 2, 2, 2)
vars2 <- c(.4, 2, 2, 2)
ns <- rep(50, 4)


vars1 <- c(10, 2, 2, 2)
vars2 <- c(.4, 2, 2, 2)
ns <- c(50, 50, 50, 50)

pooled_var2 <- sum(vars2 * (ns - 1))/(sum(ns - 1))

sqrt(pooled_var2/50 + pooled_var2/50)

# check
sqrt(vars2[1]/(2*50) + sum(vars2[2:4])/(2*50))

sqrt(vars2[3]/50 + vars2[4]/50)


# Similarly, the results for the contrasts show that Welchs t test is sometimes more powerful than Students, depending on different conditions. However, the contrast break some of the other patterns for the following reasons.
# 
# 1) Welchs and Students t tests use different groups to find the standard error. Students t test uses the variance from all four groups to estimate the common pooled variance, which is then used to estimate the standard error. Welchs t test doesn't estimate a pooled variance, and it only uses the groups included in the contrast to estimate the standard error. So Students t test can be affected by unequal variances even when the variances of the focal groups in the contrast are equal.
# 
# 2) Sample sizes are weighted differently. When there are more than two groups, Students t test weighs each group based on its sample size to estimate the pooled variance, then weighs only the groups involved the contrast based on their sample size to estimate the standard error. Welchs t test only weights the sample sizes of the groups involved in the contrast. 
# 
# In contrasts between two groups, Students t test computes the standard error by first using the variances from all four groups to find the pooled variance. The standard error will be smaller if one group has a smaller variance than the others, and it will be larger if one group has a larger variance than the others. 
# 
# from all four groups is used to find the standard error, but 
# 
# takes into account the variances from all of the groups, 










# No effects test ==== 
nsims <- 500
set.seed(2184)

no_effects <- tibble(

# specify all the conditions you want 
# these need to be the same length so it's okay if there are some duplicates within a variable as long as there aren't duplicate combinations of all variables
# don't worry about matching conditions for now
sample_ratio = c(1),
min_sample = c(30),
effect = c(0),
var_ratio = c(.2)
) %>% 
  
  # complete() will create all possible combinations of the conditions, which is why you don't have to worry about matching above
  complete(sample_ratio, min_sample, effect, var_ratio) %>% 
  
  # now we'll compute the sample sizes, means, and variances for each combination of conditions
  rowwise() %>% 
  mutate(
    ns = map2(min_sample, sample_ratio, ~ c(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio)),
    means = map(effect, ~ case_when(
      effect == 0 ~ c(6, 6, 6, 6),
      effect == .2 ~ c(6, 6.28, 6.28, 6),
      effect == .5 ~ c(6, 6.71, 6.71, 6),
      effect == .8 ~ c(6, 7.13, 7.13, 6)
    )
    ),
    vars = map(var_ratio, ~ case_when(
      var_ratio == .2 ~ c(2, 10, 10, 10),
      var_ratio == .5 ~ c(2, 4, 4, 4),
      var_ratio == 1 ~ c(2, 2, 2, 2),
      var_ratio == 2 ~ c(4, 2, 2, 2),
      var_ratio == 5 ~ c(10, 2, 2, 2)
    ))
  ) %>% 
  unnest() %>% 
  
  # group by the conditions and then nest the ns, means and vars
  # this makes it convenient to run the simulations
  group_by(sample_ratio, min_sample, effect, var_ratio) %>% 
  nest() 

mydata <- no_effects %>% 
  select(data) %>% 
  unnest()

for(i in 1:length(mydata$vars)) {
  temp = rnorm(
    n = mydata$ns[i],
    mean = mydata$means[i],
    sd = sqrt(mydata$vars[i])
  )
  if (i == 1){
    dv = temp
  } else {
    dv = c(dv, temp)
  }
}
groups <- rep(c(1, 2, 3, 4), times = mydata$ns)
true_means <- mydata$means
contrast_names = c(
  'SE (1 vs 3)', 
  'SE (2 vs 4)'
)
joint_contrasts = c('')

t.multicontrast(dv, groups, true_means, contrast_names, joint_contrasts, c(1, -1, 0, 0), c(0, 0, 1, -1), c(1, -1, 1, -1))




%>% 
  
  # run the simulations
  mutate(
    sim_results = map(data, ~ 
                        replicate(nsims, 
                                  t.multicontrast(
                                    dv = c(
                                      rnorm(n = .$ns[1], mean = .$means[1], sd = sqrt(.$vars[1])), 
                                      rnorm(n = .$ns[2], mean = .$means[2], sd = sqrt(.$vars[2])),
                                      rnorm(n = .$ns[3], mean = .$means[3], sd = sqrt(.$vars[3])),
                                      rnorm(n = .$ns[4], mean = .$means[4], sd = sqrt(.$vars[4]))
                                    ),
                                    groups = rep(c(1, 2, 3, 4), times = .$ns),
                                    true_means = .$means,
                                    contrast_names = c(
                                      'SE (1 vs 3)', 
                                      'SE (2 vs 4)'
                                    ),
                                    joint_contrasts = c(''),
                                    c(-1, 0, 1, 0), # simple effect 1 vs 3
                                    c(0, -1, 0, 1) # simple effect 2 vs 4
                                  )) %>% 
                        t() %>% 
                        as.tibble() %>% 
                        unnest()
    )) %>% 
  select(-data) %>% # remove the data so we can unnest
  unnest()

no_effects %>% 
  group_by(contrast_names) %>% 
  summarize(
    mean(ihat),
    median(ihat),
    max(ihat),
    mean(true_ihat),
    median(true_ihat),
    max(true_ihat),
    mean(se_student),
    median(se_student)
  )

no_effects %>% 
  group_by(contrast_names) %>% 
  summarize(r = cor(ihat, se_student))
  
  
nsims <- 500
set.seed(2184)
no_effects <- tibble(
  
  # specify all the conditions you want 
  # these need to be the same length so it's okay if there are some duplicates within a variable as long as there aren't duplicate combinations of all variables
  # don't worry about matching conditions for now
  sample_ratio = c(1),
  min_sample = c(30),
  effect = c(0),
  var_ratio = c(.2)
) %>% 
  
  # complete() will create all possible combinations of the conditions, which is why you don't have to worry about matching above
  complete(sample_ratio, min_sample, effect, var_ratio) %>% 
  
  # now we'll compute the sample sizes, means, and variances for each combination of conditions
  rowwise() %>% 
  mutate(
    ns = map2(min_sample, sample_ratio, ~ c(min_sample, min_sample * sample_ratio, min_sample * sample_ratio, min_sample * sample_ratio)),
    means = map(effect, ~ case_when(
      effect == 0 ~ c(6, 6, 6, 6),
      effect == .2 ~ c(6, 6.28, 6.28, 6),
      effect == .5 ~ c(6, 6.71, 6.71, 6),
      effect == .8 ~ c(6, 7.13, 7.13, 6)
    )
    ),
    vars = map(var_ratio, ~ case_when(
      var_ratio == .2 ~ c(2, 10, 10, 10),
      var_ratio == .5 ~ c(2, 4, 4, 4),
      var_ratio == 1 ~ c(2, 2, 2, 2),
      var_ratio == 2 ~ c(4, 2, 2, 2),
      var_ratio == 5 ~ c(10, 2, 2, 2)
    ))
  ) %>% 
  unnest() %>% 
  
  # group by the conditions and then nest the ns, means and vars
  # this makes it convenient to run the simulations
  group_by(sample_ratio, min_sample, effect, var_ratio) %>% 
  nest() %>% 
  
  # run the simulations
  mutate(
    sim_results = map(data, ~ 
                        replicate(nsims, 
                                  t.multicontrast(
                                    dv = c(
                                      rnorm(n = .$ns[1], mean = .$means[1], sd = sqrt(.$vars[1])), 
                                      rnorm(n = .$ns[2], mean = .$means[2], sd = sqrt(.$vars[2])),
                                      rnorm(n = .$ns[3], mean = .$means[3], sd = sqrt(.$vars[3])),
                                      rnorm(n = .$ns[4], mean = .$means[4], sd = sqrt(.$vars[4]))
                                    ),
                                    groups = rep(c(1, 2, 3, 4), times = .$ns),
                                    true_means = .$means,
                                    contrast_names = c(
                                      'SE (1 vs 3)', 
                                      'SE (2 vs 4)'
                                    ),
                                    joint_contrasts = c(
                                      'SE (1 vs 3)', 
                                      'SE (2 vs 4)'
                                    ),
                                    c(-1, 0, 1, 0), # simple effect 1 vs 3
                                    c(0, -1, 0, 1) # simple effect 2 vs 4
                                  )) %>% 
                        t() %>% 
                        as.tibble() %>% 
                        unnest()
    )) %>% 
  select(-data) %>% # remove the data so we can unnest
  unnest() %>% 
  
  group_by(sample_ratio, min_sample, effect, var_ratio) %>% # group by conditions and summarize the results
  group_by(sample_ratio, min_sample, effect, var_ratio, contrast_names) %>% 
  summarize(
    joint_reject_student = mean(joint_rejects_student > 0),
    joint_reject_welch = mean(joint_rejects_welch > 0),
    reject_student = mean(p_student < .05),
    reject_welch = mean(p_welch < .05),
    coverage_student = mean(coverage_student),
    coverage_welch = mean(coverage_welch)
  )

# save simulations
save(no_effects, file = '~/R-projects/Welch rule-ANOVA/no_effects.R')

sqrt(1/60 * sum(2, 10, 10, 10))
sqrt(1/60 * sum(2, 4, 4, 4))
sqrt(1/60 * sum(2, 2, 2, 2))
sqrt(1/60 * sum(4, 2, 2, 2))
sqrt(1/60 * sum(10, 2, 2, 2))

