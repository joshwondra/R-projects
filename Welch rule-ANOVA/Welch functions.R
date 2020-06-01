# legend name ----
legend_name <- expression(paste(
  n[A[1]*B[1]], 
  ', ',
  n[A[1]*B[2]],
  ', ',
  n[A[2]*B[1]],
  ', ',
  n[A[2]*B[2]]
  ))
  

# Functions for Rnw - Welch Manuscript - ANOVA

plot_theme <- function(){
  theme_classic() +
    theme(
      text = element_text(family = 'Helvetica'),
      title = element_text(size = 16, face = 'bold'),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14),
      strip.background = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(color = 'gray80', linetype = 'longdash')
    )
}

true_effect_plot <- function(data, which_effect) {
  
  data %>% 
    filter(effect == which_effect) %>% 
    ggplot(aes(y = y, x = factorA, group = factorB, linetype = factorB)) +
    geom_line(show.legend = FALSE) +
    geom_text(
      data = data %>% 
        filter(
          effect == which_effect, 
          !is.na(factorB_label_position)
        ),
      aes(y = factorB_label_position, label = factorB),
      size = 14*(5/14) # 5/14 is a rough conversion to make geom_text() size similar to element_text() size
    ) +
    coord_cartesian(ylim = c(5, 8)) +
    labs(title = which_effect, x = NULL, y = NULL) +
    plot_theme() +
    theme(
      panel.grid.major.y = element_blank()
    )
  
}

# useful to fill some spaces in plot_grid()
blank_plot <- function() {
  ggplot(data = data.frame()) +
    geom_blank() +
    theme(rect = element_blank())
}

blank_plot_cov <- function(cov_data) {
  cov_data %>% 
  filter(
    which_sim == 'No effects',
    contrast_names == 'Main Effect of A',
    test == 'ANOVA'
  ) %>% 
  ggplot(aes(y = coverage_rate, x = var_ratio)) +
  geom_blank() +
  scale_y_continuous (limits = c(0, .15)) +
  labs(x = NULL, y = NULL) +
  theme(
    rect = element_blank(),
    line = element_blank()
  )
}

# need to adjust simple effect labels
se_labels <- function(string) {
  string = str_replace(string, 'Simple Effect of ', 'Simple Effect of\n')
  return(string)
}

# specify what's common to all the plots
format_plot <- function(plot, ylims, ybreaks, ylabs){
  plot +
    geom_line(size = 1) +
    geom_point(size = 9, stroke = 1, color = 'white', shape = 16) + # background color
    geom_point(size = 5, stroke = 1) +
    scale_y_continuous(limits = ylims, breaks = ybreaks, labels = ylabs) +
    scale_shape_manual(
      name = legend_name,
      labels = sample_sizes,
      values = rep(c(0, 1, 2), times = 3)
    ) +
    scale_linetype_manual(
      name = legend_name,
      labels = sample_sizes,
      values = c('solid', 'twodash', 'dotted')
    ) +
    guides(linetype = FALSE, shape = FALSE) +
    facet_grid(cols = vars(contrast_names), rows = vars(test), labeller = labeller(contrast_names = se_labels)) +
    labs(y = NULL, x = NULL) +
    plot_theme()
}

# plots reject rates for a single test
reject_null_plot <- function(data, sim, which_contrast, which_test) {
  
  # set y limits and breaks
  if (
    sim == 'No effects' | 
    sim == 'Crossover interaction' & str_detect(which_contrast, 'Main') |
    sim == 'Main effects & interaction' & which_contrast == 'Simple Effect of B at A2'
  ) {
    ylims = c(0, .25)
    ybreaks = seq(0, 1, .05)
  } else if (
    sim == 'Crossover interaction' & str_detect(which_contrast, 'Int') |
    sim == 'Main effects & interaction'
  ) {
    ylims = c(0, 1)
    ybreaks = seq(0, 1, .1)
  } else {
    ylims = c(0, .5)
    ybreaks = seq(0, 1, .1)
  }
  
  # set y labels
  ylabs = paste0(ybreaks * 100, '%')
  
  temp_plot = data %>% 
    filter(
      which_sim == sim,
      contrast_names == which_contrast,
      test == which_test
    ) %>% 
    ggplot(aes(y = reject_rate, x = var_ratio, group = ns, linetype = ns, shape = ns))
  
  format_plot(temp_plot, ylims, ybreaks, ylabs)
} 

# plots difference in reject rates for a single test
reject_difference_plot <- function(data, sim, which_contrast) {
  
  ylims = c(-.2, .2)
  ybreaks = seq(-1, 1, .05) %>% round(2)
  ylabs = paste0(ybreaks * 100, '%')
  
  temp_plot = data %>% 
    filter(
      which_sim == sim,
      contrast_names == which_contrast
    ) %>% 
    spread(key = test, value = reject_rate) %>% 
    mutate(
      power_diff = `Welch's t test` - ANOVA,
      test = "Welch's - ANOVA"
    ) %>% 
    ggplot(aes(y = power_diff, x = var_ratio, group = ns, linetype = ns, shape = ns))
  
  format_plot(temp_plot, ylims, ybreaks, ylabs)
} 





coverage_plot <- function(data, sim, which_contrast, which_test) {
  
  ylims = c(.75, 1)
  ybreaks = seq(0, 1, .05)
  ylabs = paste0(ybreaks * 100, '%')
  
  temp_plot = data %>% 
    filter(
      which_sim == sim,
      contrast_names == which_contrast,
      test == which_test
    ) %>% 
    ggplot(aes(y = coverage_rate, x = var_ratio, group = ns, linetype = ns, shape = ns))
  
  format_plot(temp_plot, ylims, ybreaks, ylabs)
} 







# function to use in apply() with the table and name of the plotting function as inputs
# generates 1 plot per row of the table where the columns sim, which_contrast, and which_test are the inputs for the plotting function
generate_plots <- function(data, which_plot) {
  
  if(!which_plot %in% c('reject_null_plot', 'reject_difference_plot', 'coverage_plot')) {
    print('ERROR: Select one of the following plots: reject_null_plot, reject_difference_plot')
    return()
  }
  
  # save variables
  sim = data['sim']
  which_contrast = data['which_contrast']
  
  # create the name of the plot object 
  prefix = case_when(
    which_plot == 'reject_null_plot' ~ 'rn',
    which_plot == 'reject_difference_plot' ~ 'rd',
    which_plot == 'coverage_plot' ~ 'cov'
  )
  
  sim_name = case_when(
    sim == 'No effects' ~ 'no_effects',
    sim == 'Crossover interaction' ~ 'crossover',
    sim == 'Main effects & interaction' ~ 'me_int'
  )
  
  contrast_name = case_when(
    which_contrast == 'Main Effect of A' ~ 'meA',
    which_contrast == 'Main Effect of B' ~ 'meB',
    which_contrast == 'Interaction' ~ 'int',
    which_contrast == 'Simple Effect of B at A1' ~ 'A1cont',
    which_contrast == 'Simple Effect of B at A2' ~ 'A2cont',
    which_contrast == 'Simple Effect of A at B1' ~ 'B1cont',
    which_contrast == 'Simple Effect of A at B2' ~ 'B2cont'
  )
  
  plot_name = paste(prefix, sim_name, contrast_name, sep = '_')

  
  
  # prepare for plot
  plot_data = case_when(
    which_plot %in% c('reject_null_plot', 'reject_difference_plot') ~ 'reject_null_data',
    which_plot == 'coverage_plot' ~ 'cov_data'
  )
  
  plot_args = list(
    data = get(plot_data, envir = .GlobalEnv),
    sim = sim,
    which_contrast = which_contrast
  )
  
  # check if which_test exists and make modifications as needed
  if ('which_test' %in% names(data)) {
    which_test = data['which_test']
    test_name = case_when(
      which_test == 'ANOVA' ~ 'anova',
      which_test == "Welch's t test" ~ 'welch'
    )
    plot_name = paste(plot_name, test_name, sep = '_')
    plot_args[['which_test']] = which_test
  }  

  
  # create the plot object and assign it to an object in the global environment
  
  final_plot = do.call(
    which_plot,
    args = plot_args
  )

  assign(x = plot_name, value = final_plot, envir = .GlobalEnv)
}


# Customization for plot_grid
first_row <- function(plot) {
  plot +
  theme(
    strip.text.y = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = 'gray80', linetype = 'longdash')
  )
}
low_row <- function(plot) {
  plot +
  theme(
    strip.text = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = 'gray80', linetype = 'longdash')
  )
}
low_row_end <- function(plot) {
  plot +
  theme(
    strip.text.x = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = 'gray80', linetype = 'longdash')
  )
}

add_x_lab <- function(plot, x_lab) {
  add_sub(
    plot,
    x_lab, 
    vpadding=grid::unit(1,"lines"),
    y=.1, 
    x=0.5, 
    vjust=0
  ) %>% 
    ggdraw()
}

