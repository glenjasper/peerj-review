library(vegan)
library(ggplot2)
library(dplyr)
library(data.table)
library(tibble)
library(ggrepel)

get_box_diversity = function(data, category, columns = c(), test = 'shannon', plus = 0, project = NA, prefix = NA, width = 5, round = 2, unclassified = F) {
  label_unclassified = if_else(unclassified, '_with_unclassified', '')
  file_png = paste0('boxplot_', project, '_', if_else(is.na(prefix), '', paste0(prefix, '_')), category, '_', test, label_unclassified, '.png')
  file_svg = paste0('boxplot_', project, '_', if_else(is.na(prefix), '', paste0(prefix, '_')), category, '_', test, label_unclassified, '.svg')
  file_txt = paste0('table_', project, '_', if_else(is.na(prefix), '', paste0(prefix, '_')), category, '_', test, label_unclassified, '.txt')
  title_y = paste0('Alpha-diversity Index: ', stringr::str_to_title(test))
  
  if (category %in% colnames(data)) {
    data = data %>%
      rename(taxa = all_of(category))
    head(data)
  }
  
  if (length(columns) > 0) {
    data = data %>%
      select(taxa, all_of(columns))
    head(data)
  }
  
  if (!unclassified) {
    data = data %>%
      mutate(taxa = case_when(grepl('^Unclassified', taxa, perl = TRUE) ~ 'Unclassified',
                              grepl('Others', taxa, perl = TRUE) ~ 'Unclassified',
                              grepl('\\(no', taxa, perl = TRUE) ~ 'Unclassified',
                              taxa == '' ~ 'Unclassified',
                              TRUE ~ taxa)) %>%
      filter(taxa != 'Unclassified')
    head(data)
  }
  
  # Grouping + sum
  data = data %>%
    group_by(taxa) %>%
    summarise(across(everything(), ~ sum(.)))
  head(data)
  
  # Remove rows with zeros
  data = data %>%
    mutate(flag = rowSums(across(where(is.numeric)))) %>%
    filter(flag > 0) %>%
    select(-flag)
  head(data)
  
  # Transpose
  data = data %>%
    column_to_rownames('taxa')
  head(data)
  
  data_t = transpose(data)
  rownames(data_t) = colnames(data)
  colnames(data_t) = rownames(data)
  head(data_t)
  
  data_input = data_t
  if (test != 'chao1') {
    data_input = decostand(data_t, method = 'total')
  }
  
  data_stats = data.frame()
  switch(test,
         'shannon' = {
           # Calculate Shannon index (H)
           shannon = diversity(data_input, index = 'shannon')
           df_shannon = as.data.frame(shannon)
           colnames(df_shannon) = 'value'
           data_stats = df_shannon
         },
         'simpson' = {
           # Calculate Simpson index [1-Dominance]
           simpson = diversity(data_input, index = 'simpson')
           df_simpson = as.data.frame(simpson)
           colnames(df_simpson) = 'value'
           data_stats = df_simpson
         },
         'evenness' = {
           # Calculate Shannon index (H)
           shannon = diversity(data_input, index = 'shannon')
           
           # Calculate Evenness index (J) [Equitability | Pielou] 
           evenness = shannon/log(specnumber(data_input))
           df_evenness = as.data.frame(evenness)
           colnames(df_evenness) = 'value'
           data_stats = df_evenness
         },
         'chao1' = {
           # Calculate Richness and Chao1
           richness = estimateR(data_input)
           df_richness = as.data.frame(t(richness))
           colnames(df_richness)[2] = 'value'
           data_stats = df_richness
         },
         {
           stop(paste0("Test ", test, " does't exist."))
         }
  )
  
  data_stats = data_stats %>%
    rownames_to_column('sample')
  head(data_stats)
  
  colors = c()
  if(project == 'trindade_review') {
    data_stats = data_stats %>%
      mutate(group = case_when(sample %in% c('Crude Oil') ~ 'Soil_MC',
                               sample %in% c('PD5', 'PD6', 'PF7') ~ 'Soil_Env',
                               sample %in% c('FLU', 'HEX', 'OIL', 'PHE', 'PHE+FLU', 'PHE+HEX', 'PHE+OIL', 'PHE+PYR', 'PYR') ~ 'Water_MC',
                               sample %in% c('FAR_Island_C', 'NOR_Island_C') ~ 'Coral_Env',
                               sample %in% c('NOR_Island_W', 'PRI_Island_W', 'SAN_Island_W') ~ 'Water_Env'))
    # Factor
    data_stats$group = factor(data_stats$group, levels = c('Soil_MC', 'Soil_Env', 'Water_MC', 'Water_Env', 'Coral_Env'))
    
    colors = c('#7ABE29',
               '#E69F00',
               '#56B4E9',
               '#990099',
               '#FF68C9')
  }

  # Save file
  data_stats_txt = data_stats %>%
    # rownames_to_column('sample') %>%
    select(sample, value, group) %>%
    rename(!!quo_name(test) := value)
  head(data_stats_txt)
  
  write.table(data_stats_txt, file_txt, sep = '\t', quote = FALSE, row.names = FALSE)
  
  # Means
  means = aggregate(value ~ group, data_stats, mean)
  means$value = round(means$value, round)
  
  value_plus = mean(means$value) * plus
  
  ggplot(data_stats, aes(x = group, y = value, color = group)) + # fill = group
    geom_boxplot(size = .5, width = 0.75, outlier.shape = NA) + # don't show outliers
    geom_jitter(width = 0.25, size = 0.75) + # 0: One column # show outliers
    geom_text_repel(aes(label = sample),
                    color = 'black',
                    size = 2,
                    segment.size = 0.15,
                    max.overlaps = Inf) +
    stat_summary(fun = mean, colour = 'black', geom = 'point', 
                 shape = 18, size = 2, show.legend = FALSE) +
    geom_text(data = means, aes(label = value, y = value + value_plus),
              position = position_dodge(width = 2), size = 2, color = 'black') +
    # scale_color_manual(values = colors) +
    labs(x = NULL, y = title_y) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.key = element_blank(),
          axis.text.x = element_text(size = 6, colour = '#000000', angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 8, colour = '#000000'),
          axis.title.y = element_text(size = 10, colour = '#000000'),
          axis.line = element_line(linewidth = 0.5, colour = '#000000'),
          panel.background = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2, colour = '#EEF0F3', linetype = 'dashed'),
          panel.grid.minor = element_line(linewidth = 0.2, colour = '#EEF0F3', linetype = 'dashed'),
          panel.border = element_blank()) -> sample_plot
  
  if (length(colors) > 0) {
    sample_plot + scale_color_manual(values = colors) -> sample_plot
  }
  
  ggsave(filename = file_png, plot = sample_plot, width = width, height = 4, dpi = 800, units = 'in', device = 'png')
  ggsave(filename = file_svg, plot = sample_plot, width = width, height = 4, dpi = 800, units = 'in', device = 'svg')
}

file = 'data/data_review_trindade_order.txt'
data = read.delim(file, sep = '\t', check.names = F)
head(data)

get_box_diversity(project = 'trindade_review', prefix = 'trev', data = data, category = 'order', test = 'shannon', plus = 0.04, round = 3)
# get_box_diversity(project = 'trindade_review', prefix = 'trev', data = data, category = 'order', test = 'simpson', plus = 0.04)
# get_box_diversity(project = 'trindade_review', prefix = 'trev', data = data, category = 'order', test = 'evenness', plus = 0.04)
# get_box_diversity(project = 'trindade_review', prefix = 'trev', data = data, category = 'order', test = 'chao1', plus = 0.04)
