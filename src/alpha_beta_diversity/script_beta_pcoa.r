library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(data.table)

get_pcoa = function(data, category, columns = c(), project = NA, prefix = NA, unclassified = F, c1 = 1, c2 = 2, criteria = 'sample', k_minus = 1) {
  label_unclassified = if_else(unclassified, '_with_unclassified', '')
  label_criteria = paste0('_', criteria)
  file_png = paste0('pcoa_c', c1, c2, '_', project, '_', if_else(is.na(prefix), '', paste0(prefix, '_')), category, label_criteria, label_unclassified, '.png')
  file_svg = paste0('pcoa_c', c1, c2, '_', project, '_', if_else(is.na(prefix), '', paste0(prefix, '_')), category, label_criteria, label_unclassified, '.svg')
  file_txt = paste0('table_abundance_pcoa_', project, '_', if_else(is.na(prefix), '', paste0(prefix, '_')), category, label_criteria, label_unclassified, '.txt')
  
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
    summarise(across(everything(), ~ sum(.))) %>%
    as.data.frame
  head(data)
  
  # Remove rows with zeros
  data = data %>%
    mutate(flag = rowSums(across(where(is.numeric)))) %>%
    filter(flag > 0) %>%
    select(-flag)
  head(data)
  
  write.table(data, file_txt, sep = '\t', quote = FALSE, row.names = FALSE)
  
  # Transpose
  data = data %>%
    column_to_rownames('taxa')
  head(data)
  
  data_t = transpose(data)
  rownames(data_t) = colnames(data)
  colnames(data_t) = rownames(data)
  head(data_t)
  
  data_t = decostand(data_t, method = 'hellinger')
  
  data_log = log(data_t + 1)
  data_bray = vegdist(data_log, method = 'bray')
  data_pcoa = cmdscale(data_bray, k = nrow(data_t) - k_minus, eig = T) 
  autovalues = data_pcoa$eig
  
  explication = (autovalues/sum(autovalues)) * 100
  explication_pcoa1 = paste0('PCoA', c1, ' (', round(explication[c1], 2), '%)')
  explication_pcoa2 = paste0('PCoA', c2, ' (', round(explication[c2], 2), '%)')
  
  ###############################
  # Plot with ggplot2
  ###############################
  
  data_dims = scores(data_pcoa)[, c(c1, c2)]
  data_dims = as.data.frame(data_dims)
  
  if(project == 'trindade_review') {
    data_dims = data_dims %>%
      rownames_to_column('replica') %>%
      mutate(sample = case_when(replica %in% c('Crude Oil') ~ 'Soil_MC',
                                replica %in% c('PD5', 'PD6', 'PF7') ~ 'Soil_Env',
                                replica %in% c('FLU', 'HEX', 'OIL', 'PHE', 'PHE+FLU', 'PHE+HEX', 'PHE+OIL', 'PHE+PYR', 'PYR') ~ 'Water_MC',
                                replica %in% c('FAR_Island_C', 'NOR_Island_C') ~ 'Coral_Env',
                                replica %in% c('NOR_Island_W', 'PRI_Island_W', 'SAN_Island_W') ~ 'Water_Env'))
    
    
    # Factor
    data_dims$sample = factor(data_dims$sample, levels = c('Soil_MC', 'Soil_Env', 'Water_MC', 'Water_Env', 'Coral_Env'))
  }

  # Colors
  if(project == 'trindade_review') {
    colors = c('#7ABE29',
               '#E69F00',
               '#56B4E9',
               '#990099',
               '#FF68C9',
               '#005c3f',
               '#00ccca',
               '#fd8b73')
  }
  
  if(project == 'trindade') {
    n_groups = 1 # 3
    colors = c(rep('#DF0070', n_groups),
               rep('#00A23F', n_groups),
               rep('#ABB02B', n_groups),
               rep('#009aca', n_groups)) # cb09cb
  }

  rownames(data_dims) = NULL
  data_dims = column_to_rownames(data_dims, 'replica')
  
  o_ggplot = NA
  if (c1 == 1 & c2 == 2) {
    ggplot(data_dims, aes(x = Dim1, y = Dim2, color = sample)) -> o_ggplot
  }
  if (c1 == 2 & c2 == 3) {
    ggplot(data_dims, aes(x = Dim2, y = Dim3, color = sample)) -> o_ggplot
  }
  if (c1 == 3 & c2 == 4) {
    ggplot(data_dims, aes(x = Dim3, y = Dim4, color = sample)) -> o_ggplot
  }
  
  o_ggplot +
    geom_hline(yintercept = 0, linetype = 'dashed', color = '#9F9F9F') +
    geom_vline(xintercept = 0, linetype = 'dashed', color = '#9F9F9F') +
    geom_point() +
    geom_text_repel(label = rownames(data_dims),
                    color = 'black',
                    size = 2,
                    segment.size = 0.15,
                    max.overlaps = Inf) +
    scale_color_manual(values = colors) +
    labs(x = explication_pcoa1, y = explication_pcoa2) +
    theme(legend.title = element_blank(),
          legend.key = element_blank(),
          legend.key.width = unit(0.1, 'cm'),
          legend.key.height = unit(0.35, 'cm'),
          legend.text = element_text(size = 7),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm'), # top, right, bottom, left
          axis.text.x = element_text(size = 7, color = '#000000'),
          axis.text.y = element_text(size = 7, color = '#000000'),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.ticks = element_line(linewidth = 0.35,  color = '#000000'),
          axis.ticks.length = unit(2.5, 'pt'),
          panel.background = element_blank(),
          panel.border = element_rect(linewidth = 1, color = '#000000', fill = NA)) -> data_plot
  
  ggsave(filename = file_png, plot = data_plot, width = 6, height = 4, dpi = 600, units = 'in', device = 'png')
  ggsave(filename = file_svg, plot = data_plot, width = 6, height = 4, dpi = 600, units = 'in', device = 'svg')
}

file = 'data/data_review_trindade_order.txt'
data = read.delim(file, sep = '\t', check.names = F)
head(data)

get_pcoa(project = 'trindade_review', data = data, category = 'order', k_minus = 3)
