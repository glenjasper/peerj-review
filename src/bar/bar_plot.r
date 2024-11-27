library(ggplot2)
library(reshape2)
library(dplyr)

plot_taxonomy = function(method, tax, groups, columns, limit = 0, width = 3.5) {
  prefix = ''
  for (domain in groups) {
    i = tolower(substring(domain, 1, 1))
    prefix = paste0(prefix, i)
  }
  
  n = case_when(tax == 'Domain' ~ '1',
                tax == 'Phylum' ~ '2',
                tax == 'Class' ~ '3',
                tax == 'Order' ~ '4',
                tax == 'Family' ~ '5',
                tax == 'Genus' ~ '6',
                tax == 'Species' ~ '7',
                TRUE ~ '')
  
  file = paste0('../data/data_', tolower(tax), '.txt')
  
  png_file = paste0('barplot-', tolower(method), '-', prefix, '-', n, '-', tolower(tax), if (limit > 0) {paste0('-others-', limit)}, '.png')
  svg_file = paste0('barplot-', tolower(method), '-', prefix, '-', n, '-', tolower(tax), if (limit > 0) {paste0('-others-', limit)}, '.svg')
  txt_file_a = paste0('table_abundances_abs_', tolower(method), '-', prefix, '-', n, '-', tolower(tax), '.txt')
  txt_file_r = paste0('table_abundances_rel_', tolower(method), '-', prefix, '-', n, '-', tolower(tax), '.txt')
  
  data = read.delim(file, sep = '\t', check.names = FALSE)
  head(data)
  
  data = data %>%
    rename(Taxon = all_of(tax)) %>%
    replace(is.na(.), 0)
  head(data)
  
  # Replace 'uncultured <tax>' and '' to 'Unclassified'
  data = data %>%
    mutate(Taxon = if_else(Taxon == '', 'Unclassified', Taxon),
           Taxon = if_else(Taxon == '(en blanco)', 'Unclassified', Taxon),
           Taxon = if_else(grepl('unclassified', Taxon, perl = TRUE), 'Unclassified', Taxon)
    )
  head(data)
  # str(data)
  
  if (tax == 'Domain') {
    # Grouping + sum
    if (method %in% c('Soil', 'Water', 'All')) {
      data = data %>%
        select(-Method) %>%
        group_by(Taxon) %>%
        summarise(across(everything(), ~ sum(.)), .groups = 'drop')
      data = as.data.frame(data)
      head(data)
      
      data = data %>%
        filter(Taxon %in% groups) %>%
        select(Taxon, any_of(columns))
      head(data)
    } else {
      data = data %>%
        group_by(Method, Taxon) %>%
        summarise(across(everything(), ~ sum(.)), .groups = 'drop')
      data = as.data.frame(data)
      head(data)
      
      data = data %>%
        filter(Method == method) %>%
        filter(Taxon %in% groups) %>%
        select(Taxon, any_of(columns))
      head(data)
    }
  } else {
    # Filter groups
    if (method %in% c('Soil', 'Water', 'All')) {
      data = data %>%
        filter(Domain %in% groups)
      head(data)
    } else {
      data = data %>%
        filter(Method == method) %>%
        filter(Domain %in% groups)
      head(data)
    }
    
    data = data %>%
      mutate(Taxon = if_else(Taxon != 'Unclassified', paste0(substr(tolower(Domain), 1, 1), ':', Taxon), Taxon))
    head(data)
    
    data = data %>%
      select(Taxon, any_of(columns))
    head(data)
    
    # Grouping + sum
    data = data %>%
      group_by(Taxon) %>% # Method, Domain, Taxon
      summarise(across(everything(), ~ sum(.)), .groups = 'drop')
    data = as.data.frame(data)
    head(data)
  }
  
  # Remove rows with zeros
  data = data %>%
    mutate(flag = rowSums(across(where(is.numeric)))) %>%
    filter(flag > 0) %>%
    select(-flag)
  head(data)
  
  # Remove columns with zeros
  data = data[, colSums(data != 0) > 0]
  head(data)
  
  # Rename Others_x
  data = data %>%
    mutate(Taxon = if_else(Taxon == 'a:Others_a', 'a:Others (Archaea)', Taxon),
           Taxon = if_else(Taxon == 'b:Others_b', 'b:Others (Bacteria)', Taxon),
           Taxon = if_else(Taxon == 'f:Others_f', 'f:Others (Fungi)', Taxon))
  head(data)
  
  # Save file
  write.table(data, txt_file_a, sep = '\t', quote = FALSE, row.names = FALSE)
  
  data_long = reshape2::melt(data, id = c('Taxon'))
  head(data_long)
  
  # To relative percentage
  data_long = data_long %>%
    group_by(variable) %>%
    mutate(value = 100 * value / sum(value))
  data_long = as.data.frame(data_long)
  head(data_long)
  
  # For NaN
  data_long = data_long %>%
    replace(is.na(.), 0)
  head(data_long)
  
  # Save file
  data_casted = reshape2::dcast(data_long, Taxon ~ variable)
  write.table(data_casted, txt_file_r, sep = '\t', quote = FALSE, row.names = FALSE)
  
  if (limit > 0){
    label_others = paste0('Others (<', limit, '%)')
    
    # Set Others group
    data_long = data_long %>%
      mutate(Taxon = if_else((value < limit) & (Taxon != 'Unclassified'), label_others, Taxon)) #  & (value != 0)
    head(data_long)
    
    # Grouping Others
    data_long = data_long %>%
      group_by(Taxon, variable) %>%
      summarize(value = sum(value, na.rm = TRUE))
    head(data_long)
    
    if (tax != 'Domain'){
      # For Archaea
      names_archaea = data_long %>%
        filter(startsWith(Taxon, 'a:')) %>%
        mutate(Taxon = substring(Taxon, 3)) %>%
        select(Taxon) %>%
        distinct(Taxon) %>%
        arrange(Taxon)
      head(names_archaea)
      names_archaea = names_archaea[['Taxon']]
      
      # For Bacteria
      names_bacteria = data_long %>%
        filter(startsWith(Taxon, 'b:')) %>%
        mutate(Taxon = substring(Taxon, 3)) %>%
        select(Taxon) %>%
        distinct(Taxon) %>%
        arrange(Taxon)
      head(names_bacteria)
      names_bacteria = names_bacteria[['Taxon']]
      
      # For Fungi
      names_fungi = data_long %>%
        filter(startsWith(Taxon, 'f:')) %>%
        mutate(Taxon = substring(Taxon, 3)) %>%
        select(Taxon) %>%
        distinct(Taxon) %>%
        arrange(Taxon)
      head(names_fungi)
      names_fungi = names_fungi[['Taxon']]
      
      labels_order = c(names_archaea, names_bacteria, names_fungi, label_others)
    }
  } else {
    if (tax != 'Domain'){
      # For Archaea
      names_archaea = data %>%
        filter(startsWith(Taxon, 'a:')) %>%
        mutate(Taxon = substring(Taxon, 3)) %>%
        select(Taxon) %>%
        arrange(Taxon)
      head(names_archaea)
      names_archaea = names_archaea[['Taxon']]
      
      # For Bacteria
      names_bacteria = data %>%
        filter(startsWith(Taxon, 'b:')) %>%
        mutate(Taxon = substring(Taxon, 3)) %>%
        select(Taxon) %>%
        arrange(Taxon)
      head(names_bacteria)
      names_bacteria = names_bacteria[['Taxon']]
      
      # For Fungi
      names_fungi = data %>%
        filter(startsWith(Taxon, 'f:')) %>%
        mutate(Taxon = substring(Taxon, 3)) %>%
        select(Taxon) %>%
        arrange(Taxon)
      head(names_fungi)
      names_fungi = names_fungi[['Taxon']]
      
      labels_order = c(names_archaea, names_bacteria, names_fungi)
    }
  }
  
  data_long = data_long %>%
    mutate(Taxon = if_else(startsWith(Taxon, 'a:') | startsWith(Taxon, 'b:') | startsWith(Taxon, 'f:'), substring(Taxon, 3), Taxon))
  head(data_long)
  
  last_index = length(unique(data_long$Taxon))
  
  flag_unclassified = FALSE
  flag_others_bacteria = FALSE
  flag_others_fungi = FALSE
  flag_others = FALSE
  
  # Redefine Factor
  data_long$variable = factor(data_long$variable, levels = columns)
  if (tax != 'Domain'){
    # Flags
    if ('Others (Bacteria)' %in% unique(data_long$Taxon)) {
      flag_others_bacteria = TRUE
      
      labels_order = labels_order[labels_order != 'Others (Bacteria)']
      labels_order = c(labels_order, 'Others (Bacteria)')
    }
    if ('Others (Fungi)' %in% unique(data_long$Taxon)) {
      flag_others_fungi = TRUE
      
      labels_order = labels_order[labels_order != 'Others (Fungi)']
      labels_order = c(labels_order, 'Others (Fungi)')
    }
    if (limit > 0) {
      if (label_others %in% unique(data_long$Taxon)) {
        flag_others = TRUE
        
        labels_order = labels_order[labels_order != label_others]
        labels_order = c(labels_order, label_others)
      }
    }
    if ('Unclassified' %in% unique(data_long$Taxon)) {
      flag_unclassified = TRUE
      labels_order = c(labels_order, 'Unclassified')
    }
    
    data_long$Taxon = factor(data_long$Taxon, levels = labels_order)
  }
  
  colors = c('#005c3f',
             '#00668a',
             '#00ccca',
             '#fd8b73',
             '#b27d93',
             '#aaccaa',
             '#ab992b',
          #   '#571c34',
             '#cb09cb',
             '#e6dd50',
             '#6aa99c',
             '#443890',
             '#009aca',
             '#af062b',
             '#07a16c',
          #   '#f8184a',
             '#f0a911',
             '#DF0070',
             '#006e4f',
             '#ff9999',
             '#009798',
             '#aa00aa',
             '#ff2726',
             '#ffcc00',
             '#00d400',
             '#DF7770',
             '#73f69d',
             '#99bbaa',
             '#0933aa',
             '#506f00',
             '#980100',
             '#b5b517',
             '#aa66ff',
             '#ffff20',
             '#e65757',
             '#00ddee',
             '#7d4f20',
             '#f4a460',
             '#7548d0',
             '#3cb371',
             '#770077',
             '#274965',
             '#228b22',
             '#09332f',
             '#d4c1ff',
             '#6F6FFF',
             '#A680FF',
             '#BF7831',
             '#DDDDDD',
             '#AAAAAA',
             '#484848')
  
  # Set colors: Unclassified and Others
  color_others_bf = '#DDDDDD'
  color_others = '#AAAAAA'
  color_unclassified = '#484848'
  
  if ((flag_others_bacteria | flag_others_fungi)) {
    if (flag_others) {
      if (flag_unclassified) {
        colors[last_index - 2] = color_others_bf
        colors[last_index - 1] = color_others
        colors[last_index] = color_unclassified
      } else {
        colors[last_index - 1] = color_others_bf
        colors[last_index] = color_others
      }
    } else {
      if (flag_unclassified) {
        colors[last_index - 1] = color_others_bf
        colors[last_index] = color_unclassified
      } else {
        colors[last_index] = color_others_bf
      }
    }
  } else {
    if (flag_others) {
      if (flag_unclassified) {
        colors[last_index - 1] = color_others
        colors[last_index] = color_unclassified
      } else {
        colors[last_index] = color_others
      }
    } else {
      if (flag_unclassified) {
        colors[last_index] = color_unclassified
      } else {
        # NULL
      }
    }
  }
  
  #################
  # Stack to Fill
  #################
  ggplot(data_long, aes(x = variable, y = value, fill = Taxon)) +
    # geom_bar(position = position_fill(reverse = TRUE), stat = "identity", colour = "black", linewidth = 0.3) +
    geom_bar(position = "fill", stat = "identity", colour = "black", linewidth = 0) +
    scale_y_continuous(name = "Relative abundance (%)",
                       expand = c(0, 0),
                       breaks = seq(0, 1, .1),
                       # position = "right",
                       labels = scales::percent_format()) +
    scale_fill_manual(values = colors) +
    labs(x = NULL) +
    theme(legend.title = element_blank(),
          legend.key = element_blank(),
          # legend.key.size = unit(0.2, "cm"),
          legend.key.width = unit(0.125, "cm"),
          legend.key.height = unit(0.2, "cm"),
          legend.text = element_text(size = 4), # family = "mono"
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm'), # top, right, bottom, left
          axis.text.x = element_text(size = 3.5, colour = "#000000", angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 5, colour = "#000000"),
          # axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.ticks = element_line(linewidth = 0.35,  colour = "#000000"),
          axis.ticks.length = unit(2.5, 'pt'),
          axis.line = element_line(linewidth = 0.35,  colour = "#000000"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank()) +
    guides(fill = guide_legend(override.aes = list(color = NA))) -> data_plot
  
  ggsave(filename = png_file, plot = data_plot, width = width, height = 2, dpi = 800, units = 'in', device = 'png')
  ggsave(filename = svg_file, plot = data_plot, width = width, height = 2.25, dpi = 800, units = 'in', device = 'svg')
}

method = 'All'
groups = c('Archaea', 'Bacteria')

columns = c('Crude Oil','PHE','OIL','PHE+OIL','PYR','PHE+PYR','FLU','PHE+FLU','HEX','PHE+HEX','PD5','PD6','PF7','NOR_Island_W','PRI_Island_W','SAN_Island_W','FAR_Island_C','NOR_Island_C')
# plot_taxonomy(method = method, tax = 'Domain', groups = groups, columns = columns, limit = 0, width = 2.65)
plot_taxonomy(method = method, tax = 'Phylum', groups = groups, columns = columns, limit = 1, width = 3.5)
# plot_taxonomy(method = method, tax = 'Class', groups = groups, columns = columns, limit = 2, width = 3.65)
# plot_taxonomy(method = method, tax = 'Order', groups = groups, columns = columns, limit = 3, width = 3.825)
# plot_taxonomy(method = method, tax = 'Family', groups = groups, columns = columns, limit = 3, width = 3.65)

columns = c('PHE','OIL','PHE+OIL','PYR','PHE+PYR','FLU','PHE+FLU','HEX','PHE+HEX','PD5','PD6','PF7','NOR_Island_W','PRI_Island_W','SAN_Island_W','FAR_Island_C','NOR_Island_C')
plot_taxonomy(method = method, tax = 'Genus', groups = groups, columns = columns, limit = 1.5, width = 2.65)
# plot_taxonomy(method = method, tax = 'Species', groups = groups, columns = columns, limit = 1, width = 3.5)

groups = c('Fungi')

columns = c('Crude Oil','PHE','OIL','PHE+OIL','PYR','PHE+PYR','FLU','PHE+FLU','HEX','PHE+HEX','PD5','PD6','PF7','NOR_Island_W','PRI_Island_W','SAN_Island_W','FAR_Island_C','NOR_Island_C')
# plot_taxonomy(method = method, tax = 'Domain', groups = groups, columns = columns, limit = 0, width = 1.75)
plot_taxonomy(method = method, tax = 'Phylum', groups = groups, columns = columns, limit = 1, width = 2)
# plot_taxonomy(method = method, tax = 'Class', groups = groups, columns = columns, limit = 1, width = 2.15)
# plot_taxonomy(method = method, tax = 'Order', groups = groups, columns = columns, limit = 3, width = 2.85)
# plot_taxonomy(method = method, tax = 'Family', groups = groups, columns = columns, limit = 3, width = 2.775)
plot_taxonomy(method = method, tax = 'Genus', groups = groups, columns = columns, limit = 3, width = 2.6)
# plot_taxonomy(method = method, tax = 'Species', groups = groups, columns = columns, limit = 3, width = 3.2)
