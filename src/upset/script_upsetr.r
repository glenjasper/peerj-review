library(UpSetR)
library(dplyr)

plot_upset = function(category, groups, columns, width, height) {
  file = paste0('../data/data_', tolower(category), '.txt')
  plus = tolower(stringr::str_replace(stringr::str_replace(toString(groups), ' ', ''), ',', '_'))
  
  txt_file = paste0('upset_', tolower(category), '_', plus, '.txt')
  svg_file = paste0('upset_', tolower(category), '_', plus, '.svg')
  
  data = read.delim(file, sep = '\t', check.names = FALSE)
  head(data)
  
  data = data %>%
    rename(Taxon = all_of(category)) %>%
    replace(is.na(.), 0)
  head(data)
  
  # Remove: Unclassified <> & <> (no <> in NCBI)
  data = data %>%
    mutate(Taxon = case_when(grepl('^[Uu]nclassified', Taxon, perl = TRUE) ~ 'Unclassified',
                             grepl('^Others', Taxon, perl = TRUE) ~ 'Unclassified',
                             grepl('\\(no', Taxon, perl = TRUE) ~ 'Unclassified',
                             grepl('(en blanco)', Taxon, perl = TRUE) ~ 'Unclassified',
                             TRUE ~ Taxon)) %>%
    filter(Taxon != 'Unclassified')
  head(data)
  
  data = data %>%
    filter(Domain %in% groups) %>%
    mutate(Taxon = paste0(tolower(stringr::str_sub(Domain, 1, 1)), ':', Taxon)) %>%
    select(Taxon, any_of(columns))
  head(data)
  
  data = data %>%
    group_by(Taxon) %>%
    summarise(across(everything(), ~ sum(.)), .groups = 'drop')
  head(data)
  
  # Remove rows with zeros
  data = data %>%
    mutate(flag = rowSums(across(where(is.numeric)))) %>%
    filter(flag > 0) %>%
    select(-flag)
  head(data)
  
  # Remove columns with zeros
  data = data[, colSums(data != 0) > 0]
  head(data)
  
  data_bin = data %>%
    mutate(across(where(is.numeric), ~ replace(.x, .x > 0, 1)))
  head(data_bin)
  
  write.table(data_bin, txt_file, sep = '\t', quote = FALSE, row.names = FALSE)
  
  list_taxon = list()
  for (column in columns) {
    sample_values = data_bin %>%
      select(Taxon, all_of(column)) %>%
      filter(if_all(where(is.numeric), ~ .x > 0)) %>%
      pull(Taxon)
    
    list_taxon[[column]] = sample_values
  }
  
  svg(filename = svg_file, width = width, height = height)
  
  if (tolower(category) == 'phylum') { label = 'Number of phyla' }
  if (tolower(category) == 'class') { label = 'Number of classes' }
  if (tolower(category) == 'order') { label = 'Number of orders' }
  if (tolower(category) == 'family') { label = 'Number of families' }
  if (tolower(category) == 'genus') { label = 'Number of genera' }
  if (tolower(category) == 'species') { label = 'Number of species' }
  
  upset(data = fromList(list_taxon),
        nsets = length(list_taxon),
        sets.bar.color = '#0080c0', # #56B4E9
        # empty.intersections = 'on',
        order.by = 'freq', # degree
        # group.by = 'sets',
        nintersects = NA,
        matrix.color = '#086e6e',
        point.size = 5,
        sets.x.label = label,
        line.size = 0.75,
        text.scale = c(1.25, # intersection size title
                       1.25, # intersection size tick labels
                       1.25, # set size title
                       1.25, # set size tick labels
                       1.25, # set names
                       1.5)  # numbers above bars
        # sets.bar.color = c('maroon', 'blue', 'orange'),
  )
  
  # dev.off()
}

columns = c('Crude Oil','PHE','OIL','PHE+OIL','PYR','PHE+PYR','FLU','PHE+FLU','HEX','PHE+HEX','PD5','PD6','PF7','NOR_Island_W','PRI_Island_W','SAN_Island_W','FAR_Island_C','NOR_Island_C')
plot_upset(category = 'Phylum', groups = c('Archaea', 'Bacteria'), columns = columns, width = 7, height = 12); dev.off()

columns = c('PD5','PD6','PF7','NOR_Island_W','PRI_Island_W','SAN_Island_W','FAR_Island_C','NOR_Island_C')
plot_upset(category = 'Genus', groups = c('Archaea', 'Bacteria'), columns = columns, width = 11, height = 6); dev.off()

columns = c('Crude Oil','PD5','PD6','PF7','NOR_Island_W','PRI_Island_W','SAN_Island_W','FAR_Island_C','NOR_Island_C')
plot_upset(category = 'Phylum', groups = c('Fungi'), columns = columns, width = 3.2, height = 6.5); dev.off()

columns = c('PD5','PD6','PF7','NOR_Island_W','PRI_Island_W','SAN_Island_W','FAR_Island_C','NOR_Island_C')
plot_upset(category = 'Genus', groups = c('Fungi'), columns = columns, width = 10, height = 6); dev.off()
