library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
library(tibble)

get_tables = function(file, category, type, columns, transform_type) {
  txt_file_info = paste0('data_is_category/table_', tolower(category), '_info_', transform_type, '_', type, '.txt')
  txt_file_transpose = paste0('data_is_category/table_', tolower(category), '_transpose_', transform_type, '_', type, '.txt')
  dir.create('data_is_category', showWarnings = FALSE)
  
  data = read.delim(file, sep = '\t', check.names = FALSE)
  head(data)
  
  data = data %>%
    replace(is.na(.), 0)
  head(data)
  
  data = data %>%
    rename(Taxon = any_of(category)) %>%
    select(Domain, Taxon, any_of(columns)) %>%
    filter(Taxon != '')
  head(data)
  
  # Remove rows with zeros
  data = data %>%
    mutate(flag = rowSums(across(where(is.numeric)))) %>%
    filter(flag > 0) %>%
    select(-flag)
  head(data)
  
  data = data %>%
    group_by(Domain, Taxon) %>%
    summarise(across(everything(), ~ sum(.)), .groups = 'drop')
  data = as.data.frame(data)
  
  data_long = reshape2::melt(data, id.vars = c('Domain', 'Taxon'), value.name = 'value', variable.name = 'variable')
  head(data_long)
  
  transform_value = 0
  if (transform_type == 'relative') { transform_value = 1 }
  if (transform_type == 'hellinger') { transform_value = .5 }
  
  # To relative percentage & Hellinger
  if (type == 'by_domain') {
    data_long = data_long %>%
      group_by(Domain, variable) %>%
      mutate(value = 100 * value / sum(value)) %>%
      mutate(value = value^transform_value)
    data_long = as.data.frame(data_long)
    head(data_long)
  } else {
    data_long = data_long %>%
      group_by(variable) %>%
      mutate(value = 100 * value / sum(value)) %>%
      mutate(value = value^transform_value)
    data_long = as.data.frame(data_long)
    head(data_long)
  }
  
  # For NaN
  data_long = data_long %>%
    replace(is.na(.), 0)
  head(data_long)
  
  # Save file
  data_casted = reshape2::dcast(data_long, Domain + Taxon ~ variable)
  head(data_casted)
  write.table(data_casted, txt_file_info, sep = '\t', quote = FALSE, row.names = FALSE)
  
  data = data_casted %>%
    mutate(Taxon = paste0(substr(tolower(Domain), 1, 1), ':', Taxon)) %>%
    select(-Domain)
  head(data)
  
  # Column to rowname
  data = data %>%
    column_to_rownames('Taxon')
  head(data)
  
  data_t = transpose(data)
  rownames(data_t) = colnames(data)
  colnames(data_t) = rownames(data)
  
  col_tax = colnames(data_t)
  
  data_t = data_t %>%
    rownames_to_column('Sample') %>%
    mutate(Sample_Type = case_when(Sample %in% c('Crude Oil') ~ 'Soil_MC',
                                   Sample %in% c('PD5', 'PD6', 'PF7') ~ 'Soil_Env',
                                   Sample %in% c('FLU', 'HEX', 'OIL', 'PHE', 'PHE+FLU', 'PHE+HEX', 'PHE+OIL', 'PHE+PYR', 'PYR') ~ 'Water_MC',
                                   Sample %in% c('FAR_Island_C', 'NOR_Island_C') ~ 'Coral_Env',
                                   Sample %in% c('NOR_Island_W', 'PRI_Island_W', 'SAN_Island_W') ~ 'Water_Env')) %>%
    mutate(Group = case_when(grepl('_MC', Sample_Type) ~ 'Microcosms',
                             grepl('_Env', Sample_Type) ~ 'Environmental')) %>%
    select(Sample, Sample_Type, Group, all_of(col_tax))
  head(data_t)
  
  # Remove rows with zeros
  data_t = data_t %>%
    mutate(flag = rowSums(across(where(is.numeric)))) %>%
    filter(flag > 0) %>%
    select(-flag)
  head(data_t)
  
  write.table(data_t, txt_file_transpose, sep = '\t', quote = FALSE, row.names = FALSE)
}

file = 'data/data_is.txt'
columns = c('Crude Oil','PD5','PD6','PF7','FLU','HEX','OIL','PHE','PHE+FLU','PHE+HEX','PHE+OIL','PHE+PYR','PYR','NOR_Island_W','PRI_Island_W','SAN_Island_W','FAR_Island_C','NOR_Island_C')

get_tables(file = file, category = 'Phylum', type = 'all', columns = columns, transform_type = 'hellinger')
get_tables(file = file, category = 'Class', type = 'all', columns = columns, transform_type = 'hellinger')
get_tables(file = file, category = 'Order', type = 'all', columns = columns, transform_type = 'hellinger')
get_tables(file = file, category = 'Family', type = 'all', columns = columns, transform_type = 'hellinger')
get_tables(file = file, category = 'Genus', type = 'all', columns = columns, transform_type = 'hellinger')
get_tables(file = file, category = 'Species', type = 'all', columns = columns, transform_type = 'hellinger')
