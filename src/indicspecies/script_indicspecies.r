library(indicspecies)
library(tibble)

get_table_network = function(category, type, transform_type) {
  file = paste0('data_is_category/table_', category, '_transpose_', transform_type, '_', type, '.txt')
  txt_file = paste0('data_network/table_network_', category, '_', transform_type, '_', type, '.txt')
  dir.create('data_network', showWarnings = FALSE)
  
  data = read.delim(file, sep = '\t', check.names = FALSE)
  head(data)
  
  abundance = data[,4:ncol(data)]
  
  cluster_type = data$Sample_Type
  # cluster_type = data$Group
  
  phi = multipatt(x = abundance, cluster = cluster_type, func = 'r.g', control = how(nperm = 9999))
  
  data_summary = phi$sign
  head(data_summary)
  
  data_summary = data_summary %>%
    rownames_to_column('Taxon') %>%
    rename('Soil_Env' = 's.Soil_Env',
           'Water_Env' = 's.Water_Env',
           'Coral_Env' = 's.Coral_Env')
  
  if (category != 'family' & category != 'genus' & category != 'species') {
    data_summary = data_summary %>%
      rename('Soil_MC' = 's.Soil_MC')
  }
  if (category != 'genus' & category != 'species') {
    data_summary = data_summary %>%
      rename('Water_MC' = 's.Water_MC')
  }
  head(data_summary)
  
  add_row_by_column = function(df, row, column) {
    df_row = data.frame()
    if (df[row, column] == 1) {
      df_row = data.frame(Source = column,
                          Target = df[row, 'Taxon'],
                          Edge_Weight = df[row, 'stat'],
                          p_value = df[row, 'p.value'])
    }
    return(df_row)
  }
  
  data_network = data.frame()
  for (row in 1:nrow(data_summary)) {
    data_network = rbind(data_network, add_row_by_column(df = data_summary, row = row, column = 'Soil_Env'))
    data_network = rbind(data_network, add_row_by_column(df = data_summary, row = row, column = 'Water_Env'))
    data_network = rbind(data_network, add_row_by_column(df = data_summary, row = row, column = 'Coral_Env'))
    
    if (category != 'family' & category != 'genus' & category != 'species') {
      data_network = rbind(data_network, add_row_by_column(df = data_summary, row = row, column = 'Soil_MC'))
    }
    if (category != 'genus' & category != 'species') {
      data_network = rbind(data_network, add_row_by_column(df = data_summary, row = row, column = 'Water_MC'))
    }
  }
  
  data_network = data_network %>%
    filter(p_value <= 0.05) %>%
    mutate(Node_Color = Source) %>%
    mutate(Edge_Color = Source)
  head(data_network)
  
  data_means = data %>%
    select(-Sample, -Group) %>%
    group_by(Sample_Type) %>%
    summarise(across(everything(), ~ mean(.)))
  data_means = as.data.frame(data_means)
  head(data_means)
  
  data_means = data_means %>%
    column_to_rownames('Sample_Type')
  head(data_means)
  
  data_network_plus = data.frame()
  for (row in 1:nrow(data_network)) {
    row_base = data_network[row,]
    Node_Degree = data_means[data_network[row, 1], data_network[row, 2]] * 10
    row_plus = cbind(row_base, Node_Degree)
    
    data_network_plus = rbind(data_network_plus, row_plus)
  }
  
  write.table(data_network_plus, txt_file, sep = '\t', quote = FALSE, row.names = FALSE)
}

get_table_network(category = 'phylum', type = 'all', transform_type = 'hellinger')
get_table_network(category = 'class', type = 'all', transform_type = 'hellinger')
get_table_network(category = 'order', type = 'all', transform_type = 'hellinger')
get_table_network(category = 'family', type = 'all', transform_type = 'hellinger')
get_table_network(category = 'genus', type = 'all', transform_type = 'hellinger')
get_table_network(category = 'species', type = 'all', transform_type = 'hellinger')
