# Useful functions for mangling abisko data frames and data tables.

# Given a data frame df that contains a sample name column, add month, year, site, core (e.g. 1,2,3) and core_name e.g. "201208_E1".
# 
# sample_col_name: the name of the column that is the sample name
# using00: Use true if the names are like "20120800_E1M", otherwise "201208_E1M" is expected. Actually the depth is not required.
# from_core_name: Derive info from the 'core_name' rather than 'sample'
#
add_core_name_metadata = function(df, sample_col_name='sample',using00=TRUE,from_core_name=FALSE){
 month_levels = c('05','06','07','08','09','10')
 year_levels = c('2010','2011','2012')
 if (from_core_name & sample_col_name=='sample'){
  sample_col_name = 'core_name'
 }
 sample_names = as.data.frame(df)[,sample_col_name]
 if (using00 & !from_core_name){
  df$site = factor(gsub(x=sample_names,  '^......00_(.).*', '\\1', perl=T), levels=c('P','S','E'))
  df$year = factor(gsub(x=sample_names,  '^(....)..00_.*', '\\1', perl=T), levels=year_levels)
  df$core = factor(gsub(x=sample_names, '^......00_.(.).*', '\\1', perl=T), levels=c('1','2','3'))
  df$month = factor(gsub(x=sample_names, '^....(..)00_.*', '\\1', perl=T), levels=month_levels)
 } else {
  df$site = factor(gsub(x=sample_names,  '^......_(.).*', '\\1', perl=T), levels=c('P','S','E'))
  df$year = factor(gsub(x=sample_names,  '^(....).._.*', '\\1', perl=T), levels=year_levels)
  df$core = factor(gsub(x=sample_names, '^......_.(.).*', '\\1', perl=T), levels=c('1','2','3'))
  df$month = factor(gsub(x=sample_names, '^....(..)_.*', '\\1', perl=T), levels=month_levels)
 }
 if (!from_core_name){
  df$month = factor(df$month, levels=unique(df$month))
 }
 unique_months = unique(df$month)
 new_month_levels = month_levels[month_levels %in% unique_months]
 df$month = factor(df$month, levels=new_month_levels)

 df$core_name = paste(df$year, df$month, '_', df$site, df$core, sep='')
 return(df)
}

# Like add_core_name_metadata, except also add the depth and site_ssplit ie. 'P','S_surface','S_anaerobic','E' based on the depth.
add_metadata = function(df, sample_col_name='sample',using00=TRUE){
 df = add_core_name_metadata(df, sample_col_name,using00)
 sample_names = as.data.frame(df)[,sample_col_name]
 if (using00){
  df$depth = factor(gsub(x=sample_names, '^......00_..(.).*', '\\1', perl=T), levels=c('S','M','D','X'))
 } else {
  df$depth = factor(gsub(x=sample_names, '^......_..(.).*', '\\1', perl=T), levels=c('S','M','D','X'))
 }
 df$site_ssplit = as.character(df$site)
 df$site_ssplit[df$site == 'S' & df$depth == 'S'] = 'S_surface'
 df$site_ssplit[df$site == 'S' & df$depth != 'S'] = 'S_anaerobic'
 df$site_ssplit = factor(df$site_ssplit, levels=c('P','S_surface','S_anaerobic','E'))
 
 return(df)
}

# Given a data frame/table with metadata, remove the may 2012 long permafrost cores and return.
remove_may2012_samples = function(df){
 dff = df[df$year != '2012' | df$month != '05',]
 dff$month = factor(dff$month, levels = unique(dff$month)) #this keeps the order, apparently
 return(dff)
}

# Add material and replicate columns the a df that has a 'sample' column.
add_ftms_metadata = function(df){
 df$material = factor(gsub(x=df$sample, '.*_(..)$', '\\1', perl=T), levels=c('PW','SP'))
 df$replicate = gsub(x=df$sample, '.*_(..)_..$', '\\1')=='R2'
 return(df)
}

# worker function for split_taxonomy
get_taxonomy_field = function(df, column, taxonomy_column){
  return(factor(gsub("^\\s+|\\s+$", "", as.character(cbind(lapply(strsplit(paste(df[[taxonomy_column]], ';',';',';',';',';',';',';',';'), ';'), '[[', column))))))
}

# Split up a semi-colon separated taxonomy field
split_taxonomy = function(df, has_root=T, taxonomy_field='taxonomy'){
  start_index = 1
  if (has_root){start_index = start_index+1}
  df$domain = get_taxonomy_field(df, start_index, taxonomy_field)
  df$phylum = get_taxonomy_field(df, start_index+1, taxonomy_field)
  df$class_name = get_taxonomy_field(df, start_index+2, taxonomy_field)
  df$order_name = get_taxonomy_field(df, start_index+3, taxonomy_field)
  df$family = get_taxonomy_field(df, start_index+4, taxonomy_field)
  df$genus = get_taxonomy_field(df, start_index+5, taxonomy_field)
  df$species = get_taxonomy_field(df, start_index+6, taxonomy_field)
  return(df)
}

