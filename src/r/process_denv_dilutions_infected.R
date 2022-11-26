
process_denv_dilutions_infected <- function(filename_denv_dilutions) {
  
  df <- readxl::read_xlsx(filename_denv_dilutions) %>% 
    rename_all(tolower) %>% 
    rename(infection_status=`infection status`) %>% 
    mutate(is_infected=if_else(infection_status=="Positive",
                               1, 0))
  lookup <- tribble(
    ~dilution, ~dilution_val,
    "1:1", 1,
    "1:5", 5,
    "1:12", 12,
    "1:20", 20,
    "1:25", 25
  )
  
  sum_df <- df %>% 
    group_by(day, tissue, dilution) %>% 
    summarise(n_dissected=n(),
              n_infected=sum(is_infected)) %>% 
    ungroup() %>% 
    left_join(lookup) %>% 
    rename(dilutions=dilution_val) %>% 
    mutate(type="single") %>% 
    mutate(virus="DENV") %>% 
    mutate(refeed_time=NA)
  
  sum_df
}