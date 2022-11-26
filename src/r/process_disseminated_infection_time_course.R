
process_disseminated_infection_time_course <- function(
    filename_disseminated_time_course) {
  
  df <- readxl::read_xlsx(filename_disseminated_time_course) %>% 
    rename_all(tolower) %>% 
    rename(infection_status=`infection status`) %>%
    rename(dissemination_status=`dissemination status`) %>%
    rename(feeding_status=`feeding status`) %>% 
    mutate(feeding_status=if_else(feeding_status=="SF",
                                  "single", "double")) %>% 
    rename(midgut=infection_status) %>% 
    rename(legs=dissemination_status) %>% 
    mutate(midgut=if_else(midgut=="Positive", 1, 0),
           legs=if_else(legs=="Positive", 1, 0))
  
  sum_df <- df %>% 
    group_by(day, feeding_status) %>% 
    summarise(n=n(),
              midgut=sum(midgut),
              legs=sum(legs)) %>% 
    pivot_longer(c(midgut, legs)) %>% 
    rename(tissue=name) %>% 
    rename(n_infected=value) %>% 
    mutate(dilution="1:5") %>% 
    mutate(refeed_time=ifelse(feeding_status=="single", NA, 3)) %>% 
    mutate(dilutions=5) %>% 
    rename(type=feeding_status) %>% 
    rename(n_dissected=n) %>% 
    mutate(virus="DENV") %>% 
    ungroup()
  
  sum_df
}