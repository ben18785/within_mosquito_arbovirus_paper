
process_damage_data <- function(filename_damage) {
  
  lookup <- tribble(
    ~name, ~time,
    "Unfed", 10, # assume unfed ~ mosquitoes recovery after 10 days
    "15 mpbm", 15 / (24 * 60),
    "24 hpbm", 1,
    "36 hpbm", 1.5,
    "48 hpbm", 2,
    "72 hpbm", 3,
    "96 hpbm", 4
  )
  
  df <- readxl::read_xlsx(filename_damage) %>% 
    mutate(iterate=seq_along(Unfed)) %>% 
    mutate_all(~gsub("\\*", "", .)) %>% 
    pivot_longer(-iterate) %>% 
    left_join(lookup) %>% 
    filter(!is.na(value)) %>% 
    mutate(value=as.numeric(value)) %>% 
    rename(chp=value)
  
  df
}