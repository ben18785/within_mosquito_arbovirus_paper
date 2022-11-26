
process_experimental_data <- function(
    filename_midgut, filename_legs) {
  
  midgut <- readxl::read_xlsx(filename_midgut)
  legs <- readxl::read_xlsx(filename_legs)
  
  both <- midgut %>%
    bind_rows(legs) %>% 
    rename_all(tolower) %>% 
    rename_all(~gsub(" ", "_", .)) %>% 
    mutate(sample=as.factor(sample)) %>% 
    mutate(sample=fct_relevel(sample, "1:1", "1:5", "1:12")) %>% 
    mutate(tissue=as.factor(tissue)) %>% 
    mutate(tissue=fct_relevel(tissue, "midgut", "legs")) %>% 
    mutate(day_individual_sample=paste0(day, "-", individual, "-", sample)) %>% 
    mutate(dilutions = ifelse(sample == "1:1", 1,
                              ifelse(sample == "1:5", 5,
                                     ifelse(sample == "1:12", 12, NA))))
  both
}