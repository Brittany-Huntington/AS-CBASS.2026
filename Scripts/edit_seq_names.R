library(here)

library(dplyr)
library(stringr)
df<-read.csv(here( "Data/seq_ids.csv"))


clean_sample_ids <- function(ids) {
  ids %>%
    # 1. Convert everything to lowercase
    tolower() %>%
    
    # 2. Drop "tut" and any immediate separator following it (e.g., "tut_", "tut-", "tut ")
    str_replace_all("^tut[-_\\s]*", "") %>%
    
    # 2. Drop "t" and any immediate separator following it (e.g., "tut_", "tut-", "tut ")
    str_replace_all("^t[-_\\s]*", "") %>%
    
    # 3. Handle cases where site/time prefix has a separator between letters and numbers 
    # e.g., "t 1" or "t-1" or "occ 12" -> "t1", "occ12"
    str_replace_all("^([a-z]+)[-_\\s]+(\\d+)", "\\1\\2") %>%
    
    # 4. Standardize remaining separators (spaces, hyphens, multiple underscores) into a single "_"
    str_replace_all("[-_\\s]+", "_") %>%
    
    # 5. Trim any leading/trailing underscores
    str_trim() %>%
    str_replace_all("^_+|_+$", "")
}


df <- df %>%
  mutate(cleaned_id = clean_sample_ids(Sample))

print(df)
write.csv(df, here ("Data", "seq_ids_clean.csv"))
