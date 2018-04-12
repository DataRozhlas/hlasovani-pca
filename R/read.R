library(readr)

snemovna2017 <- read_delim("../data/raw/hl-2017ps/hl2017h1.unl", "|", col_names=c("voter_id", "vote_event_id", "option"))
