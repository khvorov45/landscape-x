library(tidyverse)

hanam <- read_csv("data-raw/hanam/HI_hanam.csv", col_types = cols())

read_fasta_as_table <- function(filename) {
  fasta_string <- read_file(filename)
  fasta_string_split <- str_split(fasta_string, ">")[[1]] %>% `[`(. != "")
  fasta_string_meta_seq_split <- str_split(fasta_string_split, "\r\n|\n", n = 2)
  fasta_string_meta_seq_split_only2 <- fasta_string_meta_seq_split %>% `[`(map_lgl(., ~ length(.x) == 2))
  tibble(
    meta_og = fasta_string_meta_seq_split_only2 %>% map_chr(~ .x[[1]]),
    seq_og = fasta_string_meta_seq_split_only2 %>% map_chr(~ .x[[2]]),
  )
}

hanam_seq <- read_fasta_as_table("data-raw/hanam/LScape virus HA protein_no signal.fas") %>%
  bind_rows(read_fasta_as_table("data-raw/hanam/Kansas_14_2017eggNoSignal.fap")) %>%
  bind_rows(read_fasta_as_table("data-raw/hanam/Switz_8060_2017egg_3c2a2_HA.fas"))

hanam_seq %>%
  select(meta_og) %>% distinct() %>% arrange(meta_og) %>% print(n = 100)

hanam_viruses <- hanam %>% 
  mutate(virus_eggcell = if_else(Egg_Cell == 1, "Cell", "Egg")) %>%
  select(
    virus_short = Short_Name, virus_abbrv = Virus_Abbrv, 
    virus_number = virus, virus_eggcell, virus_cluster = Cluster, virus_year = Year,
    virus_clade = Clade, virus_clade_code = CladeCode
  ) %>% 
  distinct() %>%
  left_join(
    hanam_seq %>%
      mutate(virus_short = recode(
        meta_og,
        "Bilthoven_16190_1968" = "Bilt/16190/68",
        "Bilthoven_21793_1972" = "Bilt/21793/72",
        "Bilthoven_1761_1976" = "Bilt/1761/76",
        "Bilthoven2271_1976" = "Bilt/2271/76",
        "Brisbane_10_2007" = "Bris/10/07",
        "Brisbane_60_2018" = "Brisbane/60/18",
        "Fujian_411_2002" = "Fujian/411/02",
        "HongKong_4801_2014Egg" = "H_Kong/4801/14e",
        "Hanam_EL134_273S03_08_08_s1_ls" = "Hanam/EL134/08",
        "Hanam_EL201_145S01_05_09_s2_ls_Vax" = "Hanam/201/09",
        "Hanam_EL444_181S01_09_10_s4_ls_Vax" = "Hanam/444/10",
        "Hanam_EL14437_062S03_05_14_s9" = "Hanam/14437/14",
        "Kansas_14_17" = "Kansas/14/17",
        "Kansas_14_2017egg.pro Translate DNA Sequence Untitled Seq #1(1,1714)" = "Kansas/14/17e",
        "Michigan_15_2014" = "Michigan/15/14",
        "NewCaledonia_104_2014" = "N_Caled/104/14",
        "Newcastle_30_2016" = "N_Castle/30/16",
        "Netherlands_233_1982" = "NethLand/233/82",
        "Netherlands_620_1989" = "NethLand/620/89",
        "Netherlands_823_1992" = "NethLand/823/92",
        "Netherlands_178_1995" = "NethLand/178/95",
        "Netherlands_301_1999" = "NethLand/301/99",
        "Netherlands_179_1993" = "NethLand/179/93",
        "NewYork_55_2004Egg" = "N_York/55/04e",
        "Perth_16_2009cell" = "Perth/16/09",
        "Perth_16_2009Egg" = "Perth/16/09e",
        "Philippines_472_2002" = "Philippine/472/02",
        "Philippines_2_1982" = "Philippine/2/82",
        "Switzerland_9715293_2013" = "Switz/9715293/13",
        "Switzerland_9715293_2013Egg" = "Switz/9715293/13e",
        "Sw/8060/17" = "Switz/8060/17",
        "Tasmania_1_1997" = "Tasmania/1/97",
        "Thailand_409_2005" = "Thai/409/05",
        "Townsville_2_1999" = "Townsville/2/99",
        "Texas_50_2012" = "Texas/50/12",
        "Texas_50_2012Egg" = "Texas/50/12e",
        "Uruguay_716_2007Egg" = "Urug/716/07e",
        "Victoria_511_2004" = "Vic/511/04",
        "Victoria_361_11" = "Vic/361/11",
        "Victoria_361_2011Egg" = "Vic/361/11e",
        "Wisconsin_67_2005egg" = "Wisc/67/05e",
      )),
    "virus_short"
  )

hanam_viruses %>%
  mutate(seq_aligned = align_mafft(seq_og))


