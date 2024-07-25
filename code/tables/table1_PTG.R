# Felzartamab Trial 
## Data pull

# Get data clinical 
library(haven)
# d <- read_sav(
#   here::here("data", "Generalfile_Felzartamab SPSS.sav")
# )
d <- read_sav("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/Generalfile_Felzartamab SPSS.sav")


# Exclude Patient STUDY_EVALUATION_ID 15 and 18 

d<- d |> dplyr::filter(STUDY_EVALUATION_ID != c(15, 18))

# Inspect dataset
library(tidyverse)
library(summarytools) # install.packages("summarytools")
library(tidyverse)
library(flextable)
#Font einstellen einmalig: 
set_flextable_defaults(font.family = "Arial", 
                       font.size = 12)
library(officer)
library(finalfit) #install.packages("finalfit")



# Select variables to be changed to factor
factor_variables <- Hmisc::.q(Female_gender, LD_Tx, Prior_Tx, PRA_anyType_above_equal10, PreTx_DSA, IndexBx_ABMR_BANFF19, IndexBx_Active_ABMR_Banff19, IndexBx_Chronic_active_ABMR_Banff19, 
                              IndexBx_C4d_positive_ABMR_Banff19, IndexBx_Borderline_Banff, Index_Early_stage_ABMR_Archetype, Index_Fully_developed_ABMR_Archetype, 
                              Index_Late_stage_ABMR_Archetype,Tac_Screening, Dual_Therapy_Screening, Triple_Therapy_Screening, DSA_HLA_class_I_Only_Screening, DSA_HLA_class_II_Only_Screening, 
                              DSA_HLA_class_I_and_II_Screening, HLA_DQ_DSA_Screening, Ethnic_origin_White, Ethnic_origin_Asian, ATG_Induction, IL2R_Induction, Tac_Initiation, CyA_Initiation, 
                              MPA_Initiation, Aza_Initiation, 
                              SOC_Infections_and_infestations_10021881, PT_Urinary_tract_infection_10046571, PT_Nasopharyngitis_10028810, PT_Bronchitis_10006451, PT_Pneumonia_bacterial_10060946, PT_Paronychia_10034016, 
                              PT_COVID_19_10084268, PT_Respiratory_syncytial_virus_infection_10061603, PT_Herpes_simplex_10019948, PT_Cytomegalovirus_viraemia_10058881, PT_Herpessimplexviraemia_10080365, 
                              PT_Clostridium_difficile_colitis_10009657, SOC_Gastrointestinal_disorders_10017947, PT_Nausea_10028813, PT_Vomiting_10047700, PT_Drymouth_10013781, PT_Toothache_10044055, PT_Tooth_abscess_10044016, PT_Diarrhoea_10012735, 
                              PT_Abdominal_distension_10000060, PT_Abdominaldiscomfort_10000059, PT_Abdominal_pain_upper_10000087, PT_Abdominal_pain_10000081, PT_Umbilical_hernia_10045458, SOC_General_disorders_and_administration_siteconditions_10018065, 
                              PT_Infusion_related_reaction_10051792, PT_Influenza_10022000,
                              PT_Oedema_peripheral_10030124, PT_Pyrexia_10037660, PT_Chills_10008531, PT_Fatigue_10016256, PT_Asthenia_10003549, PT_Malaise_10025482, PT_Swelling_10042674, PT_Chest_pain_10008479, 
                              SOC_Skin_and_subcutaneous_tissue_disorders_10040785,  PT_Pruritus_10037087, PT_Subcutaneous_hematoma_10042346, PT_Rash_10037844, PT_Tinea_versicolour_10056131, SOC_Vascular_disorders_10047065, 
                              PT_Hypotension_10021097, PT_Hypertension_10020772, SOC_Musculoskeletal_and_connective_tissue_disorders_10028395, PT_Muscle_spasms_10028334, 
                              PT_Back_pain_10003988, PT_Pain_in_extremity_10033425, PT_Arthralgia_10003239, PT_Sciatica_10039674, PT_Gouty_arthritis_10018634, PT_Neck_pain_10028836, SOC_Nervous_system_disorders_10029205, 
                              PT_Headache_10019211, SOC_Reproductive_system_and_breast_disorders_10038604, 
                              PT_Polymenorrhoea_10036086, SOC_Respiratory_thoracic_and_mediastinal_disorders_10038738,  PT_Cough_10011224, PT_Throat_irritation_10043521, PT_Dyspnoea_10013968, PT_Rhinitis_10039083, PT_Oropharyngeal_pain_10068319,
                              SOC_Blood_and_lymphatic_system_disorders_10005329,  PT_Leukopenia_10024384, 
                              PT_Anaemia_10002034, PT_Nephrogenic_anaemia_10058116, SOC_Injury_poisoning_and_procedural_complications_10022117,  PT_Radius_fracture_10037802, PT_Muscle_injury_10028314, PT_Contusion_10050584, PT_Skin_injury_10061364,
                              SOC_Cardiac_disorders_10007541,  PT_Palpitations_10033557, SOC_Renal_and_urinary_disorders_10038359, 
                              PT_Acute_kidney_injury_10069339, PT_Dysuria_10013990, PT_Renal_pain_10038490, SOC_Eye_disorders_10015919, PT_Photophobia_10034960, PT_Keratitis_10023332, SOC_Metabolism_and_nutrition_disorders_10027433, 
                              PT_Metabolic_acidosis_10027417, PT_Hypocalcemia_10020949, 
                              PT_Hyponatraemia_10021036, PeriTx_IA, CDC_PRA_above_equal10)



# change variables to factor 
d[,factor_variables] <- lapply(d[,factor_variables], factor)

# Create variable MFI >5000
### Variable DSA MFI of immunodominant DSA >5000

d<- d %>% mutate(mfi_immunodominant = case_when(DSA_highest_Category_Screening > 3~ 1, 
                                                TRUE ~ 0)) %>% 
  mutate(mfi_immunodominant = factor(mfi_immunodominant))


# Create table 1 
# Select variables 

# Copy cont. variables to be calculated as mean 
d<- d |> mutate(age_at_tx_mean = Age_at_Tx, 
                donor_age_mean = Donor_Age, 
                age_at_screening_mean = Age_at_Screening, 
                screening_egfr_mean = Screening_GFR, 
                screening_prot_krea_mean = Screening_Prot_Krea)

# Calculate median age with range for all age related variables 
age_tx_summary <- tbl1_data %>%
  group_by(Felzartamab) %>%
  summarise(Age_at_Tx_summary = paste0(round(median(Age_at_Tx, na.rm = TRUE)), " (", 
                                       round(min(Age_at_Tx, na.rm = TRUE)), "-", 
                                       round(max(Age_at_Tx, na.rm = TRUE)), ")")) %>%
  pivot_wider(names_from = Felzartamab, values_from = Age_at_Tx_summary)

donor_age_summary <- tbl1_data %>%
  group_by(Felzartamab) %>%
  summarise(Donor_Age_summary = paste0(round(median(Donor_Age, na.rm = TRUE)), " (", 
                                       round(min(Donor_Age, na.rm = TRUE)), "-", 
                                       round(max(Donor_Age, na.rm = TRUE)), ")")) %>%
  pivot_wider(names_from = Felzartamab, values_from = Donor_Age_summary)

age_screening_summary <- tbl1_data %>%
  group_by(Felzartamab) %>%
  summarise(Age_at_Screening_summary = paste0(round(median(Age_at_Screening, na.rm = TRUE)), " (", 
                                              round(min(Age_at_Screening, na.rm = TRUE)), "-", 
                                              round(max(Age_at_Screening, na.rm = TRUE)), ")")) %>%
  pivot_wider(names_from = Felzartamab, values_from = Age_at_Screening_summary)



# Calculate for total 
total_summary <- tbl1_data %>%
  summarise(
    Age_at_Tx_summary = paste0(round(median(Age_at_Tx, na.rm = TRUE)), " (", 
                               round(min(Age_at_Tx, na.rm = TRUE)), "-", 
                               round(max(Age_at_Tx, na.rm = TRUE)), ")"),
    Donor_Age_summary = paste0(round(median(Donor_Age, na.rm = TRUE)), " (", 
                               round(min(Donor_Age, na.rm = TRUE)), "-", 
                               round(max(Donor_Age, na.rm = TRUE)), ")"),
    Age_at_Screening_summary = paste0(round(median(Age_at_Screening, na.rm = TRUE)), " (", 
                                      round(min(Age_at_Screening, na.rm = TRUE)), "-", 
                                      round(max(Age_at_Screening, na.rm = TRUE)), ")")
  )

# make tibble to be add to table 
age_tx_tibble <- tibble(
  label = "Median recipient age (range) - yr",
  `0` = age_tx_summary$`0`,
  `1` = age_tx_summary$`1`,
  Total = total_summary$Age_at_Tx_summary
)

donor_age_tibble <- tibble(
  label = "Median donor age (range) – yr",
  `0` = donor_age_summary$`0`,
  `1` = donor_age_summary$`1`,
  Total = total_summary$Donor_Age_summary
)

age_screening_tibble <- tibble(
  label = "Median age of study patients (range) – yr",
  `0` = age_screening_summary$`0`,
  `1` = age_screening_summary$`1`,
  Total = total_summary$Age_at_Screening_summary
)




tbl1_data <- d %>% select(Female_gender, Age_at_Tx, age_at_tx_mean, LD_Tx, Donor_Age, donor_age_mean, MM_ABDR, 
                          Age_at_Screening, age_at_screening_mean, Years_to_ScreeningVisit, Screening_GFR, screening_egfr_mean, Screening_Prot_Krea, screening_prot_krea_mean,
                          Dual_Therapy_Screening, Triple_Therapy_Screening,
                          IndexBx_Active_ABMR_Banff19 ,IndexBx_Chronic_active_ABMR_Banff19, 
                          IndexBx_Borderline_Banff,IndexBx_TCMR_Banff,  DSA_HLA_class_I_Only_Screening, DSA_HLA_class_II_Only_Screening, DSA_HLA_class_I_and_II_Screening,
                          HLA_DQ_DSA_Screening, mfi_immunodominant, DSA_Number_Screening,
                          Felzartamab) 



explanatory <- Hmisc::.q(Female_gender, age_at_tx_mean, LD_Tx, donor_age_mean, MM_ABDR, 
                         age_at_screening_mean, Years_to_ScreeningVisit, Screening_GFR, screening_egfr_mean, Screening_Prot_Krea, screening_prot_krea_mean,
                         IndexBx_Active_ABMR_Banff19 ,IndexBx_Chronic_active_ABMR_Banff19, 
                         IndexBx_Borderline_Banff, DSA_HLA_class_I_Only_Screening, DSA_HLA_class_II_Only_Screening, 
                         DSA_HLA_class_I_and_II_Screening, HLA_DQ_DSA_Screening, mfi_immunodominant, DSA_Number_Screening )    





dependent <- ("Felzartamab")  

# Create table 1 

tbl1<- tbl1_data %>% mutate(Felzartamab = factor(Felzartamab), 
                            Female_gender = ff_label(Female_gender, "Female sex – no. (%)"), 
                            #Age_at_Tx = ff_label(Age_at_Tx, "Median recipient age (IQR) – yr"),
                            age_at_tx_mean = ff_label(age_at_tx_mean, "Mean recipient age (SD) - yr"), 
                            LD_Tx = ff_label(LD_Tx , "Living donor – no. (%)"), 
                            # Donor_Age = ff_label(Donor_Age, "Median donor age (IQR) – yr"), 
                            donor_age_mean = ff_label(donor_age_mean,"Mean donor age (SD) – yr" ),
                            MM_ABDR = ff_label(MM_ABDR, "Median HLA (A, B, DR) mismatch (IQR)"), 
                            #Age_at_Screening = ff_label(Age_at_Screening, "Median age of study patients (IQR) – yr"), 
                            age_at_screening_mean = ff_label(age_at_screening_mean, "Mean age of study patients (SD) – yr"), 
                            Years_to_ScreeningVisit = ff_label(Years_to_ScreeningVisit, "Median time to inclusion in the trial (IQR) – yr"), 
                            Screening_GFR= ff_label(Screening_GFR, "Median eGFR (IQR) – mL/min/1.73 m2"), 
                            screening_egfr_mean= ff_label(screening_egfr_mean, "Mean eGFR (SD) – mL/min/1.73 m2"), 
                            Screening_Prot_Krea = ff_label(Screening_Prot_Krea, "Median protein/creatinine ratio (IQR) –  mg/g"), 
                            screening_prot_krea_mean = ff_label(screening_prot_krea_mean, "Mean protein/creatinine ratio (SD) –  mg/g"), 
                            IndexBx_Active_ABMR_Banff19 = ff_label(IndexBx_Active_ABMR_Banff19, "Active ABMR"), 
                            IndexBx_Chronic_active_ABMR_Banff19 = ff_label(IndexBx_Chronic_active_ABMR_Banff19, "Chronic active ABMR"), 
                            IndexBx_Borderline_Banff = ff_label(IndexBx_Borderline_Banff, "Additional borderline lesion"), 
                            DSA_HLA_class_I_Only_Screening = ff_label(DSA_HLA_class_I_Only_Screening, "HLA class I DSA only – no. (%)"), 
                            DSA_HLA_class_II_Only_Screening = ff_label(DSA_HLA_class_II_Only_Screening, "HLA class II DSA only – no. (%)"), 
                            DSA_HLA_class_I_and_II_Screening = ff_label(DSA_HLA_class_I_and_II_Screening, "HLA class I and II DSA – no. (%)"), 
                            # Dual_Therapy_Screening = ff_label(Dual_Therapy_Screening, "Dual immunosuppression – no. (%)"), 
                            #Triple_Therapy_Screening = ff_label(Triple_Therapy_Screening, "Triple immunosuppression – no. (%)"), 
                            HLA_DQ_DSA_Screening = ff_label(HLA_DQ_DSA_Screening, "Anti-DQ DSA – no. (%)"), 
                            mfi_immunodominant = ff_label(mfi_immunodominant, "Peak MFI of DSA >10,000 – no. (%)"), 
                            DSA_Number_Screening = ff_label(DSA_Number_Screening, "Median DSA (IQR) – n")) |> 
  summary_factorlist(dependent, explanatory, 
                     total_col = T, 
                     #cont = "median", 
                     na_include = F,
                     p=F,
                     cont_cut = F, 
                     digits = c(0, 0, 3, 1, 0)  , 
                     cont_nonpara = c(5,7,8,10,13,20)) %>% 
  ff_remove_ref(only_binary = F) |> 
  select(-levels)


# Add the extra rows
tbl1<- 
  tbl1 %>% 
  tibble::add_row(label = "Recorded at transplantation", .after = 0) %>% 
  tibble::add_row(age_tx_tibble, .after= 2) |> 
  tibble::add_row(donor_age_tibble, .after= 5) |> 
  tibble::add_row(label = "Recorded at trial inclusion", .after = 8) %>%
  tibble::add_row(age_screening_tibble, .after= 9) |> 
  tibble::add_row(label = "Banff 2019 ABMR phenotypes (baseline biopsies) – no. (%)", .after = 16) %>%
  tibble::add_row(label = "DSA characteristics (screening visit)", .after = 20) 


# Save table 1 

flextable::qflextable(tbl1) %>% 
  set_table_properties(layout = "autofit") %>%
  set_header_labels(label = "Variables", Total = "Total (N=20)", 
                    "0" = "Placebo (n=10)", "1" = "Felzartamab (n=10)") %>% 
  flextable::bold(i = NULL, part = "header") %>% 
  flextable::bold(i = c(1,9), j=1) %>%
  flextable::padding(i= c(18:20, 22:27), j=1, padding.left=20) %>% 
  set_caption(caption="Table 1. Demographic and Characteristics of the Participants at Baseline.") %>%
  add_footer_lines("ABMR, antibody-mediated rejection; DSA, donor-specific antibody; eGFR, estimated glomerular filtration rate; HLA, human leukocyte antigen; IQR, interquartile range; MFI, mean fluorescence intensity; SD, standard deviation") %>%
  hline_bottom(part = "body", border = fp_border(color = "black", width = 1)) |> 
  hline_bottom(part = "header", border = fp_border(color = "black", width = 1)) |> 
  hline_top(part = "header", border = fp_border(color = "black", width = 1)) |> 
  flextable::save_as_docx(path =  here::here("tables", "table1.docx"), 
                          pr_section = prop_section(page_size = page_size(orient = "landscape"))) 




