### Load Libraries =================================
library(R0)
library(tidyverse)
library(zoo)
library(lubridate)
library(ggthemes)
library(ggdark)
library(scales)



### Custom functions =============================
## not in a vector
`%nin%` = Negate(`%in%`)



## Custom function for R0 estimation under gamma distributed generation times
get_R0 <- function(mean_gi = mean_gi, sd_gi = sd_gi, rho = rho) {
  exp(((mean_gi ^ 2) * log(1 + ((sd_gi ^ 2) * rho / mean_gi))) / (sd_gi ^ 2))
}



### Download Dataset ===============================
download.file(url = "https://api.covid19india.org/csv/latest/districts.csv",
              destfile = "district.csv")

### Load dataset ===================================
district <- read_csv("district.csv")



### Serial Interval Estimates =====================
sample_si <- 3924
mean_si <- 5.17  # days
lci_si <- 4.89   # days
uci_si <- 5.45   # days
sd_si <- (mean_si - lci_si)/ 1.96



### Generation Interval ==========================
set.seed(2021)
gen_time <- generation.time(type = "gamma", val = c(mean_si, sd_si))



### Time series for India ========================

# Preparing India Specific Dataset
india <- district %>% 
  group_by(Date) %>% 
  summarise(Confirmed = sum(Confirmed, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(Incidence = Confirmed - lag(Confirmed, orger_by = Date)) %>% 
  dplyr::select(-Confirmed) %>% 
  filter(!is.na(Incidence))

est_len <- length(india$Date) - 1

# Estimating Rt for each day
est_r0 <- est.R0.TD(epid = india$Incidence, gen_time, t = india$Date, begin = 1, end = est_len, nsim = 1000)

# Teasing out Rt values from the Rt object
Rt <- est_r0[["R"]]

Date <- as.Date(names(Rt))

names(Rt_india) <- NULL
# A tibble with Rt and dates
india <- as_tibble(data.frame(Date, Rt))


# Smoothed Rt (7-day SMA) and final dataset
india <- india %>% 
  mutate(smR = rollmean(Rt, k = 7, align = "right", fill = NA)) %>% 
  filter(Date >= dmy("15-05-2020")) %>% 
  dplyr::select(Date, smR)

# Writing the dataset for website
write.csv(india, "India.csv", row.names = FALSE)



### Time series for States ========================

# Preparing dataset for states
states <- district %>% 
  group_by(Date, State) %>% 
  summarise(Confirmed = sum(Confirmed, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(State) %>% 
  mutate(Incidence = Confirmed - lag(Confirmed, order_by = Date)) %>% 
  ungroup() %>% 
  dplyr::select(-Confirmed)


# Taking number of records for each state =============
x <- states %>% 
  filter(!is.na(Incidence))
table(x$State)


# Vector of state names ============================
state_names <- unique(states$State)


# Number of states included in the analysis ========
nStates <- length(state_names)

# Storage dataframe ========================
State <- "Placeholder"
Date <- as.Date("2016-10-10")
smR <- 1.00

perm_set <- as_tibble(data.frame(State = State, Date = Date, smR = smR))


# For loop for state-level analysis ==============
for(s in 1:nStates) {
  
# Dataset for one state
temp_state <- states %>% 
  filter(State == state_names[s]) %>% 
  filter(!is.na(Incidence))

# First day of case report for that state  
first_outbreak <- min(temp_state$Date[temp_state$Incidence > 0])

temp_state <- temp_state %>% 
  filter(Date >= first_outbreak)

# Change all zeros to 1
temp_state$Incidence[temp_state$Incidence <= 0] <- 1
temp_state$Incidence[is.na(temp_state$Incidence)] <- 1

# Duration of Rt estimation
est_len <- length(temp_state$Date) - 1

# Rt estimation
a <- Sys.time()

est_r0 <- est.R0.TD(epid = temp_state$Incidence, gen_time, t = temp_state$Date, begin = 1, end = est_len, nsim = 1000)

print(Sys.time() - a)


# Consolidating the numbers
Rt <- est_r0[["R"]]
Date <- as.Date(names(Rt))

names(Rt) <- NULL # removing metadata

# Preparing data frame
temp_set <- as_tibble(data.frame(Date, Rt))

temp_set <- temp_set %>% 
  mutate(smR = rollmean(Rt, k = 7, align = "right", fill = NA),
         State = state_names[s]) %>% 
  dplyr::select(State, Date, smR)

perm_set <- bind_rows(perm_set, temp_set)

print(paste0("Estimation completed for ", state_names[s], ". ", s, " of ",  nStates, " completed!"))

}

# Preparing final dataset ===================
perm_set <- perm_set %>% 
  filter(!is.na(smR)) %>% 
  rename(Region = State,
         Rt = smR)

# Write out CSV ====================
write.csv(perm_set, "States.csv", row.names = FALSE)





### Estimation for Districts ==============================

# Vector of noise in districts =====================
noise_districts <- c("Unknown", "Other State", "Other Region", "State Pool", "Italians", "Evacuees", "Railway Quarantine", "Foreign Evacuees", "Others", "Airport Quarantine", "Capital Complex")

# Prepare dataset ===================================
districts <- district %>% 
  filter(District %nin% noise_districts) %>%
  unite("id", c(State, District), sep = "_", remove = FALSE) %>% 
  group_by(id) %>% 
  mutate(Incidence = Confirmed - lag(Confirmed, order_by = Date)) %>% 
  filter(!is.na(Incidence)) %>% 
  dplyr::select(id, State, District, Date, Incidence)

# Vector of state names ============================
district_names <- unique(districts$id)


# Number of states included in the analysis ========
nDistricts <- length(district_names)

# Storage dataframe ========================
State <- "Placeholder"
District <- "Placeholder"
Date <- as.Date("2016-10-10")
smR <- 1.00

dist_set <- as_tibble(data.frame(State = State, District = District, Date = Date, smR = smR))


# For loop for state-level analysis ==============
for(d in 1:nDistricts) {
  
  # Dataset for one district
  temp_district <- districts %>% 
    filter(id == district_names[d]) %>% 
    filter(!is.na(Incidence))
  
  # State and District names
  sName <- temp_district$State[1]
  dName <- temp_district$District[1]
  
  # First day of case report for that district  
  first_outbreak <- min(temp_district$Date[temp_district$Incidence > 0])
  
  temp_district <- temp_district %>% 
    filter(Date >= first_outbreak)
  
  # Change all zeros to 1
  temp_district$Incidence[temp_district$Incidence <= 0] <- 1
  temp_district$Incidence[is.na(temp_district$Incidence)] <- 1
  
  # Duration of Rt estimation
  est_len <- length(temp_district$Date) - 1
  
  
  # Rt estimation
  a <- Sys.time()
  est_r0 <- est.R0.TD(epid = temp_district$Incidence, gen_time, t = temp_district$Date, begin = 1, end = est_len, nsim = 500)
  print(Sys.time() - a)
  
  # Consolidating the numbers
  Rt <- est_r0[["R"]]
  Date <- as.Date(names(Rt))
  
  names(Rt) <- NULL # removing metadata
  
  # Preparing data frame
  temp_set <- as_tibble(data.frame(Date, Rt))
  
  temp_set <- temp_set %>% 
    mutate(smR = rollmean(Rt, k = 7, align = "right", fill = NA),
           State = sName,
           District = dName) %>% 
    dplyr::select(State, District, Date, smR)
  
  dist_set <- bind_rows(dist_set, temp_set)
  
  print(paste0("Estimation completed for ", dName, " in ", sName, ". ", d, " of ",  nDistricts, " completed!"))
  
}

# Preparing final dataset ===================
dist_set <- dist_set %>% 
  filter(!is.na(smR))

# Write out CSV ====================
write.csv(dist_set, "Districts.csv", row.names = FALSE)



########





### Filter out Kerala ==============================
kerala <- district %>% 
  filter(State == "Kerala") %>% 
  filter(District %nin% noise_districts) %>% 
  group_by(District) %>% 
  mutate(Incidence = Confirmed - lag(Confirmed, order_by = Date)) %>% 
  dplyr::select(Date, District, Confirmed, Incidence) %>% 
  unite("id", c(Date, District), sep = " ", remove = FALSE)

### Start Date =============================================
begin_date <- ymd("2021-03-01") - 10

### Duration of estimation period in number of days ===========================
est_length <- as.numeric(Sys.Date() - begin_date - 1) #in days

### Vector of district names =====================================
district_names <- unique(kerala$District)

### Number of districts ==========================================
nDistricts <- length(district_names)

### Data frame for storing data from loop =======================
Rt <- matrix(NA, nrow = est_length, ncol = nDistricts)

### For Loop =======================================
for(d in 1:nDistricts){
dist_name <- district_names[d]

jilla <- kerala %>% 
  filter(District == dist_name) %>% 
  filter(Date >= begin_date) # Facilitating for smoothing. 

est_r0 <- est.R0.TD(epid = jilla$Incidence, gen_time, t = jilla$Date, begin = 1, end = est_length)

x <- est_r0[["R"]]
names(x) <- NULL

Rt[, d] <-  x

p <- length(x)
print(paste0(p, " days estimated for ", d, " of ", nDistricts, " :", district_names[d], "! Estimation Completed"))
rm(p, dist_name, x, est_r0, jilla)
}

### Preparing Dataset =================================
colnames(Rt) <- district_names

### Dataset and backup dataset ========================
Rt <- as_tibble(as.data.frame(Rt))
Rt_2 <-  Rt


# Rt <- Rt_2

### Preparing dataset for graphing ========================
Rt <- Rt %>% 
  mutate(Date = seq(begin_date, length.out = est_length, by = 1)) %>% 
  pivot_longer(cols = -Date, names_to = "District", values_to = "R") %>%
  group_by(District) %>% 
  mutate(smR = rollmean(R, k = 7, align = "right", fill = NA)) %>% 
  filter(Date >= dmy("15-03-2021"))

### Upper bound for y-axis for plotting ========================
y <- round(max(Rt$smR), 0)
upper_y <- as.integer(if(y > max(Rt$smR)){y} else {y + 1})


### Plotting district plots ============================


for(d in 1:nDistricts){
  district_name <- district_names[d]
      Rt %>% 
          filter(District %in% district_name) %>% 
          ggplot(aes(x = Date, y = smR)) +
          geom_hline(aes(yintercept = 1), colour = "#c9c9c9", linetype = "dashed") +
          geom_line(colour = "#0eede9") +
          geom_text(aes(x = dmy("12-04-2021"), y = 1), label = "Increasing spread", size = 3, vjust = -0.9, alpha = 0.3, colour = "#c9c9c9") +
          geom_text(aes(x = dmy("12-04-2021"), y = 1), label = "Getting contained", size = 3, vjust = 1.7, alpha = 0.3, colour = "#c9c9c9") +
          scale_x_date(breaks = breaks_pretty(length(unique(Rt$Date))/7)) +
          ylim(0, upper_y) +
          dark_theme_linedraw() +
          labs(title = paste0("Effective Reproductive Number for ", district_name, " District in Kerala, India"),
                 subtitle = "7-day simple moving average from 15 March 2021 till 6 May 2021",
                 caption = "data source: covid19india.org | @NishanShaik",
               y = "Effective R (smoothed)")

  ggsave(paste0(d,".jpeg"), 
         device = "jpeg", 
         dpi = "retina", 
         width = 35, 
         height = 25, 
         units = "in", 
         scale = 0.4)
  
  print(d)
}




### Single Facet Plot ================================
Rt %>% 
  ggplot(aes(x = Date, y = smR)) +
  geom_hline(aes(yintercept = 1), colour = "#c9c9c9", linetype = "dashed") +
  geom_line(colour = "#0eede9") +
  scale_x_date(breaks = breaks_pretty(length(unique(Rt$Date))/7)) +
  facet_wrap(~District, strip.position = "top") +
  dark_theme_classic() +
  theme(strip.background = element_rect(fill = "#383838")) +
  theme(strip.text.x = element_text(colour = "white", size = 9)) +
  labs(title = paste0("Effective Reproductive Number for Districts in Kerala, India"),
       subtitle = "Trends in 7-day simple moving averages from 15 March 2021 till 6 May 2021",
       caption = "data source: covid19india.org | @NishanShaik",
       y = "Effective R (smoothed)")
  
ggsave("01.jpeg", device = "jpeg", dpi = "retina", width = 35, height = 25, units = "in", scale = 0.4)


### State as a whole ==========================
state <- district %>% 
  filter(State == "Kerala") %>% 
  group_by(Date) %>% 
  summarise(Confirmed = sum(Confirmed, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(Incidence = Confirmed - lag(Confirmed, order_by = Date)) %>% 
  filter(Date >= begin_date)

est_r0 <- est.R0.TD(epid = state$Incidence, gen_time, t = state$Date, begin = 1, end = est_length)

Date <- names(est_r0[["R"]])
R <- est_r0[["R"]]
names(R) <-  NULL

Rts <- as_tibble(data.frame(Date, R))

Rts <- Rts %>% 
  mutate(smR = rollmean(R, k = 7, align = "right", fill = NA),
         Date = ymd(Date)) %>% 
  filter(Date >= dmy("15-03-2021"))

Rts %>% 
  ggplot(aes(x = Date, y = smR)) +
  geom_hline(aes(yintercept = 1), colour = "#c9c9c9", linetype = "dashed") +
  geom_line(colour = "#0eede9") +
  geom_text(aes(x = dmy("12-04-2021"), y = 1), label = "Increasing spread", size = 3, vjust = -0.9, alpha = 0.3, colour = "#c9c9c9") +
  geom_text(aes(x = dmy("12-04-2021"), y = 1), label = "Getting contained", size = 3, vjust = 1.7, alpha = 0.3, colour = "#c9c9c9") +
  scale_x_date(breaks = breaks_pretty(length(unique(Rts$Date))/7)) +
  ylim(0, upper_y) +
  dark_theme_linedraw() +
  labs(title = paste0("Effective Reproductive Number for Kerala, India"),
       subtitle = "7-day simple moving average from 15 March 2021 till 6 May 2021",
       caption = "data source: covid19india.org | @NishanShaik",
       y = "Effective R (smoothed)")

ggsave("001.jpeg", device = "jpeg", dpi = "retina", width = 35, height = 25, units = "in", scale = 0.4)
