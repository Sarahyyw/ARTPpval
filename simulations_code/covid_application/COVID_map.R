# After getting the individual p-value for the 5 time period, here we want to first plot the time map for all region and also the combined p-value plots

################ Geog_Region #########################
UScovid=read.csv("/Users/sarah/Desktop/ARTP/Real Application/time_series_covid19_confirmed_US.csv",header = F)
Day=as.Date(as.character(UScovid[1,12:815]),format = "%m/%d/%y")
UScovid=UScovid[-1,-c(1,2,3,4,5,8)] #3342 counties in total and 59 states, 804 days.
colnames(UScovid)=c("county","state","Lat","Long","location",as.character(Day))
state=unique(UScovid$state)

Northeast=c("Connecticut", 
            "Maine", 
            "Massachusetts", 
            "New Hampshire", 
            "New Jersey", 
            "New York", 
            "Pennsylvania", 
            "Rhode Island", 
            "Vermont",
            "District of Columbia",
            "Delaware",
            "Maryland")

#location_Northeast=rep("Northeast",length(Northeast))
#shortname_Northeast=c("CT","ME","MA","NH","NJ","NY","PA","RI","VT","DC","DE","MD")

Midwest=c(
  "Illinois", 
  "Indiana",
  "Iowa", 
  "Kansas", 
  "Michigan", 
  "Minnesota", 
  "Missouri", 
  "Nebraska", 
  "North Dakota", 
  "Ohio", 
  "South Dakota", 
  "Wisconsin" 
)
#shortname_Midwest=c("IL","IN","IA","KS","MI","MN","MO","NE","ND"<"OH","SD","WI")

Southeast=c(
  "Alabama",
  "Arkansas",
  "Florida",
  "Georgia", 
  "Kentucky", 
  "Louisiana", 
  "Mississippi", 
  "North Carolina", 
  "South Carolina", 
  "Tennessee", 
  "Virginia",
  "West Virginia"
)
#shortname_Southeast=c("AL","AR","FL","GA","KY","LA","MS","NC","SC","TN","VA","WV")

Southwest=c(
  "Arizona", 
  "New Mexico", 
  "Oklahoma", 
  "Texas" 
)
#shortname_Southwest=c("AZ","NM","OK","TX")

West=c(
  "Alaska",
  "California", 
  "Colorado",
  "Hawaii",
  "Idaho",
  "Montana",
  "Nevada",
  "Oregon",
  "Utah",
  "Washington",
  "Wyoming"
)


UScovid_input = data.frame(region=rep("NA",dim(UScovid)[1]),UScovid)
UScovid_input$region[which(UScovid$state%in% Northeast)]="Northeast"
UScovid_input$region[which(UScovid$state%in% Midwest)]="Midwest"
UScovid_input$region[which(UScovid$state%in% Southeast)]="Southeast"
UScovid_input$region[which(UScovid$state%in% Southwest)]="Southwest"
UScovid_input$region[which(UScovid$state%in% West)]="West"

pval_by_location = readRDS("/Users/sarah/Desktop/ARTP/Real Application/P_individual_Geog_region.rds")
regions_list <- list(
  Midwest = Midwest,
  Northeast = Northeast,
  Southeast = Southeast,
  Southwest = Southwest,
  West = West
)

# Prepare county data
us_states <- map_data("state")
us_counties <- map_data("county")
colnames(us_counties) <- c("long","lat","group","order","state","county")

# Loop over regions
Region_df_list <- list()

for (i in seq_along(regions_list)) {
  
  Region_name <- names(regions_list)[i]
  Region_states <- regions_list[[i]]
  
  pval_by_timeWindow <- pval_by_location[[i]]  # matching order to regions_list
  UScovid_region <- UScovid_input %>% filter(region == Region_name)
  
  countiesP <- cbind(
    UScovid_region[, 2:6],
    T1 = pval_by_timeWindow[[1]],
    T2 = pval_by_timeWindow[[2]],
    T3 = pval_by_timeWindow[[3]],
    T4 = pval_by_timeWindow[[4]],
    T5 = pval_by_timeWindow[[5]]
  )
  
  countiesP[is.na(countiesP)] <- 1
  countiesP$state <- tolower(countiesP$state)
  countiesP$county <- tolower(countiesP$county)
  #us_counties$county[is.na(match(us_counties$county, countiesP$county))] <- "unassigned"
  
  Region_states <- tolower(Region_states)
  us_countiesP <- left_join(us_counties, countiesP)
  
  cols <- paste0("T", 1:5)
  for (col in cols) {
    us_countiesP[[col]][us_countiesP[[col]] <= 10^-3] <- 10^-3
  }
  
  Region_df <- us_countiesP %>% filter(state %in% Region_states)
  Region_df_list[[Region_name]] <- Region_df
}


Region_long <- Region_df_list$Southeast %>%
  pivot_longer(cols = starts_with("T"),
               names_to = "Test",
               values_to = "pval") %>%
  mutate(logp = -log10(as.numeric(pval)))

ggplot(Region_long, aes(x = long, y = lat, group = group, fill = logp)) +
  geom_polygon(color = "gray90", size = 0.05) +
  coord_equal() +
  scale_fill_gradient(
    low = "white", high = "#99000D", na.value = "white"
    #limits = c(0, 3)
  ) +
  labs(fill = "-log10 p-value") +
  facet_wrap(~ Test, nrow = 1) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "gray90", color = NA),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "bottom",
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(1.0, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )


Region_long_all <- bind_rows(
  lapply(names(Region_df_list), function(region_name) {
    Region_df_list[[region_name]] %>%
      mutate(Region = region_name) %>%
      pivot_longer(cols = starts_with("T"),
                   names_to = "Test",
                   values_to = "pval") %>%
      mutate(logp = -log10(as.numeric(pval)))
  })
)

ggplot(Region_long_all, aes(x = long, y = lat, group = group, fill = logp)) +
  geom_polygon(color = "gray90", size = 0.05) +
  scale_fill_gradient(low = "white", high = "#99000D",na.value = "white") +
  labs(fill = "-log10 p-value") +
  facet_grid(rows = vars(Region), cols = vars(Test), scales = "free") +  # allow zoom per region
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "gray90", color = NA),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(color = "black", face = "bold"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(1.0, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )



################ EBA region #################
UScovid=read.csv("~/Desktop/ARTP/Real Application/time_series_covid19_confirmed_US.csv",header = F)
Day=as.Date(as.character(UScovid[1,12:815]),format = "%m/%d/%y")
UScovid=UScovid[-1,-c(1,2,3,4,5,8)] #3342 counties in total and 59 states, 804 days.
colnames(UScovid)=c("county","state","Lat","Long","location",as.character(Day))

state=unique(UScovid$state)
# numbers=table(UScovid$state)
# numbers=numbers[which(numbers>=80)]

NewEngland=c(
  "Connecticut",
  "Maine",
  "Massachusetts",
  "New Hampshire",
  "Rhode Island",
  "Vermont"
)

GreatLakes=c(
  "Wisconsin", 
  "Michigan", 
  "Illinois", 
  "Indiana", 
  "Ohio"
)

Plains=c("Iowa",
         "Kansas",
         "Minnesota",
         "Missouri",
         "Nebraska",
         "North Dakota",
         "South Dakota"
)

RockyMountain=c(
  "Colorado",
  "Idaho",
  "Montana",
  "Utah",
  "Wyoming"
)

Southeast=c(
  "Alabama",
  "Arkansas",
  "Florida",
  "Georgia", 
  "Kentucky", 
  "Louisiana", 
  "Mississippi", 
  "North Carolina", 
  "South Carolina", 
  "Tennessee", 
  "Virginia",
  "West Virginia"
)

Southwest=c(
  "Arizona", 
  "New Mexico", 
  "Oklahoma", 
  "Texas" 
)

Mid_Atlantic =c("Delaware","District of Columbia","Maryland","New Jersey","New York","Pennsylvania")
Farwest=c("Alaska","California","Hawaii","Nevada","Oregon","Washington")

UScovid_input = data.frame(BEA=rep("NA",dim(UScovid)[1]), UScovid)
UScovid_input$BEA[which(UScovid$state%in% NewEngland)]="NewEngland"
UScovid_input$BEA[which(UScovid$state%in% Mid_Atlantic)]="Mid_Atlantic"
UScovid_input$BEA[which(UScovid$state%in% Plains)]="Plains"
UScovid_input$BEA[which(UScovid$state%in% GreatLakes)]="GreatLakes"
UScovid_input$BEA[which(UScovid$state%in% Southeast)]="Southeast"
UScovid_input$BEA[which(UScovid$state%in% Southwest)]="Southwest"
UScovid_input$BEA[which(UScovid$state%in% RockyMountain)]="RockyMountain"
UScovid_input$BEA[which(UScovid$state%in% Farwest)]="Farwest"

table(UScovid_input$BEA)
location=list(NewEngland=NewEngland, Mid_Atlantic=Mid_Atlantic,Plains=Plains,
              GreatLakes=GreatLakes,RockyMountain=RockyMountain,Farwest=Farwest)
location=names(location)


pval_by_location = readRDS("/Users/sarah/Desktop/ARTP/Real Application/P_individual_BEA_region.rds")
regions_list <- list(
  NewEngland = NewEngland,
  Mid_Atlantic = Mid_Atlantic,
  Plains = Plains,
  GreatLakes = GreatLakes,
  RockyMountain = RockyMountain,
  Farwest = Farwest
)

# Prepare county data
us_states <- map_data("state")
us_counties <- map_data("county")
colnames(us_counties) <- c("long","lat","group","order","state","county")

# Loop over regions
Region_df_list <- list()

for (i in seq_along(regions_list)) {
  
  Region_name <- names(regions_list)[i]
  Region_states <- regions_list[[i]]
  
  pval_by_timeWindow <- pval_by_location[[i]]  # matching order to regions_list
  UScovid_region <- UScovid_input %>% filter(BEA == Region_name)
  
  countiesP <- cbind(
    UScovid_region[, 2:6],
    T1 = pval_by_timeWindow[[1]],
    T2 = pval_by_timeWindow[[2]],
    T3 = pval_by_timeWindow[[3]],
    T4 = pval_by_timeWindow[[4]],
    T5 = pval_by_timeWindow[[5]]
  )
  
  countiesP[is.na(countiesP)] <- 1
  countiesP$state <- tolower(countiesP$state)
  countiesP$county <- tolower(countiesP$county)
  #us_counties$county[is.na(match(us_counties$county, countiesP$county))] <- "unassigned"
  
  Region_states <- tolower(Region_states)
  us_countiesP <- left_join(us_counties, countiesP)
  
  cols <- paste0("T", 1:5)
  for (col in cols) {
    us_countiesP[[col]][us_countiesP[[col]] <= 10^-3] <- 10^-3
  }
  
  Region_df <- us_countiesP %>% filter(state %in% Region_states)
  Region_df_list[[Region_name]] <- Region_df
}


Region_long_all <- bind_rows(
  lapply(names(Region_df_list), function(region_name) {
    Region_df_list[[region_name]] %>%
      mutate(Region = region_name) %>%
      pivot_longer(cols = starts_with("T"),
                   names_to = "Test",
                   values_to = "pval") %>%
      mutate(logp = -log10(as.numeric(pval)))
  })
)

ggplot(Region_long_all, aes(x = long, y = lat, group = group, fill = logp)) +
  geom_polygon(color = "gray90", size = 0.05) +
  scale_fill_gradient(low = "white", high = "#99000D",na.value = "white") +
  labs(fill = "-log10 p-value") +
  facet_grid(rows = vars(Region), cols = vars(Test), scales = "free") +  # allow zoom per region
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "gray90", color = NA),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(color = "black", face = "bold"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(1.0, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

Region_long <- Region_df_list$NewEngland %>%
  pivot_longer(cols = starts_with("T"),
               names_to = "Test",
               values_to = "pval") %>%
  mutate(logp = -log10(as.numeric(pval)))

ggplot(Region_long, aes(x = long, y = lat, group = group, fill = logp)) +
  geom_polygon(color = "gray90", size = 0.05) +
  coord_equal() +
  scale_fill_gradient(
    low = "white", high = "#99000D", na.value = "white"
    #limits = c(0, 3)
  ) +
  labs(fill = "-log10 p-value") +
  facet_wrap(~ Test, nrow = 1) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "gray90", color = NA),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "bottom",
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(1.0, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )



############## state region ###########
UScovid=read.csv("~/Desktop/ARTP/Real Application/time_series_covid19_confirmed_US.csv",header = F)
Day=as.Date(as.character(UScovid[1,12:815]),format = "%m/%d/%y")
UScovid=UScovid[-1,-c(1,2,3,4,5,8)] #3342 counties in total and 59 states, 804 days.
colnames(UScovid)=c("county","state","Lat","Long","location",as.character(Day))

state=unique(UScovid$state)
numbers=table(UScovid$state)
numbers=numbers[which(numbers>=80)]
location=names(numbers)
UScovid_input = UScovid
table(UScovid_input$state)

pval_by_location = readRDS("/Users/sarah/Desktop/ARTP/Real Application/P_individual_state_region.rds")

# Prepare county data
us_states <- map_data("state")
us_counties <- map_data("county")
colnames(us_counties) <- c("long","lat","group","order","state","county")

# Loop over regions
regions_list = names(pval_by_location); names(regions_list) = names(pval_by_location)
Region_df_list <- list()

for (i in seq_along(regions_list)) {
  
  Region_name <- names(regions_list)[i]
  Region_states <- regions_list[[i]]
  
  pval_by_timeWindow <- pval_by_location[[i]]  # matching order to regions_list
  UScovid_region <- UScovid_input %>% filter(state == Region_name)
  
  countiesP <- cbind(
    UScovid_region[, 1:6],
    T1 = pval_by_timeWindow[[1]],
    T2 = pval_by_timeWindow[[2]],
    T3 = pval_by_timeWindow[[3]],
    T4 = pval_by_timeWindow[[4]],
    T5 = pval_by_timeWindow[[5]]
  )
  
  countiesP[is.na(countiesP)] <- 1
  countiesP$state <- tolower(countiesP$state)
  countiesP$county <- tolower(countiesP$county)
  #us_counties$county[is.na(match(us_counties$county, countiesP$county))] <- "unassigned"
  
  Region_states <- tolower(Region_states)
  us_countiesP <- left_join(us_counties, countiesP)
  
  cols <- paste0("T", 1:5)
  for (col in cols) {
    us_countiesP[[col]][us_countiesP[[col]] <= 10^-3] <- 10^-3
  }
  
  Region_df <- us_countiesP %>% filter(state %in% Region_states)
  Region_df_list[[Region_name]] <- Region_df
}

Region_long <- Region_df_list$Georgia %>%
  pivot_longer(cols = starts_with("T"),
               names_to = "Test",
               values_to = "pval") %>%
  mutate(logp = -log10(as.numeric(pval)))

ggplot(Region_long, aes(x = long, y = lat, group = group, fill = logp)) +
  geom_polygon(color = "gray90", size = 0.05) +
  coord_equal() +
  scale_fill_gradient(
    low = "white", high = "#99000D", na.value = "white"
    #limits = c(0, 3)
  ) +
  labs(fill = "-log_10 P-value") +
  facet_wrap(~ Test, nrow = 1) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "gray90", color = NA),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "bottom",
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(1.0, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  ggtitle("Georgia")






