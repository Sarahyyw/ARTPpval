############Real application: COVID19 desease survilence##############
#data is Provided by Johns Hopkins Center for Systems Science and Engineering (CSSE)

UScovid=read.csv("Real Application/time_series_covid19_confirmed_US.csv",header = F)
Day=as.Date(as.character(UScovid[1,12:815]),format = "%m/%d/%y")
UScovid=UScovid[-1,-c(1,2,3,4,5,8)] #3342 counties in total and 59 states, 804 days.
colnames(UScovid)=c("county","state","Lat","Long","location",as.character(Day))
state=unique(UScovid$state)
table(UScovid$state)
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

#location=list(Midwest=Midwest, Northeast=Northeast, 
              #Southeast=Southeast, Southwest=Southwest, West=West)

#shortname_West=c("AK","CA","CO","HI","ID","MT","NV","OR","UT","WA","WY")
#Other=c("Puerto Rico","Northern Mariana Islands","Diamond Princess","Guam","Grand Princess","American Samoa","Virgin Islands")
#shortname_Other=rep("NA",length(Other))
UScovid_input = data.frame(region=rep("NA",dim(UScovid)[1]),UScovid)
UScovid_input$region[which(UScovid$state%in% Northeast)]="Northeast"
UScovid_input$region[which(UScovid$state%in% Midwest)]="Midwest"
UScovid_input$region[which(UScovid$state%in% Southeast)]="Southeast"
UScovid_input$region[which(UScovid$state%in% Southwest)]="Southwest"
UScovid_input$region[which(UScovid$state%in% West)]="West"
#UScovid_input$region[which(UScovid$state%in% Other)]="Other"
#saveRDS(UScovid_input,file = "UScovid_original_Geog_region.rds")

table(UScovid_input$region); #table(UScovid_input$region[UScovid_input$region!="NA"])
#UScovid_input[which(UScovid_input$region=="NA"),1:6]
location=c("Midwest","Northeast","Southeast","Southwest","West")


m=7
baseline_len=14
pval_by_location=list()

for(i in 1:length(location)){
  input_data=UScovid_input[UScovid_input$region==location[i],]
  input_data=input_data[,-(1:(6))]
  baseline=input_data[,(1:baseline_len)]
  case=input_data[,21:ncol(input_data)]
  baseline=log(apply(baseline,c(1,2),as.numeric)+1)
  case=log(apply(case,c(1,2),as.numeric)+1)
  
  #len_case=floor(dim(case)[2]/m)
  pval_by_timeWindow=list()
  print(paste0("i=",i, " location=",location[i]))
  
  for(j in 1:5){
    
    case_test=case[,(((j-1)*m+1):((j)*m))]
    #print(sum(case_test))
    pval_vec = sapply(1:dim(case)[1],function(ind){
      set.seed(12345)
      if(mean(as.numeric(baseline[ind,]))==0 & mean(as.numeric(case_test[ind,]))==0){
        return(1)
      }else{
        case_count=as.numeric(case_test[ind,]) + rnorm(length(case_test[ind,]), 0, 1e-10)
        base_count=as.numeric(baseline[ind,]) + rnorm(length(baseline[ind,]), 0, 1e-10)
        tryCatch((t.test(case_count, base_count, alternative = "greater")$p.value), error = function(e) {
          return(1)
        })}
    })
    
    print(mean(pval_vec<=0.05))
    print(sum(pval_vec>1))
    pval_by_timeWindow[[j]]=pval_vec
  }
  pval_by_location[[i]]=pval_by_timeWindow
}

names(pval_by_location)=location
#saveRDS(pval_by_location,file = "P_individual_Geog_region.rds")



