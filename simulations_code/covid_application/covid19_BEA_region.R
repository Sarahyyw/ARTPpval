############Real application: COVID19 desease survilence##############
#data is Provided by Johns Hopkins Center for Systems Science and Engineering (CSSE)
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

m=7
baseline_len=14
pval_by_location=list()

for(i in 1:length(location)){
  input_data=UScovid_input[UScovid_input$BEA==location[i],]
  input_data=input_data[,-(1:6)]
  baseline=input_data[,(1:baseline_len)] #first 14 days as baseline
  case=input_data[,21:ncol(input_data)]
  
  baseline=log(apply(baseline,c(1,2),as.numeric)+1)
  case=log(apply(case,c(1,2),as.numeric)+1)
  pval_by_timeWindow=list()
  print(paste0("i=",i, " location=",location[i]))
  
  ## first five overlapping windows
  for(j in 1:5) { 
    case_test=case[,((j-1)*m+1):((j)*m)]
    
    pval_vec = sapply(1:dim(case)[1],function(ind){
      if(mean(as.numeric(baseline[ind,]))==0 & mean(as.numeric(case_test[ind,]))==0){
        return(1)
      }else{
        case_count=as.numeric(case_test[ind,])+rnorm(length(case_test[ind,]), 0, 1e-10)
        base_count=as.numeric(baseline[ind,])+rnorm(length(baseline[ind,]), 0, 1e-10)
        tryCatch((t.test(case_count, base_count, alternative = "greater")$p.value), error = function(e) {
          return(1)
        })}
    })
    print(mean(pval_vec<=0.05))
    print(sum(pval_vec>1))
    pval_by_timeWindow[[j]]=pval_vec
    #print(j)
  }
  pval_by_location[[i]]=pval_by_timeWindow
  #print(paste0("i=",i, " location=",location[i]))
}

names(pval_by_location)=location
#saveRDS(pval_by_location,file = "P_individual_BEA_region.rds")


