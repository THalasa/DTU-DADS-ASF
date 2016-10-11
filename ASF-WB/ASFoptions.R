#### Creating ASF list of options for ASF() version 0.15.1
#### Includes the following functions:
####
#### ASFoptions
####
#### The output is a list of parameters.
#### The output of the model is produced by the summary function.

ASFoptions<-function(...){
  args<-list(...)

  ## List of default values:
  defaults<-list(
                 ############################################################
                 n=1,                   # Number of simulated epidemics
                 runID=NULL,            # ID used for output and temp. files
                                        #  using a random number if NULL
                 stepInFile=NULL,       # File used to step into an outbreak
                                        #  and start all iterations from
                                        #  same setup.
                 maxTime=365,           # Maximum length of each outbreak
                 ##DF=3,
                 Tstoch=FALSE,          # Number of infected animals in newly
                                        #  infected herd: stochastic (T) or
                                        #  exactly 1 (F)
                 RFstoch=FALSE,         # Number of intraherd disease
                                        #  transmissions:
                                        #  binomial chain model (T) or
                                        #  Reed-Frost model (F)
                 seed=NULL,             # Seed for random number generator.
                                        #  If negative the seed is found as:
                                        #  set.seed(iteration+abs(seed))
                 pause=0,               # Time in seconds to pause between
                                        #  new graphs.
                 delaySteps=1,          # Vector of starting times for time
                                        #  dependent diagnosis delays
                 delayTimes=expression(round(rpert(n,1,2,4))), # Vector of
                                        #  diagnosis delays for starting times
                                        #  in delaySteps. These should be
                                        #  expressions or text that can be
                                        #  parsed to expressions
                 ignoreStatus=TRUE,     # T/F value for ingoring disease
                                        #  status information in input file
                                        #  (all reset to 1 if True)
                 indexHerdFunction="selectIndexHerd", # Name of function used
                                        #  to select index herd for each
                                        #  simulated epidemic (iteration)
                 indexHerdSelect=list(herdType=1:18), # Argument to function 
                 #basicScenario=TRUE,
                                        #  used to select index herd for each
                                        #  simulated epidemic (iteration)
                 indexDirect=FALSE,     # Index herd infected by direct (T)
                                        #  or indirect (F) contact
                 depopTeams=Inf,        # Number of locations that can be
                                        #  depoped at a time. not used anymore
                 Capacity=c(4800),      # The culling capacity per day
                                        #   4800 swine
                 TracePeriod=30,        # Tracing will go back to the defined number of days; default 30 days   
                 #CapSurvay=c(450),     # Surveillance capacity; number of herds/day
                 
                 RepVisSurvZone=14,     # How often the visit within the surveillance zone be repeated
                 traceBack=FALSE,       # Should traceback information be saved
                 traceBackDelay=1,
                 rateBefD=1,            # Abattoir rate before detection first 
                                        # infected herd for swine when there is
                                        # a movement from swine herds to slaughter
                 rateAftD=2,            # Abattoir rate after detection first 
                                        # infected herd when there is
                                        # a movement from the herds to slaughter (less frequent visits despite of a higher value for rate
                                        # because the Exp function will compress it more ;-))
                 FirstDetPara=0.0255,   # proportion of sick and dead animals
                 ProbSelDiag=1,         # The probability of diagnosing a selected herd for diagnosis
                 Detailed=FALSE,        # should detailed surveillance output be printed
                 DumpData=1,            # Number of data lines that can be reached before
                                        # data about survyed herds can be dumpped in the output file
                 ProbSelPV1=0,          # proportion of herds that will be tested (PCR) during first protection zone visit
                 ProbSelSV1=0,          # proportion of herds that will be tested (PCR) during first surveillance zone visit
                 ProbSelSV2=0,          # proportion of herds that will be tested during second surveillance zone visit
                 ProbSelTIDC=0.1,       # proportion of traced herds from indirect contacts that will be tested (PCR) visit
                 SecSurVisitOLSZ=0,     # number of days. herds in overlapping surveillance zones will get a new visit every SecSurVisitOLSZ days.
                 DelayStartVisitOLPZ=0, # number of days. herds in overlapping protection zone will get a new visit every DelayStartVisitOLPZ once the
                                        # they continue in the protection zone and the time of second PZ visit has passed 
                 SecSurVisit=40,        # surveillance visit. called second here because the first is not mandatory.
                 SecProtVisit=45,       # second visit in protection zone
                 firstSurvVisit=FALSE,  # Allow first surveillance visit(Yes/No)
                 DelayStartVisit=2,     # Number of days before the visiting of herds for surveillance would start
                 DelayVisit=7,          # delay for the extra visits for herds in overlapping zones.
                 MortalityIncrease=2,   # level of increase in mortality before potential detection.
                 MortalityIncreaseZone=1.5,   # level of increase in mortality before potential detection.
                 InfPropSub=0.1,        # risk of infection from subclinical animals (before clinical signs appeared)
                 PerDeadAnim=0.95,      # percentage of animals that die following infection
                 DaysDead=5,            # number of past days to be used to determine infectiousness of leftovers of dead animals
                 #ReqSampSiz=30,         # sample size when herds are tested.
                 DeadImpact=1,          # parameter to address the impact of leftovers on disease spread
                 ImpDeadTime=1,         # parameter to address uncertainty of survivability of virus in leftovers
                 PZoneDuration=50,      # protection Zones duration should be 50 days at the start of each iteration
                 SZoneDuration=45,      # Surveillance Zones duration should be 45 days at the start of each iteration
              SerologyTesting=c(1,5,6), # The type of visit (1=PV2, 2=PV1, 3=SV1, 4=SV2, 5=trace IDC, 6=Trace DC) where serology testing will be applied.
                 PCRTesting=c(5,6),     # The type of visit (1=PV2, 2=PV1, 3=SV1, 4=SV2, 5=trace IDC, 6=Trace DC) where PCR testing will be applied.
                 NumDeadAnimFirst=5,    # Number of dead animals in the herd for first detection
                 NumDeadAnimAftFirst=1, # Number of dead animals in the herds for detection after first detection occured
                 NumDeadAnimSurv=1,     # Number of dead animals in the herd for detection through surveillance.
                 DaysSurDead=7,         # Number of past days to be used to survay dead animals                                
                 numTestDead=5,         # Number of dead animals tested 
                 DelaySubDeadSamp=1,    # Delay on the submited samples to arrive to the laboratory

ToTracedIDC=2:4,

probSelectTIDC=c(0,0.838,0.2,0.125), # the probability that a movement will not be forgetten and it will be traced and visited,
LocSpLim=c(0.1,0.5,1,2) # cutoffs in km for local spread distance probabilities
,
DistList=list(                               ## a list that includes the probability of infection (distProb) through local spread given the distance in km (distCat)
  distCat=c(1,2,3,4,5),                      ## from the infectious herd
  distProb=c(0.1,0.006,0.002,0.000015,0)     ## these probabilities are based on Boklund et al. (2009) for CSF after reduction by 50% to include
  )                                          ## the lower infectivity of ASF as carried out by Nigsch et al., 2013.
  , 
probList=list(   ## Default distributions
  DistCat=c(0,1,3,10,15,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,300), ## distance categories in km
  distcat=c(0,1,3,10,15,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,300)

  ),
LambdaWB=c(0.10824,0.055,0)# the categories represent distances of highRisk, lowRisk and noRisk herds
                 ,
LocalSPWB=c(0.0006,0.0002,0.0000015,0),# the categories represent distances of highRisk (within 0.5km), medium Risk (from 0.5 to 1km) 
                                       # low risk (from 1 to 2km) and no Risk herds (>2km)
                                       # these values are from Boklund et al. 2008 and reduced 20 fold, due to 
                                       # lower spread of ASF than CSF (Nigsch et al., 2013) and smaller herd sizes
                                       # of wild boar
                 
                 newInfFunctions=c(     # Vector of functions used to make new infections (including parameters).
                   "DIRinf3('LamAll',MovSwProb,'pMatAll','RiskDC',MovMatAll,restMovedSize=35,label=1)", 
                   "DIRinf3('LambdaWeaners',MovWeProb,'pMatWea','RiskDC',MovMatWean,restMovedSize=10,label=1)", 
                   "INDflex('LamAb',SwMovAbProb,'relDC','pMatMovAb','RiskAb',probMatrix=MovAb,Reduction=0.5,Abattoir=TRUE,label=2)",
                   "INDflex('LamMRC',MedRiskMovProb,'relIMC','pMatMRC','RiskMRC',Reduction=1,label=3)",
                   "INDflex('LamLRC',LowRiskMovProb,'relILC','pMatLRC','RiskLRC',Reduction=1,label=4)",
                   "LASinf(localsize=2,label=5)",
                   "WildBoar(relCont='relWB',RiskVar='RiskCatToWB2',ProbCont='RiskWB',label=6)",
                   "LSWildBoar(RiskVar='LocSpWB',label=7)"),

                 controlFunctions=c(    # Vector of functions used for movement controls, tracing and surveillance
                   "controlAll(effectDC='rpert(n,0.95,0.98,1)',label='SS')",
                   "controlDiag(effectDC=1,effectIMC='rpert(n,0.7,0.8,0.95)',effectILC='rpert(n,0.95,0.98,1)',label='CD')",
                   "SurvZone(size=10,effectDC='rpert(n,0.95,0.98,1)',effectIMC='rpert(n,0.7,0.8,0.95)',effectILC='rpert(n,0.2,0.3,0.5)',label='SZ')",
                   "ProtZone(size=3,effectDC='rpert(n,0.95,0.98,1)',effectIMC='rpert(n,0.7,0.8,0.95)',effectILC='rpert(n,0.2,0.3,0.5)',label='PZ')",
                   "traceDC(prob=0.99,probdetect=0.95,delay='round(rpert(n,1,2,3))',tracetime='round(runif(n,0,2))',duration=30,label='traceDirect')",
                   "traceIDC(timetotrace='round(runif(n,0,4))',delayvisitMed='round(rpert(n,0,1,2))',delayvisitLow='round(rpert(n,0,2,4))',duration=30,label='traceInDirect')",
                   "SurvZonesHerds()"),
                  # "SurvDead()"),
                   
                                  
                 ############################################################
                 ## Files
                 infofile="DataDADSASFWB.csv",    # File with herd locations and type 
                 typesfile="typesfile.csv",  # Definitions of type parameters
                 runfile="",                 # File used for additional output
                fileMovMatAll="MovMatAll.csv", # bla bla
                fileMovMatWean="MovMatWean.csv",
                fileMovAb="MovAb.csv",
                fileMovSwProb="MovSwProb.csv",
                fileMovWeProb="MovWeProb.csv",
                fileSwMovAbProb="SwMovAbProb.csv",
                fileMedRiskMovProb="MedRiskMovProb.csv",
                fileLowRiskMovProb="LowRiskMovProb.csv",

                chroniclefile=FALSE,   # File name or FALSE to use runID.

                 ############################################################
                 ## Output, Graphs and text
                 hideMap=NULL,          # Hide map (T) or show map (F)
                                        #  while running
                 itupdate=10,           # Update period (number of iterations)
                                        #  for summary graphs
                 tornCol=NULL,          # Typesfile column to reduce for
                                        #  tornado plot (default NULL)
                 tornMult=0.9,          # Multiplier for tornCol column of
                                        #  typesfile
                 hidePlots=FALSE,       # Hide (T) or show (F) summary plots,
                                        #  including final risk map

                 summaryFunction="sumTh", # Name of function used for
                                        #  summaries. It is called with
                                        #  arguments: "init", "day", "iter",
                                        #  and "final"
                 verbose=FALSE,         # Make verbose output while running

                 ############################################################
                 ## movement control
                 ## this is done in the initialization function
                 #gDaysUntilBaseline=eval(parse('round(rpert(n,18,21,23))')), # Number of days until baseline
                                        #  controls go into effect or NA to
                                        #  get random detection using a
                                        #  distribution in the types file.
                 
                 ############################################################
                 interventionFunctions=c( # Vector of functions used to make 
                                          #  interventions (including parameters)
                   "DummyInter()"
                   ),                   #  Format: strings or expressions.
                 #DayStartcull=14,       # default day to start pre-emptive culling
                 cullTypes=c(1:18),     # herd types to be culled
                 CullDays=1:365,        # which days culling should be considered
                 ###########################################################
                 ## Economic Values. Costs in DKK
                 CostAnimPCR= 530,        # costs of PCR testing per animal
                 CostAnimSer= 68,         # costs of serology testing per animal
                 CosIndsendelse=179,      # costs of sending a package of max 1kg. from postDanmark 
                 CostsSamTesSer=137,      # costs for sampling (traveling time, sampling time, sending the samples and materials) and testing (serology) per sample from a herd. time used 3 hours
                 CostsSamTesPCR=599,      # costs for sampling (traveling time, sampling time, sending the samples and materials) and testing (serology) per sample from a herd. time used 3 hours

                 #Vet/hour = 800kr., technician/hour=400kr., 3 hours for 60 animals and 2 hours for 30 animals (including traveling), materials 9kr/sample
                 #CostsSamTesSer60=8220,   # costs for sampling (traveling time, sampling time, sending the samples and materials) and testing (serology) 60 animals from a herd. time used 3 hours
                 #CostsSamTesPCR60=35940,  # costs for sampling (traveling time, sampling time, sending the samples and materials) and testing (serology) 60 animals from a herd. time used 3 hours
                 #CostsSamTesSer30=4710,   # costs for sampling (traveling time, sampling time, sending the samples and materials) and testing (serology) 30 animals from a herd. time used 2 hours 
                 #CostsSamTesPCR30=18570,  # costs for sampling (traveling time, sampling time, sending the samples and materials) and testing (serology) 30 animals from a herd. time used 2 hours
                 CostsVisFarSw=1565,      # costs of clinical surveillence swine farms
                 logCosSw=92956,          # logistic costs per herd for swine herds
                 clCosSw=935.77,          # cleaning and disinfection costs per animal for swine herds
                 clCosFin=230.5,          # cleaning and disinfection costs per animal for swine herds
                 clCosPig=1000000,        # cleaning and disinfection costs per herd swine herd
                 CESPSow=9.99,            # costs of empty sow stable per day
                 CESPFin=0.55,            # costs of empty finisher stable per day
                 CosWSFinishers=750,      # costs of welfare slaughter of finishers
                 FactWSFinishers=0.011,   # factor for finishers welfare slaughter per sow
                 CosWSWeaners=375,        # costs of welfare salughter of weaners
                 FactWSWeaners=0.067,     # factor for weaners welfare slaughter per sow
                 comSows=3365.34,         # compensation per sow
                 comFin=521.91,           # compensation per finisher pig
                 WeanValue=4.2,           # value of wanears for compensation
                 expBanAftCEU=0,          # export ban after culling the last infected herd and lifting the zone in relation to EU countries
                 expBanAftCNEU=90,        # export ban after culling the last infected herd and lifting the zone to non EU countries
                 totExpSwNEU=28286093,    # total returns on export of live swine and pork products to non-EU countries per day
                 totExpLivSwEU=15662784,  # total returns on export of live swine to EU countries per day
                 totExpSwProdEU=55500143, # total returns on export of pork products to EU countries per day
                 reducedPrice=0.25        # reduction in the price for products to non-EU
                  )

  if (length(args)>0){
    ## Check which values to change
    changes<-match(names(args),names(defaults))
    ## Stop if trying to change nonexisting parameter
    if (any(is.na(changes)))
      stop(paste( paste(names(args)[is.na(changes)],collapse=", ") ,
                 " is/are not valid parameter names. \n See help file."))
    ## Change the requested values and return
    for (i in 1:length(changes)){
      defaults[[changes[i]]]<-args[[i]]
    }
  }##EndOf if (length(args))
  ## Converting strings to functions
  defaults$indexHerdFunction<-match.fun(defaults$indexHerdFunction)
  defaults$summaryFunction<-match.fun(defaults$summaryFunction)
  if (defaults$verbose) cat("Leaving ASFoptions. ")
  return(defaults)
}##EndOf ASFoptions
