# A function that will calculate our decision table for TLM or SEAM, I hope
#mod.select:      The model you are running, "TLM" and "SEAM" are current options
#data:            The model output
#catch.scenarios  The catch removals you want to build the table for, a vector of numbers is what it needs
#n.sims           The number of simulations to run, default is 1e6 which runs in seconds
#TRP              The target reference point, if it exists
#USR              The Upper Stock Reference point, if it exists
#LRP              The Limit(lower) Reference Point, if it exists
#RR               The Removal Reference Point, if it exists. In % so use 10 not 0.1!
dec.tab <- function(mod.select = "SEAM",data = NULL, catch.scenarios = seq(0,10000,by=50),n.sims = 1e6,TRP = NULL,USR = NULL,LRP = NULL, RR = NULL)
{
library(SEBDAM)

# Number of catch scenarios
n.catch.scenarios <- length(catch.scenarios)


pred.proc <- get_processes(data)

# Growth data.
gs <- data$obj$env$data$gI[length(data$obj$env$data$gI)]
gRs <- data$obj$env$data$gR[length(data$obj$env$data$gR)]

if(mod.select == "TLM")
{
  B.log <- pred.proc$log_processes$log_B[length(pred.proc$log_processes$log_B)-1]
  B.log.se <- pred.proc$log_processes$se_log_B[length(pred.proc$log_processes$se_log_B)-1]
  R.log <- pred.proc$log_processes$log_R[length(pred.proc$log_processes$log_R)-1]
  R.log.se <- pred.proc$log_processes$se_log_R[length(pred.proc$log_processes$se_log_R)-1]
  m.log <- pred.proc$log_processes$log_m[length(pred.proc$log_processes$log_m)-1]
  m.log.se <- pred.proc$log_processes$se_log_m[length(pred.proc$log_processes$se_log_m)-1]
}

if(mod.select != "TLM")
{
  B.log <- pred.proc$log_tot_frame$log_totB[length(pred.proc$log_tot_frame$log_totB)-1]
  B.log.se <- pred.proc$log_tot_frame$se_log_totB[length(pred.proc$log_tot_frame$se_log_totB)-1]
  R.log <- pred.proc$log_tot_frame$log_totR[length(pred.proc$log_tot_frame$log_totR)-1]
  R.log.se <- pred.proc$log_tot_frame$se_log_totR[length(pred.proc$log_tot_frame$se_log_totR)-1]
  m.log <- pred.proc$log_tot_frame$log_mean_m[length(pred.proc$log_tot_frame$log_mean_m)-1]
  m.log.se <- pred.proc$log_tot_frame$se_log_mean_m[length(pred.proc$log_tot_frame$se_log_mean_m)-1]
}


# Initialize some data...
decision.table <- data.frame(catch = rep(NA,n.catch.scenarios), exploit = rep(NA,n.catch.scenarios),
                             per.diff = rep(NA,n.catch.scenarios), B.change = rep(NA,n.catch.scenarios), 
                             prob.decline = rep(NA,n.catch.scenarios))

for(i in 1:n.catch.scenarios)
{
  # So now sample from the log normals using the above data
  Bs.tot <- rlnorm(n.sims,B.log,B.log.se)
  Rs.tot <- rlnorm(n.sims,R.log,R.log.se)
  ms     <- rlnorm(n.sims,m.log,m.log.se)
  # Run the model through the catch scenarios
  Bst <- (exp(-ms))*gs*(Bs.tot-catch.scenarios[i])
  Rst <- (exp(-ms))*gRs*(Rs.tot) 
  B2 <- Bst + Rst
  exploit <- 100* (catch.scenarios[i]/(B2+ catch.scenarios[i]))
  B.diff <- B2 - exp(B.log) 
  B.per.diff <- 100*((B2 - exp(B.log) )/exp(B.log) )
  decision.table$catch[i] <- as.numeric(catch.scenarios[i])
  decision.table$per.diff[i] <- signif(median(B.per.diff),digits=3)
  decision.table$B.change[i]   <- signif(median(B.diff),digits=3)
  decision.table$prob.decline[i] <- signif(length(B2[B2 < exp(B.log)]) / n.sims,digits = 2)
  # Making these so decimal places look nice
  if(median(exploit) == 0) decision.table$exploit[i]  <- "0.0"
  if(median(exploit) < 1 & median(exploit) <= 0.95 & median(exploit) != 0) decision.table$exploit[i]  <- as.character(signif(median(exploit),digits=1))
  if(median(exploit) < 1 & median(exploit) >=0.95) decision.table$exploit[i]  <- "1.0"
  if(median(exploit) >= 1) decision.table$exploit[i]  <- as.character(signif(median(exploit),digits=2))
  # If we have reference points add these to the table.
  if(!is.null(LRP)) 
  {
    raw.LRP <- length(B2[B2 < LRP]) / n.sims
    decision.table$prob.below.LRP[i] <- signif(raw.LRP,digits = 2)
    if(raw.LRP > 0.995) decision.table$prob.below.LRP[i] <- "> 0.99"
    if(raw.LRP < 0.00995) decision.table$prob.below.LRP[i] <- "< 0.01"
    if(raw.LRP >= 0.00995 & raw.LRP < 0.0995) decision.table$prob.below.LRP[i] <- signif(raw.LRP, digits = 1)
  }
  if(!is.null(USR)) 
  {
    raw.USR <- length(B2[B2 < USR])/ n.sims
    decision.table$prob.below.USR[i] <- signif(raw.USR,digits = 2)
    if(raw.USR > 0.995) decision.table$prob.below.USR[i] <- "> 0.99"
    if(raw.USR < 0.00995) decision.table$prob.below.USR[i] <- "< 0.01"
    if(raw.USR >= 0.00995 & raw.USR < 0.0995) decision.table$prob.below.USR[i] <- signif(raw.USR, digits = 1)
  }
  if(!is.null(TRP)) 
  {
    raw.TRP <- length(B2[B2 < TRP]) / n.sims
    decision.table$prob.below.TRP[i] <- signif(raw.TRP, digits = 2)
    if(raw.TRP > 0.995) decision.table$prob.below.TRP[i] <- "> 0.99"
    if(raw.TRP < 0.00995) decision.table$prob.below.TRP[i] <- "< 0.01"
    if(raw.TRP >= 0.00995 & raw.TRP < 0.0995) decision.table$prob.below.TRP[i] <- signif(raw.TRP, digits = 1)
  }
  if(!is.null(RR)) 
  {
    raw.RR <- length(exploit[exploit < RR])/ n.sims
    decision.table$prob.below.RR[i] <- signif(raw.RR, digits = 2)
    if(raw.RR > 0.995) decision.table$prob.below.RR[i] <- "> 0.99"
    if(raw.RR < 0.00995) decision.table$prob.below.RR[i] <- "< 0.01"
    if(raw.RR >= 0.00995 & raw.RR < 0.0995) decision.table$prob.below.RR[i] <- signif(raw.RR, digits = 1)
  }
  
  
}  # End table loop
ndt <- c("Catch (tonnes)", "Exploitation (%)", "Biomass change (%)", "Biomass change (tonnes)", "Probability of Decline")
if(!is.null(LRP)) ndt <- c(ndt,"Probability biomass is below LRP")
if(!is.null(USR)) ndt <- c(ndt,"Probability biomass is below USR")
if(!is.null(TRP)) ndt <- c(ndt,"Probability biomass is below TRP")
if(!is.null(RR)) ndt <- c(ndt,"Probability exploitation is below RR")

names(decision.table) <- ndt

return(dt = decision.table)

} # end function
