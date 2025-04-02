library(BOIN)
library(dplyr)

sim.one.trial = function(trial.id = 1,
                         target = 0.25,
                         p.true = c(0.1, 0.3, 0.5), 
                         ncohort = 10,
                         cohortsize = 3,
                         n.earlystop = 100,
                         startdose = 1,
                         titration = FALSE,
                         p.saf = 0.6 * target, 
                         p.tox = 1.4 * target,
                         cutoff.eli = 0.95,
                         extrasafe = FALSE, 
                         offset = 0.05,
                         boundMTD=FALSE,
                         n.cap = 12,
                         end.backfill = TRUE,
                         n.per.month = 3,
                         dlt.window = 1,
                         p.response.true = c(1, 1, 1),
                         three.plus.three = FALSE){

  
  
  
  ########### get Weibull parameters
  weibull.shape = (log(-log(1 - p.true)) - log(-log(1 - p.true / 2))) / (log(dlt.window) - log(dlt.window / 2))
  weibull.scale = exp(log(dlt.window) - log(-log(1 - p.true)) / weibull.shape)
  
  ########### decision table ##################
  ndose = length(p.true)
  npts = ncohort * cohortsize
  
  if (cohortsize > 1) {
    temp = get.boundary(target, ncohort + ndose * ceiling(n.cap / cohortsize),
                        cohortsize, n.earlystop=ncohort*cohortsize + ndose * n.cap,
                        p.saf, p.tox, cutoff.eli, extrasafe)$full_boundary_tab
  } else {
    temp = get.boundary(target, ncohort + ceiling(n.cap / cohortsize),
                        cohortsize, n.earlystop=ncohort*cohortsize + n.cap,
                        p.saf, p.tox, cutoff.eli, extrasafe)$boundary_tab
  }
  
  dselect = 0
  
  b.e = temp[2, ]
  b.d = temp[3, ]
  b.elim = temp[4, ]
  
  ######### start trial #######################
  y = rep(0, ndose)
  n = rep(0, ndose)
  earlystop = 0
  d = startdose
  elimi = rep(0, ndose)
  
  ##  current dose escalation cohort
  dlt.calendar.cohort = rep(NA, cohortsize)
  dlt.event.cohort = rep(0, cohortsize)
  
  ## backfill patients
  dlt.backfill.calendar = matrix(NA, nrow = n.cap, ncol = length(p.true))
  dlt.backfill.event = matrix(0, nrow = n.cap, ncol = length(p.true))
  
  first.response = rep(Inf, length(p.true))
  first.response.end.cohort = rep(Inf, length(p.true))
  
  clock = 0 
  n.cohort = 1
  
  #check.i = 1
  pat.i = 1
  
  while (n.cohort <= ncohort){  
  #while (n.cohort <= ncohort && check.i <= 4){  
    
    #print(check.i)
    #check.i <- check.i + 1
    
    ## new patient
    arrival.time = rexp(1, rate = n.per.month)
    clock = ifelse(pat.i == 1, 0, clock + arrival.time)
    pat.i = pat.i + 1
    
    ## if current cohort is still open to recruitment
    if (any(is.na(dlt.calendar.cohort))){ 
      # patient is assigned to dose level d 
      
      cohort_pos = max(c(0,which(!is.na(dlt.calendar.cohort)))) + 1
      
      # simulate time to DLT
      dlt.time.star = rweibull(1, shape = weibull.shape[d], scale = weibull.scale[d])
      dlt.event = dlt.time.star < dlt.window
      dlt.time = min(dlt.time.star, dlt.window)
      
      # record calendar time of DLT
      dlt.calendar.cohort[cohort_pos] = clock + dlt.time
      
      # record event
      dlt.event.cohort[cohort_pos] = dlt.event
      
      # simulate response
      if (rbinom(1,1,p.response.true[d]) == 1) first.response[d] = min(first.response[d], clock + dlt.window)
      if (cohort_pos == cohortsize && first.response[d] < Inf) first.response.end.cohort[d] = min(first.response.end.cohort[d], clock + dlt.window)
      
    }
    ## else if current cohort is ready for evaluation...
    else if (max(dlt.calendar.cohort) < clock){
      
      # Officially add the new data to n and y
      n[d] = n[d] + cohortsize
      y[d] = y[d] + sum(dlt.event.cohort)
      
      
      # ... Also collect the information from backfilling
      n.backfill = colSums(dlt.backfill.calendar[,1:d, drop = F] < clock, na.rm = TRUE)
      y.backfill = colSums(dlt.backfill.event[,1:d, drop = F] & dlt.backfill.calendar[,1:d, drop = F] < clock, na.rm = TRUE)
      
      y.total = y.backfill + y[1:d]
      n.total = n.backfill + n[1:d]
      
      
      # ... Also check number of backfill RECRUITED (NOT NECESSARILY REACHED EVALUATION)
      #n.backfill.rec = colSums(dlt.backfill.calendar[,1:d, drop = F] < Inf, na.rm = TRUE)
      n.backfill.rec = n.backfill
      n.total.rec = n.backfill.rec + n[1:d]
      
      
      
      ## ...what does the decision rule say for each dose level
      decision = ifelse(y.total <= b.e[n.total], "escalate",
                        ifelse(y.total >= b.elim[n.total], "eliminate",
                               ifelse(y.total  >= b.d[n.total], "deescalate", "stay")))
      
      
      # special case when 'three.plus.three = TRUE)'
      if (three.plus.three && n.total.rec[d] == 3 && y.total[d] == 1){
         decision[d] = "stay"
      }
      
      
      ## special case d = 1 
      if (d == 1){
        if (decision[d] == "escalate" ){
          if (elimi[d + 1] == 0){
            d = d + 1
          }
          else {
            d = d
            if(n.total.rec[d] >= n.earlystop) break
          }
        }
        else if (decision[d] == "stay"){
          d = d
          if(n.total.rec[d] >= n.earlystop) break
        }
        else if (decision[d] == "eliminate"){
          ## eliminate 
          earlystop = 1
          break
        }
        else if (extrasafe){
          if (d == 1 && n.total[1] >= 3) {
            if (1 - pbeta(target, y.total[1] + 1, n.total[1] - y.total[1] +
                          1) > cutoff.eli - offset) {
              earlystop = 1
              break
            }
          }
        }
        else {
          d = d
          if(n.total.rec[d] >= n.earlystop) break
        }
        
      }
      ## check for conflicting decisions
      else if (decision[d] == "escalate" && all(decision[1:(d-1)] == "escalate")){
        ## escalate
        if (d != ndose && elimi[d + 1] == 0) {
          d = d + 1
        }
        else{
          d = d
          if(n.total.rec[d] >= n.earlystop) break
        }
      }
      else if (decision[d] == "stay" && all(decision[1:(d-1)] %in% c("stay", "escalate"))){
        ## stay
        d = d
        if(n.total.rec[d] >= n.earlystop) break
      }
      else if (decision[d] == "deescalate" && all(decision[1:(d-1)] %in% c("stay", "escalate"))){
        ## deescalate
        d = d - 1
      }
      else if (decision[d] == "eliminate" && all(decision[1:(d-1)] %in% c("stay", "escalate"))){
        ## eliminate
        elimi[d:ndose] = 1
        d = d - 1
        
      }
      else if (decision[d] == "escalate"){
        ## conflict
        ## look for the highest backfilled dose whose data conflict with d
        b.star = max(which(decision[1:(d-1)] != "escalate"))
        
        y.pooled = sum(y.total[b.star:d])
        n.pooled = sum(n.total[b.star:d])
        
        ## if ... then ...
        if (y.pooled  <= b.e[n.pooled]){
          if (d != ndose && elimi[d + 1] == 0) {
            d = d + 1
          }
          else{
            d = d
            if(n.total.rec[d] >= n.earlystop) break
          }
        }
        else if (y.pooled >= b.d[n.pooled]){
          ## de-escalate to the highest dose k that is deemed safe with q_k < lambda_d..
          d.safe = cumsum(y.total[b.star:d])  < b.d[cumsum(n.total[b.star:d])]
          
          if(!any(d.safe)){
            d = b.star - 1
            if (d == 0){
              if (y.total[1] >= b.elim[n.total[1]]){
                earlystop = 1
                break
              }
              else if (extrasafe && (1 - pbeta(target, y[1] + 1, n[1] - y[1] +
                                                1) > cutoff.eli - offset)){
                earlystop = 1
                break
              }
              else{
                d = 1
              }
            }
          }
          else if (b.star - 1 + max(which(d.safe)) == d){
            d = d
            if(n.total.rec[d] >= n.earlystop) break
          }
          else {
            d = b.star - 1 + max(which(d.safe))
          }
          
          
        }
        else {
          d =  d
          if(n.total.rec[d] >= n.earlystop) break
        }
        
      }
      else if (decision[d] %in% c("deescalate", "eliminate", "stay")){
        ## conflict
        if (decision[d] == "eliminate"){
          ## eliminate
          elimi[d:ndose] = 1
        }
        ## look for the highest backfilled dose whose data conflict with d
        b.star = max(which(decision[1:(d-1)] %in% c("deescalate", "eliminate")))
        
        y.pooled = sum(y.total[b.star:d])
        n.pooled = sum(n.total[b.star:d])
        
        ## if ... then ...
        if (y.pooled  <= b.e[n.pooled]){
          if (d != ndose && elimi[d + 1] == 0) {
            d = d + 1
          }
          else{
            d = d
            if(n.total.rec[d] >= n.earlystop) break
          }
        }
        else if (y.pooled >= b.d[n.pooled]){
          ## de-escalate to the highest dose k that is deemed safe with q_k < lambda_d..
          d.safe = cumsum(y.total[b.star:d])  < b.d[cumsum(n.total[b.star:d])]
          
          if(!any(d.safe)){
            d = b.star - 1
            if (d == 0){
              if (y.total[1] >= b.elim[n.total[1]]){
                earlystop = 1
                break
              }
              else if (extrasafe && (1 - pbeta(target, y[1] + 1, n[1] - y[1] +
                                               1) > cutoff.eli - offset)){
                earlystop = 1
                break
              }
              else{
                d = 1
              }
            }
          }
          else if (b.star - 1 + max(which(d.safe)) == d){
            d = d
            if(n.total.rec[d] >= n.earlystop) break
          }
          else {
            d = b.star - 1 + max(which(d.safe))
          }
          
          
        }
        else {
          d =  d
          if(n.total.rec[d] >= n.earlystop) break
        }
      }
      
      n.cohort= n.cohort + 1
      ## Stop because max number of cohorts has been reached
      if (n.cohort > ncohort) break;
  
      
      
      ## start new cohort 
      dlt.calendar.cohort = rep(NA, cohortsize)
      dlt.event.cohort = rep(0, cohortsize)
      
      cohort_pos = max(c(0, which(!is.na(dlt.calendar.cohort)))) + 1
      
      ## simulate DLT time
      # simulate time to DLT
      dlt.time.star = rweibull(1, shape = weibull.shape[d], scale = weibull.scale[d])
      dlt.event = dlt.time.star < dlt.window
      dlt.time = min(dlt.time.star, dlt.window)
      
      # record calendar time of DLT
      dlt.calendar.cohort[cohort_pos] = clock + dlt.time
      
      # record event
      dlt.event.cohort[cohort_pos] = dlt.event
      
      # simulate response
      if (rbinom(1,1,p.response.true[d]) == 1) first.response[d] = min(first.response[d], clock + dlt.window)
      if (cohort_pos == cohortsize && first.response[d] < Inf) first.response.end.cohort[d] = min(first.response.end.cohort[d], clock + dlt.window)
      
      
    }
    else if (d > 1 && clock > min(first.response[1:(d-1)])){
      
      ## check which cohorts are potentially open to backfill
      
      n.backfill = colSums(dlt.backfill.calendar[,1:(d-1), drop = F] < clock, na.rm = TRUE)
      y.backfill = colSums(dlt.backfill.event[,1:(d-1), drop = F] & dlt.backfill.calendar[,1:(d-1), drop = F] < clock, na.rm = TRUE)
      
      n.d = sum(dlt.calendar.cohort < clock, na.rm = TRUE)
      y.d = sum(dlt.event.cohort & dlt.calendar.cohort < clock, na.rm = TRUE)
      
      y.total = y.backfill + y[1:(d-1)]
      n.total = n.backfill + n[1:(d-1)]
      
      y.total.plus = y.total + c(y.backfill[-1], y.d) + y[2:d]
      n.total.plus = n.total + c(n.backfill[-1], n.d) + n[2:d]
      
      ## what do we need to be eligible for backfill????
      cond.1 = first.response[1:(d-1)] < clock
      cond.2 = y.total  < b.d[n.total]
      cond.3 = y.total.plus < b.d[n.total.plus]
      cond.4 = colSums(dlt.backfill.calendar[,1:(d-1), drop = F] < Inf, na.rm = TRUE) + n[1:(d-1)] <= n.cap
      
      safe.doses = cond.2 | cond.3
      if (any(safe.doses == FALSE)) safe.doses[min(which(safe.doses == FALSE)):(d-1)] <- FALSE
      
      open.doses = cond.1 & cond.4 & safe.doses
      
      
      ## allocate to the highest dose open for backfill
      if (any(open.doses)){
        d.backfill = max(which(open.doses))
        
        ## simulate DLT time
        # simulate time to DLT
        dlt.time.star = rweibull(1, shape = weibull.shape[d.backfill], scale = weibull.scale[d.backfill])
        dlt.event = dlt.time.star < dlt.window
        dlt.time = min(dlt.time.star, dlt.window)
        
        # record calendar time of DLT
        backfill.pos = max(c(0, which(!is.na(dlt.backfill.calendar[, d.backfill])))) + 1
        
        dlt.backfill.calendar[backfill.pos, d.backfill] = clock + dlt.time
        dlt.backfill.event[backfill.pos, d.backfill] = dlt.event
        
      }
      
      
      
    }
    else{
      ## refer away the patient
      ## do nothing
    }
    
    
  }
  
  #######################
  ## collect all the data
  n.backfill = colSums(dlt.backfill.calendar[drop = F] < max(dlt.calendar.cohort), na.rm = TRUE)
  #n.backfill = colSums(dlt.backfill.calendar[drop = F] < Inf, na.rm = TRUE)
  y.backfill = colSums(dlt.backfill.event[drop = F] & dlt.backfill.calendar[drop = F] < max(dlt.calendar.cohort), na.rm = TRUE)
  #y.backfill = colSums(dlt.backfill.event[drop = F] & dlt.backfill.calendar[drop = F] < Inf, na.rm = TRUE)
  
  y.total = y.backfill + y
  n.total = n.backfill + n
  
  #######################
  ## apply the mtd logic
  if (earlystop == 1) {
    dselect = 99
  }
  else {
    dselect = select.mtd(target, n.total, y.total, cutoff.eli,
                         extrasafe, offset, boundMTD = boundMTD, p.tox=p.tox)$MTD
  }
  
  
  max.t = max(dlt.calendar.cohort)
  
  out.df <- data.frame(trial.id = trial.id,
                       d = 1:ndose,
                       n = n.total,
                       y = y.total,
                       dselect = dselect,
                       max.t = max.t)
  
  return(out.df)
  
  
}

#############################################################################
get.oc.bf <- function(ntrial = 1000,
                      seed = 3262,
                      target = 0.25,
                      p.true = c(0.1, 0.3, 0.5), 
                      ncohort = 10,
                      cohortsize = 3,
                      n.earlystop = 100,
                      startdose = 1,
                      titration = FALSE,
                      p.saf = 0.6 * target, 
                      p.tox = 1.4 * target,
                      cutoff.eli = 0.95,
                      extrasafe = FALSE, 
                      offset = 0.05,
                      boundMTD=FALSE,
                      n.cap = 12,
                      end.backfill = TRUE,
                      n.per.month = 3,
                      dlt.window = 1,
                      p.response.true = c(1, 1, 1),
                      three.plus.three = FALSE){
  
  ############ Sanity #############
  if (target < 0.05) {
    stop("the target is too low!")
  }
  if (target > 0.6) {
    stop("the target is too high!")
  }
  if ((target - p.saf) < (0.1 * target)) {
    stop("the probability deemed safe cannot be higher than or too close to the target!")
  }
  if ((p.tox - target) < (0.1 * target)) {
    stop("the probability deemed toxic cannot be lower than or too close to the target!")
  }
  if (offset >= 0.5) {
    stop("the offset is too large!")
  }
  if (n.earlystop <= 6) {
    warning("the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18.")
  }
  ### BF-BOIN specific ####
  if (cohortsize < 3) {
    stop("not considering cohortsize less than 3 for this implementation of BF-BOIN")
  }
  if (titration == TRUE) {
    stop("not considering titration for this implementation of BF-BOIN")
  }
  if (three.plus.three == TRUE && target != 0.25){
    stop("3 plus 3 option is only available when the target is 0.25")
  }
  if (three.plus.three == TRUE && cohortsize != 3){
    stop("3 plus 3 option is only available when the cohortsize is 3")
  }

  set.seed(seed)
 
  res = purrr::map_df(1:ntrial, sim.one.trial,
                      target = target,
                      p.true = p.true, 
                      ncohort = ncohort,
                      cohortsize = cohortsize,
                      n.earlystop = n.earlystop,
                      startdose = startdose,
                      titration = titration,
                      p.saf = p.saf, 
                      p.tox = p.tox,
                      cutoff.eli = cutoff.eli,
                      extrasafe = extrasafe, 
                      offset = offset,
                      boundMTD = boundMTD,
                      n.cap = n.cap,
                      end.backfill = end.backfill,
                      n.per.month = n.per.month,
                      dlt.window = dlt.window,
                      p.response.true = p.response.true,
                      three.plus.three = three.plus.three)
  
  
  ndose = length(p.true)
  
  res_sum = res |> group_by(d) |> summarise(nptsdose = mean(n), ntoxdose = mean(y))
  nptsdose = res_sum |> pull(nptsdose)
  ntoxdose = res_sum |> pull(ntoxdose)
  dselect = res |> group_by(trial.id) |> slice(1) |> pull(dselect)  
  max.t = res |> group_by(trial.id) |> slice(1) |> pull(max.t)  
  selpercent = rep(0, ndose)
  for (i in 1:ndose) {
    selpercent[i] = sum(dselect == i)/ntrial * 100
  }
  
  out = list(selpercent = selpercent, 
             npatients = nptsdose,
             percentpatients = nptsdose / sum(nptsdose) * 100,
             ntox = ntoxdose, 
             totaltox = sum(res$y) / ntrial, 
             totaln = sum(res$n)/ntrial,
             percentstop = sum(dselect == 99) / ntrial * 100,
             duration = mean(max.t))
  
  out
  
}

#############################################################################
get.oc.bf(ntrial = 1000,
          seed = 9,
          target = 0.25,
          p.true = c(0.001, 0.001), 
          ncohort = 10,
          cohortsize = 3,
          n.earlystop = 6,
          startdose = 1,
          titration = FALSE,
          cutoff.eli = 0.95,
          extrasafe = FALSE, 
          offset = 0.1,
          boundMTD=FALSE,
          n.cap = 12,
          end.backfill = TRUE,
          n.per.month = 1,
          dlt.window = 1,
          p.response.true = c(0.001, 0.001))
  

#############################################################################
get.oc.bf(ntrial = 1000,
          seed = 9,
          target = 0.25,
          p.true = c(0.1, 0.5), 
          ncohort = 10,
          cohortsize = 3,
          n.earlystop = 9,
          startdose = 1,
          titration = FALSE,
          cutoff.eli = 0.95,
          extrasafe = TRUE, 
          offset = 0.1,
          boundMTD=FALSE,
          n.cap = 12,
          end.backfill = TRUE,
          n.per.month = 1,
          dlt.window = 1,
          p.response.true = c(0.001, 0.001))


#############################################################################
get.oc.bf(ntrial = 1000,
          seed = 9,
          target = 0.25,
          p.true = c(0.1, 0.5), 
          ncohort = 10,
          cohortsize = 3,
          n.earlystop = 9,
          startdose = 1,
          titration = FALSE,
          cutoff.eli = 0.95,
          extrasafe = TRUE, 
          offset = 0.1,
          boundMTD=FALSE,
          n.cap = 12,
          end.backfill = TRUE,
          n.per.month = 1,
          dlt.window = 1,
          p.response.true = c(0.001, 0.001),
          three.plus.three = TRUE)










