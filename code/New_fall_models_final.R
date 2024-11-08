#final set of new fall models


####CDDP####
#Chilling Degree Day model with photoperiod adapted from PTT
#Existing model from Schadel et al., 2023

CDDP <- function(par, data){
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1]) # int
  T_base <- par[2]
  F_crit <- par[3]
  
  # create forcing/chilling rate vector
  # forcing
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  Rf <- ((1 - data$Li) / 24) * Rf
  Rf[1:t0,] <- 0
  
  # DOY of fall senescence criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}


####CDDs####
#Chilling Degree Day (CDD) adapted from
#with a sigmoidal temperature response
#Existing model from Schadel et al., 2023

CDDs <- function(par, data){
  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  t0 <- round(par[1])
  b <- par[2]
  c <- par[3]
  F_crit <- par[4]
  
  # sigmoid temperature response
  Rf <- 1 / (1 + exp(-b * (data$Tmini - c)))
  Rf[1:t0,] <- 0
  
  # DOY of fall senescence criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####CDD2####
#t0 = SOS

CDD2 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 2){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base <- par[1]
  F_crit <- par[2]
  
  # create forcing/chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any precip before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    Rf[1:SOS,i] <- 0
  }
  
  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}



####CDD3####
#t0 = SOS + dt

CDD3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base <- par[1]
  F_crit <- par[2]
  dt <- round(par[3])
  
  # create forcing/chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any precip before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}



####CDDs3####
#t0 = SOS + dt

CDDs3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  b <- par[1]
  c <- par[2]
  F_crit <- par[3]
  dt <- round(par[4])
  
  # sigmoid temperature response
  Rf <- 1 / (1 + exp(-b * (data$Tmini - c)))
  
  #Set any precip before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  
  # DOY of fall senescence criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}



####SOS####
#EOS is certain number of days after SOS

SOS <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 1){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  dt <- round(par[1])
  
  #SOS plus # of days
  doy <- data$SOS_dates + dt
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}



####DD_W3####
#brown-down occurs after enough consecutive dry days (DD)
#t0 = SOS + dt

DD_W3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  P_crit <- par[1]
  P_base <- par[2]
  dt <- round(par[3])
  
  #Precip requirement
  Rw <- data$Pi
  
  ##Calculate cdd
  
  #If precip = 0, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw > P_base, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)
  
  #force days before t0 to be 0 so not counted in cdd
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      k_cdd[1:t0,i] <- 0}
  }
  
  
  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }
  
  #Pull out years
  col_names <- names(k_cdd)
  
  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))
  
  #loop years (columns) through function
  for (i in 1:length(col_names)) {
    
    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output
    
  }
  
  # DOY of budburst criterium
  doy <- apply(cdd_matrix,2, function(xt){
    data$doy[which(xt >= P_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####DD_SM3####
#t0 = SOS + dt

DD_SM3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  SM_base <- par[1]
  SM_crit <- par[2]
  dt <- round(par[3])
  
  #Precip requirement
  Rw <- data$SM
  
  ##Calculate cdd
  
  #If SM < SM_base, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw > SM_base, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)
  
  #force days before t0 to be 0 so not counted in cdd
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      k_cdd[1:t0,i] <- 0}
  }
  
  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }
  
  #Pull out years
  col_names <- names(k_cdd)
  
  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))
  
  #loop years (columns) through function
  for (i in 1:length(col_names)) {
    
    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output
    
  }
  
  # DOY of budburst criterium
  doy <- apply(cdd_matrix,2, function(xt){
    data$doy[which(xt >= SM_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####DD_VPD3#### 
#t0 = SOS + dt

DD_VPD3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  VPD_base <- par[1]
  VPD_crit <- par[2]
  dt <- round(par[3])
  
  #Precip requirement
  Rw <- data$VPDi
  
  ##Calculate cdd
  
  #If precip = 0, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw < VPD_base, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)
  
  #force days before t0 to be 0 so not counted in cdd
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      k_cdd[1:t0,i] <- 0}
  }
  
  #function to add up cdd, resetting when VPD above VPD base
  cdd_func <- function(col){
    as.matrix(with(k_cdd, ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }
  
  #Pull out years
  col_names <- names(k_cdd)
  
  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))
  
  #loop years (columns) through function
  for (i in 1:length(col_names)) {
    
    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output
    
  }
  
  # DOY of budburst criterium
  doy <- apply(cdd_matrix,2, function(xt){
    data$doy[which(xt >= VPD_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}




####W3####
#t0 = SOS + dt

W3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  a <- par[1]
  P_crit <- par[2]
  dt <- round(par[3])
  dp <- round(par[4])
  
  #Make output matrix
  rows <- nrow(data$Pi)
  cols <- ncol(data$Pi)
  Rw <- matrix(0,rows,cols)
  k <- c(1:floor(dp))
  W <- 0
  
  #Make precip matrix
  Pi_matrix <- data$Pi
  
  #Set any precip before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  # loop through to find precip over last "dp" days with decay function 
  
  for (i in 1:cols){
    
    for (j in 1:rows){
      
      t0 <- SOS_day[i] + dt
      
      if (j <= t0) {
        next
        
      }else{
        
        if(j > dp){
          
          (my_vector <- Pi_matrix[j-k,i]*(a^k))
          W <- sum(my_vector) + Pi_matrix[j,i]
        }
      }
      
      if (W <= P_crit){
        Rw[j,i] <- 1000
        break
      }
    }
  }
  
  
  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(xt == 1000)[1]]
  })
  
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####SM3####
#t0 = SOS + dt

SM3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  SM_base <- par[1]
  SM_crit <- par[2]
  dt <- round(par[3])
  
  # create forcing/chilling rate vector
  Rw <- data$SM - SM_base
  Rw[Rw > 0] <- 0
  
  #Set any precip before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rw[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) <= SM_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}



####VPD3####
#t0 = SOS + dt

VPD3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  VPD_base <- par[1]
  VPD_crit <- par[2]
  dt <- round(par[3])
  
  # create forcing/chilling rate vector
  Rw <- data$VPDi - VPD_base
  Rw[Rw < 0] <- 0
  
  #Set any VPD before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if(t0 <= 365){
      Rw[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= VPD_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}



####CDD_DD_W3####
#t0 = SOS + dt

CDD_DD_W3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  T_base <- par[1]
  F_crit <- par[2]
  P_crit <- par[3]
  P_base <- par[4]
  dt <- round(par[5])
  
  #Precip requirement
  Rw <- data$Pi
  
  ##Calculate cdd
  
  #If precip = 0, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw > P_base, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)
  
  #force days before t0 to be 0 so not counted in cdd
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      k_cdd[1:t0,i] <- 0}
  }
  
  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }
  
  #Pull out years
  col_names <- names(k_cdd)
  
  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))
  
  #loop years (columns) through function
  for (i in 1:length(col_names)) {
    
    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output
    
  }
  
  # create forcing/chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any temp before t0 to 0
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy_1 <- apply(cdd_matrix,2, function(xt){
    data$doy[which(xt >= P_crit)[1]]})
  
  doy_2 <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  #choose earlier date for each year
  doy <- pmin(doy_1, doy_2, na.rm=TRUE)
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####CDD_DD_SM3####
#t0 = SOS + dt

CDD_DD_SM3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  T_base <- par[1]
  F_crit <- par[2]
  SM_base <- par[3]
  SM_crit <- round(par[4])
  dt <- round(par[5])
  
  #Precip requirement
  Rw <- data$SM
  
  ##Calculate cdd
  
  #If precip = 0, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw > SM_base, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)
  
  #force days before t0 to be 0 so not counted in cdd
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      k_cdd[1:t0,i] <- 0}
  }
  
  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }
  
  #Pull out years
  col_names <- names(k_cdd)
  
  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))
  
  #loop years (columns) through function
  for (i in 1:length(col_names)) {
    
    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output
    
  }
  
  # create forcing/chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any temp before t0 to 0
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy_1 <- apply(cdd_matrix,2, function(xt){
    data$doy[which(xt >= SM_crit)[1]]})
  
  doy_2 <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  #choose earlier date for each year
  doy <- pmin(doy_1, doy_2)
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####CDD_DD_VPD3####
#t0 = SOS + dt

CDD_DD_VPD3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  T_base <- par[1]
  F_crit <- par[2]
  VPD_base <- par[3]
  VPD_crit <- round(par[4])
  dt <- round(par[5])
  
  #Precip requirement
  Rw <- data$VPDi
  
  ##Calculate cdd
  
  #If precip = 0, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw < VPD_base, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)
  
  #force days before t0 to be 0 so not counted in cdd
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      k_cdd[1:t0,i] <- 0}
  }
  
  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }
  
  #Pull out years
  col_names <- names(k_cdd)
  
  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))
  
  #loop years (columns) through function
  for (i in 1:length(col_names)) {
    
    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output
    
  }
  
  # create forcing/chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any temp before t0 to 0
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy_1 <- apply(cdd_matrix,2, function(xt){
    data$doy[which(xt >= VPD_crit)[1]]})
  
  doy_2 <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  #choose earlier date for each year
  doy <- pmin(doy_1, doy_2)
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####CDD_W3####
#t0 = SOS + dt

CDD_W3 <- function(par, data ){
  
  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  T_base <- par[1]
  a <- par[2]
  F_crit <- par[3]
  P_crit <- par[4]
  dt <- round(par[5])
  dp <- round(par[6])
  
  # create chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any temp before t0 to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  
  #Calculate precip#
  
  #Make output matrix
  rows <- nrow(data$Pi)
  cols <- ncol(data$Pi)
  Rw <- matrix(0,rows,cols)
  k <- c(1:floor(dp))
  W <- 0
  
  #Make precip matrix
  Pi_matrix <- data$Pi
  
  # loop through to find precip over last "dp" days with decay function 
  
  for (i in 1:cols){
    
    for (j in 1:rows){
      
      t0 <- SOS_day[i] + dt
      
      if (j <= t0) {
        next
        
      }else{
        
        if(j > dp){
          
          (my_vector <- Pi_matrix[j-k,i]*(a^k))
          W <- sum(my_vector) + Pi_matrix[j,i]
        }
      }
      
      if (W <= P_crit){
        Rw[j,i] <- 1000
        break
      }
    }
  }
  
  # DOY of budburst criterium
  doy_1 <- apply(Rw,2, function(xt){
    data$doy[which(xt == 1000)[1]]
  })
  
  # DOY of budburst criterium
  doy_2 <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  #choose earlier date for each year
  doy <- pmin(doy_1, doy_2, na.rm=TRUE)
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####CDD_SM3####

CDD_SM3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  T_base <- par[1]
  F_crit <- par[2]
  SM_base <- par[3]
  SM_crit <- par[4]
  dt <- round(par[5])
  
  # create VPD forcing rate vector
  Rw <- data$SM - SM_base
  Rw[Rw > 0] <- 0
  
  #Set any VPD before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rw[1:t0,i] <- 0}
  }
  
  
  # create chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any temp before t0 to 0
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy_1 <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) <= SM_crit)[1]]
  })
  
  doy_2 <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  #choose earlier date for each year
  doy <- pmin(doy_1, doy_2, na.rm=TRUE)
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####CDD_VPD3####
#t0 = SOS + dt

CDD_VPD3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  T_base <- par[1]
  F_crit <- par[2]
  VPD_base <- par[3]
  VPD_crit <- par[4]
  dt <- round(par[5])
  
  # create VPD forcing rate vector
  Rw <- data$VPDi - VPD_base
  Rw[Rw < 0] <- 0
  
  #Set any VPD before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rw[1:t0,i] <- 0}
  }
  
  
  # create chilling rate vector
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any temp before t0 to 0
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  # DOY of budburst criterium
  doy_1 <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= VPD_crit)[1]]
  })
  
  doy_2 <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  #choose earlier date for each year
  doy <- pmin(doy_1, doy_2, na.rm=TRUE)
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}


####W_T_bal3####
#t0 = SOS + dt


W_T_bal3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  a <- par[1]
  F_crit <- par[2]
  dt <- round(par[3])
  
  
  #Precip requirement
  Rw <- data$Pi - (a*data$Ti)
  
  #Set any precip before SOS to 20 (so not triggered)
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rw[1:t0,i] <- 0}
  }
  
  
  # DOY of fall senescence criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) <= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####CDDP3####
#t0 = SOS + dt
#adjusted so same direction of influence on Rf for darker and colder

CDDP3 <- function(par, data){
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base <- par[1]
  F_crit <- par[2]
  dt <- round(par[3])
  
  # create forcing/chilling rate vector
  # forcing
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  #Set any precip before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  #temper Rf based on Li
  Rf <- -1 * ((1-(data$Li/24)) * Rf)
  
  # DOY of fall senescence criterium
  doy <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



####SMP3####
#t0 = SOS + dt

SMP3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  SM_base <- par[1]
  SM_crit <- par[2]
  dt <- round(par[3])
  
  # create forcing/chilling rate vector
  Rw <- data$SM - SM_base
  Rw[Rw > 0] <- 0
  
  #Set any precip before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rw[1:t0,i] <- 0}
  }
  
  #temper Rf based on Li
  Rw <- (1-(data$Li/24)) * Rw
  
  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) <= SM_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}



####VPDP3####
#VPD & photoperiod

VPDP3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  VPD_base <- par[1]
  F_crit <- par[2]
  dt <- round(par[3])
  
  # create forcing/chilling rate vector
  Rw <- data$VPD - VPD_base
  Rw[Rw < 0] <- 0
  
  #Set any VPD before SOS to 0
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rw[1:t0,i] <- 0}
    
  }
  
  
  #temper Rw based on Li
  Rw <- (1-(data$Li/24)) * Rw
  
  # DOY of budburst criterium
  doy <- apply(Rw,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
  
}


####CDDP_DD_W3####
#t0 = SOS + dt

CDDP_DD_W3 <- function(par, data){
  
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  #parameters
  T_base <- par[1]
  F_crit <- par[2]
  P_crit <- par[3]
  P_base <- par[4]
  dt <- round(par[5])
  
  #Precip requirement
  Rw <- data$Pi
  
  ##Calculate cdd
  
  #If precip = 0, then assign a "1" to the cell
  k_cdd <- data.frame(ifelse(Rw > P_base, 0, 1))
  colnames(k_cdd) <- paste0(data$site, "_", data$year)
  
  #force days before t0 to be 0 so not counted in cdd
  SOS_day <- as.vector(data$SOS_dates)
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      k_cdd[1:t0,i] <- 0}
  }
  
  #function to add up cdd, resetting when rain occurs
  cdd_func <- function(col){
    as.matrix(with(k_cdd, ave(k_cdd[[col]], cumsum(k_cdd[[col]] == 0), FUN = cumsum)))
  }
  
  #Pull out years
  col_names <- names(k_cdd)
  
  #Make empty matrix
  cdd_matrix <- matrix(NA, nrow = 365, ncol = length(col_names))
  
  #loop years (columns) through function
  for (i in 1:length(col_names)) {
    
    year <- col_names[i]
    output <- cdd_func(col = year)
    cdd_matrix[,i] <- output
    
  }
  
  # create forcing/chilling rate vector
  # forcing
  Rf <- data$Tmini - T_base
  Rf[Rf > 0] <- 0
  
  
  #Set any precip before SOS to 0
  
  for (i in 1:length(SOS_day)){
    
    SOS <- SOS_day[i]
    t0 <- SOS + dt
    
    if (t0 < 365){
      Rf[1:t0,i] <- 0}
  }
  
  #temper Rf based on Li
  Rf <- -1 * ((1-(data$Li/24)) * Rf)
  
  # DOY of budburst criterium
  doy_1 <- apply(cdd_matrix,2, function(xt){
    data$doy[which(xt >= P_crit)[1]]})
  
  doy_2 <- apply(Rf,2, function(xt){
    data$doy[which(cumsum(xt) >= F_crit)[1]]
  })
  
  #choose earlier date for each year
  doy <- pmin(doy_1, doy_2, na.rm=TRUE)
  
  # set export format, either a rasterLayer
  # or a vector
  shape_model_output(data = data, doy = doy)
}



