library(readr)
library(data.table)
library(dplyr)
    
rec_num = 150

coordinates_rec_latlon_id_gcy <- read_delim("../coordinates_rec_latlon_id_gcy.csv", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
coordinates_receivers_recnum <- read_delim("../coordinates_receivers_recnum.csv","\t", escape_double = FALSE, col_names = FALSE,col_types = cols(X1 = col_skip()), trim_ws = TRUE)

coordinates_rec_latlon_id_gcy <- data.table(coordinates_rec_latlon_id_gcy)
coordinates_receivers_recnum <- data.table(coordinates_receivers_recnum)

coordinates_rec_latlon_id_gcy <- coordinates_rec_latlon_id_gcy[, !c('X3','X4')]
coordinates_receivers_recnum <- coordinates_receivers_recnum[, !c('X3','X4')]

#setkey(coordinates_receivers_recnum,X2)
#setkey(coordinates_rec_latlon_id_gcy,X2)

x <- seq(0,360,10)
Result <- left_join(coordinates_receivers_recnum,coordinates_rec_latlon_id_gcy)

# select n random samples 
#load('recs.RData')
rn = sample(seq(1,289),rec_num)
Result = Result[rn,]

result <- Result[rep(seq_len(nrow(Result)), each = length(x)-1), ]
result

# Dummy SNR data
for (dr in dir(pattern = "stack")) {
  y <- read_csv(sprintf("%s/avgsnr.txt", dr), col_names = FALSE)
  y = y[rn,]
  
  #y <- replicate(length(x)-1, runif(rec_num))
  
  snr_list <- as.vector(t(y))
  length(snr_list)
  result$snr = snr_list
  result$radius = 0.7
  result <- within(result, radius[is.nan(snr)] <- 0.0)
  result <- within(result, snr[is.nan(snr)] <- 0.0)
  #result <- within(result, snr[snr < 1.0] <- 1.0)
  #result <- within(result, snr <- log(snr))
  
  result$startAngle <- rep(x[1:length(x)-1], each=1, times=rec_num)
  result$endAngle = rep(x[2:length(x)], each=1, times=rec_num)
  result
  print(max(result$snr))
  write.table(result, sprintf('%s/result.csv', dr), col.names = FALSE)
}