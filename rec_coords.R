library(readr)
library(data.table)

rec_num = 289

coordinates_rec_latlon_id_gcy <- read_delim("Documents/Arjun/coordinates_rec_latlon_id_gcy.csv", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
coordinates_receivers_recnum <- read_delim("Documents/Arjun/coordinates_receivers_recnum.csv","\t", escape_double = FALSE, col_names = FALSE,col_types = cols(X1 = col_skip()), trim_ws = TRUE)

coordinates_rec_latlon_id_gcy <- data.table(coordinates_rec_latlon_id_gcy)
coordinates_receivers_recnum <- data.table(coordinates_receivers_recnum)
  
coordinates_rec_latlon_id_gcy <- coordinates_rec_latlon_id_gcy[, !c('X3','X4')]
coordinates_receivers_recnum <- coordinates_receivers_recnum[, !c('X3','X4')]

setkey(coordinates_receivers_recnum,X2)
setkey(coordinates_rec_latlon_id_gcy,X2)

Result <- merge(coordinates_receivers_recnum,coordinates_rec_latlon_id_gcy, all.x=TRUE)[0:rec_num,]
result <- Result[rep(seq_len(nrow(Result)), each = length(x)-1), ]
result

# The angle sequence
x <- seq(0,360,10)

# Dummy SNR data

y <- replicate(length(x)-1, runif(rec_num)) 

snr_list <- as.vector(t(y)) 
length(snr_list)
result$snr = snr_list
result$radius = 0.4
result$startAngle <- rep(x[1:length(x)-1], each=1, times=rec_num)
result$endAngle = rep(x[2:length(x)], each=1, times=rec_num)
result
write.csv(result, '/home/shamik/Documents/Arjun/result.csv')
