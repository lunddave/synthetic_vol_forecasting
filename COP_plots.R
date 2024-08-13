
COP <- read.csv('COP.csv')


par(mfrow = c(3,1))
## Dates ##

# We will now use  2008-09-24 00:00:00
# Training period begins  2008-03-08 00:00:00
# Period to predict begins  2008-09-25 00:00:00
# Prediction period ends  2008-10-02 00:00:00
donor1 <- COP[COP$Date >= '2008-03-08' & COP$Date <= '2008-09-25',]

plot(y = donor1$Close, x = as.Date(donor1$Date), type = 'l')

abline(v = as.Date('2008-09-24'), col = "pink")

#
#
# We will now use  2014-11-26 00:00:00
# Training period begins  2014-05-10 00:00:00
# Period to predict begins  2014-11-27 00:00:00
# Prediction period ends  2014-12-04 00:00:00
donor2 <- COP[COP$Date >= '2014-05-10' & COP$Date <= '2014-11-27',]

plot(y = donor2$Close, x = as.Date(donor2$Date), type = 'l')
abline(v = as.Date('2014-11-26'), col = "pink")



# We will now use  2020-03-06 00:00:00
# Training period begins  2019-08-19 00:00:00
# Period to predict begins  2020-03-07 00:00:00
# Prediction period ends  2020-03-14 00:00:00

donor3 <- COP[COP$Date >= '2019-08-19' & COP$Date <= '2020-03-07',]

plot(y = donor3$Close, x = as.Date(donor3$Date), type = 'l')
abline(v = as.Date('2020-03-06'), col = "pink")

