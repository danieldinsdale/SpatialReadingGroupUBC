######################################################
#TEMPORARY FILE TO FIGURE OUT DATA STRUCTURE #########
######################################################
# original column names
# [1] "Geevor Tin Mine -- Simms Lode West                                              ²"
# [2] "5445 5 1 -9"                                                                      
# [3] "X co-ordinate       Y co-ordinate       Thickness (inches)  Grade lbs/ton SnO2"   
# [4] "Date code"
# read in data skipping first 4 lines (column names)
tinData <- read.table("geevor.dat", skip=4)
# manually trying to figure out name for each column
colnames(tinData) <- c("", "", "", "XCoord", "YCoord", "Thickness", "Grade", "DateCode")
# testing columns/coordinates (looks like a Donkey Kong level)
plot(x=tinData$XCoord, y=tinData$YCoord)
