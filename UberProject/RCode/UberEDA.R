library(ggplot2) # Data visualization
setwd("~/UBC/SpatialReadingGroupUBC/UberProject/Data")
# Input data files are available in the "../input/" directory.
# For example, running this (by clicking run or pressing Shift+Enter) will list the files in the input directory
# Any results you write to the current directory are saved as output.

# using code from kaggle to read in data
uber = read.csv('uber-raw-data-sep14.csv', stringsAsFactors = F)

uber$Date = sapply(strsplit(uber$Date.Time, split = " "), function(x) x[[1]][1])
uber$Date = as.Date(uber$Date, format = "%m/%d/%Y")
uber = uber[order(uber$Date),]

min_long = min(uber$Lon)
max_long = max(uber$Lon)
min_lat = min(uber$Lat)
max_lat = max(uber$Lat)

#plot(x=uber$Lat, y=uber$Lon)

# code from kaggle
m <- ggplot(data=uber,aes(Lon,Lat)) + 
  geom_point(size=0.06, color="white", alpha = 0.2) +
  scale_x_continuous(limits=c(min_long, max_long)) +
  scale_y_continuous(limits=c(min_lat, max_lat)) +
  theme(panel.background = element_rect(fill = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
        