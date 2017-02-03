# Tin Data analysis #

# Load required packages #
library(INLA)
library(RandomFields)
library(splancs)
library(geoR)
library(fields)
library(ggplot2)


# read in data #
tin = read.table('geevor.dat', skip = 4)

colnames(tin)[4:8]=c('Xcoord', 'Ycoord', 'Thick_inch',  'Grade', 'Date')
tin$Date = as.factor(tin$Date)

summary(tin) # Zero grade is found

# X coordinates and Y coordinates are sampled often at integer spacings so a fine discretized grid needed #

# Define a 'Total Tin' variable equal to the product of the vein thickness and grade #
tin$Total = tin$Thick_inch * tin$Grade

hist(tin$Grade)  
hist(tin$Thick_inch)
hist(tin$Total)
# Log transform of grade required to reduce the range, Thickness is slightly right-skew #       

tin$Grade = log(tin$Grade + 1)
tin$Thick_inch = log(tin$Thick_inch)
tin$Total = log(tin$Total + 1)

hist(tin$Grade)  
hist(tin$Thick_inch)
hist(tin$Total)
# Much better behaved #

tin_dat = as.geodata(tin,coords.col=4:5,data.col=c(2,3,6,7,8,9))
dup.coords(tin_dat) # duplicate measurements found - good estimate of nugget effect? #
tin_dat2 = jitterDupCoords(tin_dat, max = 0.1)

par(mfrow=c(1,1))
bubble(SpatialPointsDataFrame(coords = tin_dat$coords, data = data.frame(exp(as.numeric(tin_dat$data[,'Grade'])))), maxsize = 1)
bubble(SpatialPointsDataFrame(coords = tin_dat$coords, data = data.frame(exp(as.numeric(tin_dat$data[,'Thick_inch'])))), maxsize = 1)

bubble_data_grade = data.frame(x = tin_dat$coords[,1], y = tin_dat$coords[,2],
                            z = exp(as.numeric(tin_dat$data[,'Grade'])))
bubble.plot_grade = ggplot(data = bubble_data_grade) + geom_point(aes(x=x, y=y, size=z)) +
  #geom_polygon(data = data.frame(poly_x = galicia.poly[,1], poly_y = galicia.poly[,2]), aes(x = poly_x, y = poly_y, alpha = 0)) +
  scale_size_continuous(range = c(0, 5)) +
  xlab('') + ylab('') + ggtitle('Sampled Grade') +
  theme(legend.position = "none")
bubble.plot_grade
       
bubble_data_thick = data.frame(x = tin_dat$coords[,1], y = tin_dat$coords[,2],
                               z = exp(as.numeric(tin_dat$data[,'Thick_inch'])))
bubble.plot_thick = ggplot(data = bubble_data_thick) + geom_point(aes(x=x, y=y, size=z)) +
  #geom_polygon(data = data.frame(poly_x = galicia.poly[,1], poly_y = galicia.poly[,2]), aes(x = poly_x, y = poly_y, alpha = 0)) +
  scale_size_continuous(range = c(0, 5)) +
  xlab('') + ylab('') + ggtitle('Sampled Vein Thickness') +
  theme(legend.position = "none")
bubble.plot_thick

bubble_data_total = data.frame(x = tin_dat$coords[,1], y = tin_dat$coords[,2],
                               z = exp(as.numeric(tin_dat$data[,'Total'])))
bubble.plot_total = ggplot(data = bubble_data_total) + geom_point(aes(x=x, y=y, size=z)) +
  #geom_polygon(data = data.frame(poly_x = galicia.poly[,1], poly_y = galicia.poly[,2]), aes(x = poly_x, y = poly_y, alpha = 0)) +
  scale_size_continuous(range = c(0, 5)) +
  xlab('') + ylab('') + ggtitle('Sampled Total tin') +
  theme(legend.position = "none")
bubble.plot_total

# Compute Sample variograms #
vario_Grade<-variog(coords = tin_dat$coords, data = as.numeric(tin_dat$data[,'Grade']))
plot(vario_Grade,pch=19) 

vario_Thick<-variog(coords = tin_dat$coords, data = as.numeric(tin_dat$data[,'Thick_inch']))
plot(vario_Thick,pch=19)

vario_Total<-variog(coords = tin_dat$coords, data = as.numeric(tin_dat$data[,'Total']))
plot(vario_Total,pch=19)

# All three variogram plots show two interesting features :
# a) The replications allow for the estimation of the nugget effect at distance of zero #
#    Thus measurement error is present and large.
# b) The shapes of the variograms are all concave thus a rough spatial surface is predicted 
#    Kappa (roughness parameter) of the Matern will be 0.5 or 1 - equivalently exponential
#    correlation function.

# Fit the model with Maximum likelihood - Warning takes a long time #
fit_ML_Grade<-likfit(coords = tin_dat2$coords, data = as.numeric(tin_dat2$data[,'Grade']),cov.model="mat",kappa=0.5,ini=c(0.2047,0.0595),nugget=0.4)
summary(fit_ML_Grade)
# sigmasq (partial sill) =  0.3864, phi (range parameter)  =  0.0601, (estimated) nugget =  0.9438, beta =  3.3413 
# re-run to check convergence
fit_ML_Grade<-likfit(coords = tin_dat$coords, data = as.numeric(tin_dat$data[,'Grade']),cov.model="mat",kappa=0.5,ini=c(0.8,0.2),nugget=0.1)
summary(fit_ML_Grade)

tausq_Grade<-fit_ML_Grade$tausq
sigmasq_Grade<-fit_ML_Grade$sigmasq
phi_Grade<-fit_ML_Grade$phi
vario.x_Grade<-seq(from = 0, to = 2300, length.out = 100000)
vario.y_Grade<-tausq_Grade+sigmasq_Grade*(1-exp(-vario.x_Grade/phi_Grade))
plot(vario_Grade,pch=19,xlab="u",ylab="V(u)", main = "Variogram of Grade data")


fit_ML_Thick<-likfit(coords = tin_dat$coords, data = as.numeric(tin_dat$data[,'Thick_inch']),cov.model="mat",kappa=0.5,ini=c(0.2047,0.0595),nugget=0.4)
summary(fit_ML_Thick)
# re-run to check convergence
fit_ML_Thick<-likfit(coords = tin_dat$coords, data = as.numeric(tin_dat$data[,'Thick_inch']),cov.model="mat",kappa=0.5,ini=c(0.8,0.2),nugget=0.1)
summary(fit_ML_Thick)

fit_ML_Total<-likfit(coords = tin_dat$coords, data = as.numeric(tin_dat$data[,'Total']),cov.model="mat",kappa=0.5,ini=c(0.2047,0.0595),nugget=0.4)
summary(fit_ML_Total)
# re-run to check convergence
fit_ML_Total<-likfit(coords = tin_dat$coords, data = as.numeric(tin_dat$data[,'Total']),cov.model="mat",kappa=0.5,ini=c(0.8,0.2),nugget=0.1)
summary(fit_ML_Total)


########################### INLA ########################
###define the dimension of the grid
nrow=100
ncol=100
n=nrow*ncol
xrange=range(tin$Xcoord)
yrange=range(tin$Ycoord)
cutpointsx=seq((xrange[1]-0.1),(xrange[2]+0.1),length.out=nrow)
cutpointsy=seq((yrange[1]-0.1),(yrange[2]+0.1),length.out=ncol)

###work on lead97
grid1.97 <- as.numeric(cut(tin$Xcoord,cutpointsx))
grid2.97 <- as.numeric(cut(tin$Ycoord,cutpointsy))
grid97=cbind(grid1.97,grid2.97)
node97=numeric(length(tin$Xcoord))
for(i in 1: length(tin$Xcoord))
  node97[i]=inla.lattice2node(grid97[i,1],grid97[i,2],nrow=nrow,ncol=ncol)

# Put duplicated coordinates into a new data set #
count = 0
#node = rep(0 times = length(tin$Xcoord))
#X = rep(0, times = length(tin$Xcoord))
#Y = rep(0, times = length(tin$Ycoord))

node_temp = node97
X_temp = tin$Xcoord
Y_temp = tin$Ycoord
ind=T
count = 1
y_total = matrix(rep(NA,n), nrow=n, ncol=1)
temp_Total = tin$Total
while(ind == T)
{

  if(count == 1)
  {
    node = node_temp[duplicated(node_temp) == F] # only select unique values #
    X = X_temp[duplicated(node_temp) == F]
    Y = Y_temp[duplicated(node_temp) == F]
    
    y_total[node+(n*(count-1))] = temp_Total[duplicated(node_temp) == F]
  }
  if (count != 1)
  {
    #browser()
    y_total = rbind(y_total, matrix(NA, nrow=n, ncol=dim(y_total)[2] )) ### ERROR
    #node = c(node,node_temp[duplicated(node_temp) == F]) # only select unique values #
    node = node_temp[duplicated(node_temp) == F]
    X = c(X,X_temp[duplicated(node_temp) == F])
    Y = c(Y,Y_temp[duplicated(node_temp) == F])

    y_total = cbind(y_total, rep(NA,(count)*n))
    y_total[node+(n*(count-1)),count] = temp_Total[duplicated(node_temp) == F]
  }
  node_temp = node_temp[duplicated(node_temp) == T] # choose only the nonunique values
  X_temp = X_temp[duplicated(node_temp) == T]
  Y_temp = Y_temp[duplicated(node_temp) == T]
  temp_Total = temp_Total[duplicated(node_temp) == T]
  
  ind = sum(duplicated(node_temp) == T) > 0
  count = count + 1
}

##### Count - 1 now denotes the number of replicated datasets ######

### Build the data set for use with INLA ##
count = count-1 # denotes number of replicates 

mu0 = c(rep(1,count*n))
ii = rep(c(1:n),count)
r = rep(1:count, each = n)

data = list(y=y_total,mu0=mu0,ii=ii,jj=jj,kk=kk)

####### Second method #######
formula = y ~ mu0 -1 + 
  f(ii, model = "matern2d", nrow=nrow, ncol=ncol, nu = 1,
    initial = c(2, log(20)), fixed=c(FALSE,FALSE),
    param = c(1,0.001, NA,NA), constr=TRUE, replicate=r)

res = inla(formula, family = c("gaussian","gaussian","gaussian"), data =data, verbose=TRUE,
           control.inla= list(strategy = "gaussian"))
summary(res)

res_summary = res$summary.random
res_means = res_summary$ii$mean + summary(res)$fixed[1,1]



plot.means1 = inla.vector2matrix(res_means[1:10000], nrow, ncol)
plot.means2 = inla.vector2matrix(res_means[10000+1:10000], nrow, ncol)

# Spatial plot means based on the actual replicated observed values.
image.plot(x=cutpointsx, y=cutpointsy, z=exp(plot.means2), xlab = "", ylab = "", axes = T, col=rev(heat.colors(10)))
points(x=tin$Xcoord, y=tin$Ycoord, pch = 15, cex=0.3)
title(main = "", font.main = 4)

threshold_field = res_means[10000+1:60000]
threshold_field[threshold_field<7 & threshold_field>5] = NA
min_field = rep(NA,n)
max_field=rep(NA,n)
med_field =rep(NA,n)
for (ind in 1:n)
{
    min_field[ind] = min(res_means[ind], threshold_field[10000+ind], threshold_field[20000+ind], threshold_field[30000+ind], threshold_field[40000+ind], threshold_field[50000+ind], na.rm=T)
    med_field[ind] = median(c(res_means[ind], threshold_field[10000+ind], threshold_field[20000+ind], threshold_field[30000+ind], threshold_field[40000+ind], threshold_field[50000+ind]), na.rm=T)
    max_field[ind] = max(res_means[ind], threshold_field[10000+ind], threshold_field[20000+ind], threshold_field[30000+ind], threshold_field[40000+ind], threshold_field[50000+ind], na.rm=T)
}

plot.min = inla.vector2matrix(min_field, nrow, ncol)
image.plot(x=cutpointsx, y=cutpointsy, z=plot.min, xlab = "", ylab = "", axes = T, col=rev(heat.colors(10)), zlim = c(0.1,9))
points(x=tin$Xcoord, y=tin$Ycoord, pch = 15, cex=0.3)
title(main = "", font.main = 4)

plot.med = inla.vector2matrix(med_field, nrow, ncol)
image.plot(x=cutpointsx, y=cutpointsy, z=plot.med, xlab = "", ylab = "", axes = T, col=rev(heat.colors(10)), zlim = c(0.1,9))
points(x=tin$Xcoord, y=tin$Ycoord, pch = 15, cex=0.3)
title(main = "", font.main = 4)

plot.max = inla.vector2matrix(max_field, nrow, ncol)
image.plot(x=cutpointsx, y=cutpointsy, z=plot.max, xlab = "", ylab = "", axes = T, col=rev(heat.colors(10)), zlim = c(0.1,9))
points(x=tin$Xcoord, y=tin$Ycoord, pch = 15, cex=0.3)
title(main = "", font.main = 4)

# Plot the lower 95 and upper 95 field #
res_LCL = res_summary$ii$`0.025quant` + summary(res)$fixed[1,5]
plot_LCL = inla.vector2matrix(res_LCL[1:10000], nrow, ncol)
res_UCL = res_summary$ii$`0.975quant` + summary(res)$fixed[1,3]
plot_UCL = inla.vector2matrix(res_UCL[1:10000], nrow, ncol)

# Spatial plot LCL 
image.plot(x=cutpointsx, y=cutpointsy, z=plot_LCL, xlab = "", ylab = "", axes = T, col=rev(heat.colors(10)),zlim=c(1.75,9))
points(x=tin$Xcoord, y=tin$Ycoord, pch = 15, cex=0.3)
title(main = "", font.main = 4)

# Spatial plot UCL 
image.plot(x=cutpointsx, y=cutpointsy, z=plot_UCL, xlab = "", ylab = "", axes = T, col=rev(heat.colors(10)),zlim=c(1.75,9))
points(x=tin$Xcoord, y=tin$Ycoord, pch = 15, cex=0.3)
title(main = "", font.main = 4)

