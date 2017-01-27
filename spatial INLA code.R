###read data
lead97=read.table("lead97-new.txt")
names(lead97) <- c("X","Y","data")
lead00=read.table("lead00-new.txt")
names(lead00) <- c("X","Y","data")

###log-scale and scale coordinates
loglead97 <- data.frame(X=lead97$X/100000,Y=lead97$Y/100000,data=log(lead97$data))
loglead00 <- data.frame(X=lead00$X/100000,Y=lead00$Y/100000,data=log(lead00$data))

###define the dimension of the grid
nrow=100
ncol=100
n=nrow*ncol
xrange=range(c(loglead97$X,loglead00$X))
yrange=range(c(loglead97$Y,loglead00$Y))
cutpointsx=seq((xrange[1]-0.1),(xrange[2]+0.1),length.out=nrow)
cutpointsy=seq((yrange[1]-0.1),(yrange[2]+0.1),length.out=ncol)

###work on lead97
grid1.97 <- as.numeric(cut(loglead97$X,cutpointsx))
grid2.97 <- as.numeric(cut(loglead97$Y,cutpointsy))
grid97=cbind(grid1.97,grid2.97)
node97=numeric(dim(loglead97)[1])
for(i in 1: dim(loglead97)[1])
  node97[i]=inla.lattice2node(grid97[i,1],grid97[i,2],nrow=nrow,ncol=ncol)

###work on lead00
grid1.00 <- as.numeric(cut(loglead00$X,cutpointsx))
grid2.00 <- as.numeric(cut(loglead00$Y,cutpointsy))
grid00=cbind(grid1.00,grid2.00)
node00=numeric(dim(loglead00)[1])
for(i in 1: dim(loglead00)[1])
  node00[i]=inla.lattice2node(grid00[i,1],grid00[i,2],nrow=nrow,ncol=ncol)

### Build the data set for use with INLA ##
y97 = rep(NA,n)
y97[node97] =  loglead97$data
y00 = rep(NA,n)
y00[node00] =  loglead00$data

pois97 = rep(0,n)
pois97[node97] = 1

yy=matrix(NA,4*n,2)
yy[1:n,1] = y97
yy[n+1:n,1 ] = y00
yy[2*n+1:n,2]=pois97

mu0 = c(rep(1,n),rep(0,3*n))
mu1 = c(rep(0,n),rep(1,n),rep(0,2*n))
alpha = c(rep(0,2*n),rep(1,n),rep(0,n))
ii = c(1:n,1:n,rep(NA,2*n))
jj = c(rep(NA,2*n),1:n,1:n)
replicates = c(rep(1,n),rep(2,n),rep(1,n),rep(2,n))

data = list(yy=yy,mu0=mu0,mu1=mu1,alpha=alpha,ii=ii,jj=jj,replicates=replicates)

### define the model and run it
formula = yy ~ alpha + mu0 + mu1 -1 + 
  f(ii, model = "matern2d", nrow=nrow, ncol=ncol, replicate=replicates, nu = 1,
    initial = c(2, log(20)), fixed=c(FALSE,FALSE),
    param = c(1,0.001, NA,NA), constr=TRUE) +
  f(jj, copy="ii",replicate=replicates, fixed=FALSE, param=c(0,0.1), initial=0.1)

res = inla(formula, family = c("gaussian", "poisson"), data =data, verbose=TRUE,
           control.inla= list(strategy = "gaussian"))
summary(res)

### Uncheck to obtain a better estimate of the model hyperparameters
#h.res=inla.hyperpar(res,dz=1,diff.logdens=15,verbose=T)
#summary(h.res)

# repeat just for the 1997 data #
yy = matrix(NA, 2*n, 2)
yy[1:n,1] = y97
yy[n+(1:n),2] = pois97
mu_97 = c(rep(1,times = n), rep(0, times = n))
alpha_97 = c(rep(0, times =  n), rep(1, times = n))
ii = c(1:n, rep(NA, n))
jj = c(rep(NA,n), 1:n)
data_97 = list(yy=yy,mu_97=mu_97,alpha_97=alpha_97,ii=ii,jj=jj)

formula_97 = yy ~ alpha_97 + mu_97 -1 + 
  f(ii, model = "matern2d", nrow=nrow, ncol=ncol, nu = 1,
    initial = c(2, log(20)), fixed=c(FALSE,FALSE),
    param = c(1,0.001, NA,NA), constr=TRUE) +
  f(jj, copy="ii", fixed=FALSE, param=c(0,0.1), initial=0.1)

res_97 = inla(formula_97, family = c("gaussian", "poisson"), data =data_97, verbose=TRUE,
              control.inla= list(strategy = "gaussian"))

summary(res_97)
### This is to obtain a better estimate of the model hyperparameters
#h.res_97 = inla.hyperpar(res_97,dz=1,diff.logdens=15,verbose=T)
#summary(h.res_97)

# repeat just for the 2000 data #
yy = y00
mu_00 = c(rep(1,times = n))
ii = c(1:n)
data_00 = list(yy=yy,mu_00=mu_00,ii=ii)

formula_00 = yy ~ mu_00 -1 + 
  f(ii, model = "matern2d", nrow=nrow, ncol=ncol, nu = 1,
    initial = c(2, log(20)), fixed=c(FALSE,FALSE),
    param = c(1,0.001, NA,NA), constr=TRUE)

res_00 = inla(formula_00, family = c("gaussian") , data =data_00, verbose=TRUE,
              control.inla= list(strategy = "gaussian"))

summary(res_00)
### This is to obtain a better estimate of the model hyperparameters
#h.res_00 = inla.hyperpar(res_00,dz=1,diff.logdens=15,verbose=T)
#sumary(h.res_00)

# plot an image plot of the mean S surface from the joint model #
# Note since both 2000 and 1997 shared the same latent process S
# need only extract one #

pref.summary_whole = h.res$summary.random
pref.means_1997 = pref.summary_whole$ii$mean[1:(nrow*ncol)] + summary(h.res)$fixed[2,1]
plot.pref.means_1997 = inla.vector2matrix(pref.means_1997, nrow, ncol)
pref.means_2000 = pref.summary_whole$ii$mean[nrow+(1:(nrow*ncol))] + summary(h.res)$fixed[3,1]
plot.pref.means_2000 = inla.vector2matrix(pref.means_2000, nrow, ncol)

dev.off()
par(pty="s",mfrow=c(1,1))
image(plot.pref.means_1997)
title("Preferentially adjusted fitted surace")

galicia.poly2 = galicia.poly
galicia.poly2[,1] = 100* (galicia.poly2[,1] - min(galicia.poly2[,1])) / (max(galicia.poly2[,1]) - min(galicia.poly2[,1]))
galicia.poly2[,2] = 100* (galicia.poly2[,2] - min(galicia.poly2[,2])) / (max(galicia.poly2[,2]) - min(galicia.poly2[,2]))

plot_data_97 = data.frame(x_97 = rep(c(1:nrow),ncol), y_97 = rep(c(1:ncol), each = nrow),
                          z_97 = pref.means_1997)
heat.plot_97 = ggplot(data = plot_data_97) + geom_raster(aes(x=x_97, y=y_97, fill=z_97, colour = "grey50"), interpolate = T) +
  geom_polygon(data = data.frame(poly_x = galicia.poly2[,1], poly_y = galicia.poly2[,2]), aes(x = poly_x, y = poly_y), colour="black", fill=NA) +
  xlab('') + ylab('') + ggtitle('Predicted log lead distribution') +
  theme(legend.position = "none", axis.text = element_text(size = 0))
heat.plot_97
