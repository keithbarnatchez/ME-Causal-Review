#-------------------------------------------------------------------------------
### poss corr, positive z coeff
set.seed(536457) # 536457
n <- 25
corxz <- 0.5
Sigma <- matrix(c(1,corxz,corxz,1),nrow=2,byrow=T)
Xmat <- rmvnorm(n=n,mean = c(0,0),sigma = Sigma)
me_z <- 1
me_x <- .5
bx <- 2
sigy <- 1


# covariates
x <- Xmat[,1] ; z <- Xmat[,2]
x_meas <- x + me_x*rnorm(n)
z_meas <- z + me_z*rnorm(n)

# outcome
y <- rep(bx*x + 2*z + sigy*rnorm(n), 2)

# partial out z with and without measurement error
xres_true <- lm(x ~ z)$residuals
xres_err <- lm(x_meas ~ z_meas)$residuals

X <- c(xres_err,xres_true)
Measurement <- c(rep('With error',length(x)),
                 rep('True',length(x)))
pair <- rep(1:length(x),2)

df <- data.frame(y,X,Measurement,pair,corr='corr(X,Z)=0.5', zcoeff='+ z effect')
#-------------------------------------------------------------------------------
### negative corr, positive z coeff
n <- 25
corxz <- -0.5
Sigma <- matrix(c(1,corxz,corxz,1),nrow=2,byrow=T)
Xmat <- rmvnorm(n=n,mean = c(0,0),sigma = Sigma)

# covariates
x <- Xmat[,1] ; z <- Xmat[,2]
x_meas <- x + me_x*rnorm(n)
z_meas <- z + me_z*rnorm(n)

# outcome
y <- rep(bx*x + 2*z + sigy*rnorm(n), 2)

# partial out z with and without measurement error
xres_true <- lm(x ~ z)$residuals
xres_err <- lm(x_meas ~ z_meas)$residuals

X <- c(xres_err,xres_true)
Measurement <- c(rep('With error',length(x)),
                 rep('True',length(x)))
pair <- rep(1:length(x),2)

dftemp <- data.frame(y,X,Measurement,pair,corr='corr(X,Z)=-0.5', zcoeff='+ z effect')
df <- rbind(df,dftemp)
#-------------------------------------------------------------------------------
### negative corr, negative z coeff
n <- 25
corxz <- -0.5
Sigma <- matrix(c(1,corxz,corxz,1),nrow=2,byrow=T)
Xmat <- rmvnorm(n=n,mean = c(0,0),sigma = Sigma)

# covariates
x <- Xmat[,1] ; z <- Xmat[,2]
x_meas <- x + me_x*rnorm(n)
z_meas <- z + me_z*rnorm(n)

# outcome
y <- rep(bx*x - 2*z + sigy*rnorm(n), 2)

# partial out z with and without measurement error
xres_true <- lm(x ~ z)$residuals
xres_err <- lm(x_meas ~ z_meas)$residuals

X <- c(xres_err,xres_true)
Measurement <- c(rep('With error',length(x)),
                 rep('True',length(x)))
pair <- rep(1:length(x),2)

dftemp <- data.frame(y,X,Measurement,pair,corr='corr(X,Z)=-0.5', zcoeff='- z effect')
df <- rbind(df,dftemp)
#-------------------------------------------------------------------------------
### positive corr, negative z coeff
n <- 25
corxz <- 0.5
Sigma <- matrix(c(1,corxz,corxz,1),nrow=2,byrow=T)
Xmat <- rmvnorm(n=n,mean = c(0,0),sigma = Sigma)

# covariates
x <- Xmat[,1] ; z <- Xmat[,2]
x_meas <- x + me_x*rnorm(n)
z_meas <- z + me_z*rnorm(n)

# outcome
y <- rep(bx*x - 2*z + sigy*rnorm(n), 2)

# partial out z with and without measurement error
xres_true <- lm(x ~ z)$residuals
xres_err <- lm(x_meas ~ z_meas)$residuals

X <- c(xres_err,xres_true)
Measurement <- c(rep('With error',length(x)),
                 rep('True',length(x)))
pair <- rep(1:length(x),2)

dftemp <- data.frame(y,X,Measurement,pair,corr='corr(X,Z)=0.5', zcoeff='- z effect')
df <- rbind(df,dftemp)

#-------------------------------------------------------------------------------
# make final plot
pan_plot <- ggplot(df, aes(y=y,x=X,color=Measurement,group=pair)) + geom_point() + theme_bw() +
  geom_line(color='black',alpha=0.3) + facet_grid(factor(zcoeff)~factor(corr)) +
  geom_smooth(method='lm',aes(x=X,y=y,group=Measurement),se=F) +
  scale_color_manual( values = c("#00BFC4", "#F8766D")) +
  theme(legend.position="bottom") +
  labs(x = 'X - E[X|Z]', y = 'Y', title='Varying effects of measurement error') 
pan_plot
ggsave('../../output/figures/me-effects-grid.pdf',width = 6, height = 4,units = 'in')