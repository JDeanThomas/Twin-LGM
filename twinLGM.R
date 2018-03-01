
require(OpenMx)
require(MASS)
require(psych)
require(stringr)
require(stats)
require(polycor)
setwd("~/Documents/Classes/Spring 2013/Stats for Genetics II/project/updates LGM")

#-------------------------------------------------------------------

dat <- read.table("modelvars_contin_fullresponses_wide.csv",sep=",",header=TRUE,na.strings="-999")
head(dat)
summary(dat)

n <- names(dat)
n <- str_replace(n,"[.]","_")
names(dat) <- n
head(dat)


# Remove observations with NA vals for definition variables
model <- dat[is.na(dat$pnanx_1)==FALSE & is.na(dat$pnanx_2)==FALSE,]
summary(model)



#-------------------------------------------------------------------

	# Select smoking frequency variables
smk <- model[,c(2,6:9,36:41)] 			# peer

head(smk)
names(smk) <- c("zyg","sex_1","sex_2","anx_1","anx_2","smk14_1","smk14_2","smk17_1","smk17_2","smk22_1","smk22_2")
summary(smk)

smkMZ <- subset(smk,zyg==1)
smkMZ2 <- smkMZ[-1]
smkDZ <- subset(smk,zyg==2)
smkDZ2 <- smkDZ[-1]

#-------------------------------------------------------------------


selvars <- c("smk14_1","smk17_1","smk22_1","smk14_2","smk17_2","smk22_2")

# 3 Latent factors - w/quadratic

nf <- 3		# Number of latent factors (intercept+slope)



# Labels for a, c & e pathways from A, C & E factors to INTECEPT & SLOPE factors
AlLabs   <- paste("al", do.call("c", sapply(seq(1, nf), function(x){ paste(x:nf, x,sep="") })), sep="")							
ClLabs   <- paste("cl", do.call("c", sapply(seq(1, nf), function(x){ paste(x:nf, x,sep="") })), sep="")
ElLabs   <- paste("el", do.call("c", sapply(seq(1, nf), function(x){ paste(x:nf, x,sep="") })), sep="")

# Labels for LGC latent factor loadings. Matrix labels are populated down the column
FlLabs   <- matrix(NA,nv,nf)
for (i in 1:nv) {
	FlLabs[i,1:nf] <- paste("f",1:nf,"_v",i,sep="")
}


# values for factor (slope/intercept) loadings
facVals <- c(1,1,1,
		   0,3.5,8,
		   0,12.25,64)


# Matrices ac, cc, and ec to store a, c, and e path coefficients from latent factors(s) to Int & Slope 
pathAl   <- mxMatrix( type="Lower", nrow=nf, ncol=nf, free=c(T,T,T,T,T,T), values=c(0,1,.5,0,0,0), labels=AlLabs, name="al" )
pathCl   <- mxMatrix( type="Lower", nrow=nf, ncol=nf, free=c(T,T,T,T,T,T), values=c(1,1,.5,0,0,0), labels=ClLabs, name="cl" )
pathEl   <- mxMatrix( type="Lower", nrow=nf, ncol=nf, free=c(T,T,T,T,T,T), values=c(0,1,1,0,0,0), labels=ElLabs, name="el" )
# NB: No constraint on the lower & upper bounds

# Matrix f for fixed factor loadings from Intercept & Slope to observed variables
pathFl    <- mxMatrix( type="Full", nrow=nv, ncol=nf, free=FALSE, values=facVals, labels=FlLabs,name="fl" )

# Matrices as, cs, and es to store a, c, and e path coefficients for specific factors (residuals)
pathAs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=c(1,2,2), labels=AsLabs, name="as" )
pathCs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=T, values=0, labels=CsLabs, name="cs" )
pathEs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=c(1,4,4), labels=EsLabs, name="es" )


# Matrices generated to hold A, C, and E computed Variance Components
covA      <- mxAlgebra( expression=fl %&% (al %*% t(al)) + as %*% t(as), name="A" )
covC      <- mxAlgebra( expression=fl %&% (cl %*% t(cl)) + cs %*% t(cs), name="C" )
covE      <- mxAlgebra( expression=fl %&% (el %*% t(el)) + es %*% t(es), name="E" )

# Algebra to compute total variances and standard deviations (diagonal only)
covP      <- mxAlgebra( expression=A+C+E, name="V" )
matI      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
invSD     <- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD")
rowvars  <- rep('vars',nv)
colvars  <- rep(c('A','C','E','SA','SC','SE'),each=nv)
estvars  <- mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V), name="vars", dimnames=list(rowvars,colvars))

# Matrices to compute A, C & E variance components for INTERCEPT & SLOPE & Quadratic
covAis      <- mxAlgebra( expression=(al %*% t(al)), name="Ais" )
covCis      <- mxAlgebra( expression=(cl %*% t(cl)), name="Cis" )
covEis      <- mxAlgebra( expression=(el %*% t(el)), name="Eis" )
covPis      <- mxAlgebra( expression=Ais+Cis+Eis, name="Vis" )		
matIis      <- mxMatrix( type="Iden", nrow=nf, ncol=nf, name="Iis")
lat_var_SD       <- mxAlgebra( expression=solve(sqrt(Iis*Vis)), name="iSDis")
# Algebras generated to hold Parameter Estimates and Derived Variance Components for INTERCEPT & SLOPE
 rowvars2   <- rep('latentvars',nf)
 colvars2   <- rep(c('Ais','Cis','Eis','SAis','SCis','SEis'),each=nf)
 estvars2   <- mxAlgebra( expression=cbind(Ais,Cis,Eis,Ais/Vis,Cis/Vis,Eis/Vis), name="latentvars", dimnames=list(rowvars2,colvars2))



					
# Specify intercept & slope latent factor means = triangles coming off latent INT & SLOPE
MeansIS   <- mxMatrix( type="Full", nrow=nf, ncol=1, free=T, labels=c("Im","Sm","Qm"), values=c(1,4,1), name="LMeans" )

pathB1	  <- mxMatrix( type="Full", nrow=nf, ncol=1, free=T, values=c(0,-.5,-1), labels=c("B_int_sex","B_slo_sex","B_quad_sex"), name="Beta1" ) 

# Matrix for calling up definition variables
sex1   <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.sex_1"), name="sex1")
sex2   <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.sex_2"), name="sex2")

pathB2	  <- mxMatrix( type="Full", nrow=nf, ncol=1, free=T, values=c(-.02), labels=c("B_int_anx","B_slo_anx","B_quad_anx"), name="Beta2" ) 

# Matrix for calling up definition variables
anx1   <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.anx_1"), name="anx1")
anx2   <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.anx_2"), name="anx2")



# Calcualate expected means + sex effects
Means1 <- mxAlgebra( expression= ( t((fl %*% ( LMeans + (sex1 %x% Beta1) + (anx1 %x% Beta2)))) ), name="Mean1")   
Means2 <- mxAlgebra( expression= ( t((fl %*% ( LMeans + (sex2 %x% Beta1) + (anx2 %x% Beta2)))) ), name="Mean2")   
eMeans <- mxAlgebra( expression= cbind( Mean1,Mean2 ), name="Means" )	



# Expected Variance/Covariance Matrices
# NB LABEL ALL A, C & E PATHWAYS. Otherwise, MX will duplicate estimates & ~double # of parameters                                           
covMZ     <- mxAlgebra( expression= rbind( cbind( A+C+E , A+C),
                                           cbind( A+C   , A+C+E)),    name="expCovMZ" )                                                                                  
covDZ     <- mxAlgebra( expression= rbind( cbind( A+C+E    , 0.5%x%A+C),
                                           cbind(0.5%x%A+C , A+C+E)), name="expCovDZ" )

# Data objects for Multiple Groups
dataMZ   <- mxData( observed=smkMZ2, type="raw" )
dataDZ   <- mxData( observed=smkDZ2, type="raw" )

# Objective objects for Multiple Groups: Likelihood of each individual's data given the model summed
# or aggregated across the entire sample.
objMZb    <- mxFIMLObjective( covariance="expCovMZ", means="Means", dimnames=selvars)
objDZb    <- mxFIMLObjective( covariance="expCovDZ", means="Means", dimnames=selvars)

# Combine Groups				
pars	 <- c( pathAl, pathCl, pathEl, pathFl, pathAs, pathCs, pathEs, covA, covC, covE, covP, matI, invSD, estvars, MeansIS,covAis, covCis, covEis, covPis, matIis, lat_var_SD, estvars2)		
defs	 <- c(sex1, sex2, anx1, anx2, pathB1, pathB2, Means1, Means2, eMeans)		
modelMZ  <- mxModel( pars, defs, covMZ, dataMZ, objMZb, name="MZ" ) 
modelDZ  <- mxModel( pars, defs, covDZ, dataDZ, objDZb, name="DZ" ) 
minus2ll <- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj      <- mxAlgebraObjective( "m2LL" )
LgcAceModelQ <- mxModel( "LgcACE", pars, modelMZ, modelDZ, minus2ll, obj )

LgcAceFitQ <- mxRun(LgcAceModelQ)
LgcAceFitQ <- mxRun(LgcAceFitQ)
LgcAceFitQ <- mxRun(LgcAceFitQ)
LgcAceFitQ <- mxRun(LgcAceFitQ)
LgcAceFitQ <- mxRun(LgcAceFitQ)
LgcAceFitQ <- mxRun(LgcAceFitQ)
LgcAceSummQ <- summary(LgcAceFitQ)
LgcAceSummQ

