#Analysis of example series

path <- "xxx" #path for data


#Functions for computing marginal properties:

#Compute median from ECDFs (states 0-m):
f2med <- function(f){
	m <- length(f)
	if(max(is.na(f))==1){
		return(NA)
	}else if(f[m]<0.5){
		return(m)
	}else{
		return(min( (0:(m-1))[f>=0.5] ))
	}
}

#Compute mode from PMFs (states 0-m):
p2mod <- function(p){
	m <- length(p)-1
	if(max(is.na(p))==1){
		return(NA)
	}else{
		return((0:m)[which.max(p)])
	}
}

IOV <- function(f){
	4/length(f) * sum(f*(1-f))
}

skew <- function(f){
	2/length(f) * sum(f) - 1
}

maxlag <- 10
alpha <- 0.05
z2 <- qnorm(1-alpha/2)







######################
#Actual data analysis:
######################


load(paste(path, "DATA.rda", sep=""))

#For peak severity data, DATA has
#Levels: none mild moderate severe
#Then, transformation to rank counts via
data <- as.numeric(DATA)-1

#For stress data, DATA already numerically coded (Likert 0-10), so
data <- DATA

T <- length(data)


#Missing data:
Odata <- 1-is.na(data)
T-sum(Odata) #number and ...
(1:T)[Odata==0] #position of missing data
hatpi <- sum(Odata)/T #estimated observation rate

acf(Odata)


#Ordinal scale:

#For peak severity data:
states <- levels(DATA)
states #"none"     "mild"     "moderate" "severe"

#For stress data:
states <- 0:10 #Likert scale

nostates <- length(states)

plot(data, type="b", pch=19, cex=0.5, ylim=c(0,nostates-1), xlab = "t", ylab = "DATA", yaxt="n")
axis(side=2, at=c(0:(nostates-1)), labels=states, las=2)
abline(v=(1:T)[Odata==0], lty=3, lwd=1.5)




#Binarization:
bincodes <- diag(1,nostates)
databin <- bincodes[data+1,]

#Marginal properties (na.rm=TRUE for amplitude modulation):
hatp <- colMeans(databin, na.rm=TRUE) #PMF
hatf <- cumsum(hatp)[-nostates] #CDF

f2med(hatf) #Median
p2mod(hatp) #Mode
IOV(hatf)
skew(hatf)

#Plot of PMF:
plot(0:(nostates-1), hatp, type="h", lwd=3, xlab="Range", ylab="PMF", xaxt='n')
axis(side=1, at=c(0:(nostates-1)), labels=states)

#Plot of CDF:
plot(0:(nostates-1), c(hatf,1), type="h", lwd=3, xlab="Range", ylab="CDF", xaxt='n', ylim=c(0,1))
axis(side=1, at=c(0:(nostates-1)), labels=states)


#Computation of Cohen's kappa:

#Entries of enumerator of kappa's variance
hatfmat <- array(0, c(nostates-1,nostates-1))
for(i in c(1:(nostates-1))){
for(j in c(1:(nostates-1))){
	hatfmat[i,j] <- (hatf[min(i,j)]-hatf[i]*hatf[j])*(hatf[min(i,j)]+(1-2*hatpi)*hatf[i]*hatf[j])
}} #for i,j

#Asymptotic standard deviation (plug-in approach):
sda <- sqrt( 1/T/hatpi^2 * sum(hatfmat) / sum(hatf*(1-hatf))^2 )
crit2 <- z2*sda

hatbivprob <- array(0,c(nostates-1,nostates-1,maxlag))
hatbivcdf <- array(0, c(nostates-1,maxlag))
cohenord <- rep(0,maxlag)

for(k in c(1:maxlag)){
	for(i in c(1:(nostates-1))){
		for(j in c(1:(nostates-1))){
			hatbivprob[i,j,k] <- mean(databin[(k+1):T,i]*databin[1:(T-k),j], na.rm=TRUE)
			#if one of the databin =NA, then product =NA,
			#so na.rm=TRUE implements amplitude modulation approach
		}
	#diagonal of bivariate cdf
	hatbivcdf[i,k] <- sum(hatbivprob[1:i,1:i,k])
	}

	#ordinal kappa:
	cohenord[k] <- (sum(hatbivcdf[,k])-sum(hatf^2))/sum(hatf*(1-hatf))
}

#Plot of ordinal Cohen's kappa with randomly missing data:
low <- min(-1/T/hatpi-crit2, cohenord)
plot(cohenord, type="h", xlab = "h", ylab = expression(paste("Cohen's   ",kappa[ord](h))), lwd=2, ylim=c(low,1))
abline(h=c(-1/T/hatpi+crit2,-1/T/hatpi-crit2), col=gray(0.5), lwd=1, lty=2)
abline(h=-1/T/hatpi)



