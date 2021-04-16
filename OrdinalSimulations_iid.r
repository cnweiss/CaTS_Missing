#Simulate ordinal time series with missing data
path <- "resOrdinalSimulations_iid.txt"

#For simuation, we use the rank count approach

IOV <- function(f){
	4/length(f) * sum(f*(1-f))
}

skew <- function(f){
	2/length(f) * sum(f) - 1
}



reps <- 1e4
write(reps, path, append=TRUE)

maxlag <- 5
alpha <- 0.05 #level for test

#Data examples imply: 3, 0.2, 0.35 and 10, 0.45, 0.5

parN <- 3 #10
parp <- 0.15 #0.2 #0.35 #0.45

bincodes <- diag(1, parN+1) #later used for binarization

#fraction observed:
tabpi <- c(1,0.95,0.9,0.85,0.8,0.75,0.5)

tabT <- c(50,100,250,500,1000)

for(parpi in tabpi){
for(T in tabT){
	
	results <- array(NA, c(reps, parN+2+2+maxlag))
	for(r in 1:reps){
		#Generate time series with i.i.d. rank counts
		data <- rbinom(T, parN,parp)
		
		#Consider missing data:
		obsnum <- T
		obsrate <- 1
		if(parpi<1){
			tabO <- rbinom(T, 1,parpi)
			data[tabO==0] <- NA
			obsnum <- sum(tabO)
			obsrate <- obsnum/T
		}
		
		#Binarization of data:
		databin <- bincodes[data+1,] #assign row vectors
		
		hatpi <- colMeans(databin, na.rm=TRUE)
		hatf <- cumsum(hatpi)[-(parN+1)]
		
		hatIOV <- IOV(hatf)
		hatskew <- skew(hatf)
		
		#Asymptotics with randomly missing data:
		hatfmat <- array(0, c(parN,parN))
		for(i in c(1:parN)){
		for(j in c(1:parN)){
			hatfmat[i,j] <- (hatf[min(i,j)]-hatf[i]*hatf[j])*(hatf[min(i,j)]+(1-2*obsrate)*hatf[i]*hatf[j])
		}} #for i,j
		
		crit2 <- qnorm(1-alpha/2, sd=sqrt(16/parN^2/hatIOV^2/T/obsrate^2 * sum(hatfmat)))
		#gives critical values
		#-1/obsnum -+ crit2


		hatbivprob <- array(0,c(parN,parN,maxlag))
		hatbivcdf <- array(0, c(parN,maxlag))
		cohenord <- rep(0,maxlag)

		for(k in c(1:maxlag)){
			for(i in c(1:parN)){
				for(j in c(1:parN)){
					hatbivprob[i,j,k] <- mean(databin[(k+1):T,i]*databin[1:(T-k),j], na.rm=TRUE)
				}
			#diagonal of bivariate cdf
			hatbivcdf[i,k] <- sum(hatbivprob[1:i,1:i,k])
			}

			#ordinal kappa:
			cohenord[k] <- (sum(hatbivcdf[,k])-sum(hatf^2))/sum(hatf*(1-hatf))
		} #for k
		
		results[r,] <- c(hatf, hatIOV,hatskew, -1/obsnum-crit2,-1/obsnum+crit2, cohenord)
	} #for reps
	
	res <- c(parN,parp,parpi,T)
	res <- c(res, colMeans(results[,1:parN])) #mean of hatf
	res <- c(res, cov(results[,1:parN])[upper.tri(matrix(0, parN, parN), diag=TRUE)]) #covariances
	res <- c(res, mean(results[,parN+1])) #mean IOV
	res <- c(res, var(results[,parN+1])) #var IOV
	res <- c(res, mean(results[,parN+2])) #mean skew
	res <- c(res, var(results[,parN+2])) #var skew
	res <- c(res, colMeans(results[,(parN+4+1):(parN+4+maxlag)])) #mean cohen
	res <- c(res, apply(results[,(parN+4+1):(parN+4+maxlag)], 2, var)) #var cohen
	for(k in 1:maxlag){
		res <- c(res, sum(1-(results[,(parN+4+k)]>results[,(parN+3)])*(results[,(parN+4+k)]<results[,(parN+4)]))) #rejections
	}
	write(res, path, ncolumns=length(res), append=TRUE)
} #for T
} #for pi



