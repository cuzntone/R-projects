library(MASS)
###Simulate data
n <- 500
set.seed(2345678)
genExp<-x <- matrix(rnorm(n*20, mean = 1, sd = sqrt(25)), nrow=n, byrow=T)
y <- x[,1] + 2*x[,4] + 1.5*x[,6] + 3*x[,8] + 0.5*x[,9] + rnorm(n, mean = 0, sd = sqrt(0.25))

###################
c<-100
nloops<-20000
burnIn<-nloops/2 #use the first nloops/2 as burn-in
k<-1
set.seed(seed=1000)

# assume the intercept is always included
Gam<-rbinom(ncol(genExp),1,prob=0.5)

matGam<-matrix(NA,nrow=nloops,ncol=ncol(genExp))

###Main function
Gibbs<-function(Gam,genExp,y,i)
{
  int<-rep(1,nrow(genExp))
#	xGam<-genExp[,which(Gam==1)]
	xGam<-cbind(int,genExp[,which(Gam==1)])
	
	sGam<-t(y)%*%y-c/(1+c)*t(t(xGam)%*%y)%*%ginv((t(as.matrix(xGam))
		%*%as.matrix(xGam)))%*%(t(xGam)%*%y)

	GamTemp<-Gam

	if (Gam[i]==0)
	{
		p0<-(-sum(Gam)/2)*log(1+c)-nrow(genExp)/2*log(sGam)

		GamTemp[i]<-1
		xGam<-cbind(int,genExp[,which(GamTemp==1)])
#		xGam<-genExp[,which(GamTemp==1)]

		sGam<-t(y)%*%y-c/(1+c)*t(t(xGam)%*%y)%*%ginv((t(as.matrix(xGam))%*%
			as.matrix(xGam)))%*%(t(xGam)%*%y)
		p1<-(-sum(GamTemp)/2)*log(1+c)-nrow(genExp)/2*log(sGam)

	}
	if (Gam[i]==1)
	{
		p1<-(-sum(Gam)/2)*log(1+c)-nrow(genExp)/2*log(sGam)

		GamTemp[i]<-0
		xGam<-cbind(int,genExp[,which(GamTemp==1)])
#		xGam<-genExp[,which(GamTemp==1)]

		sGam<-t(y)%*%y-c/(1+c)*t(t(xGam)%*%y)%*%ginv((t(as.matrix(xGam))%*%
			as.matrix(xGam)))%*%(t(xGam)%*%y)
		p0<-(-sum(GamTemp)/2)*log(1+c)-nrow(genExp)/2*log(sGam)

	}		
		p<--log(exp(p0-p1)+1)

		judge<-runif(1)
		if (p>log(judge)) Gam[i]<-1
		else Gam[i]<-0
		return(Gam)
}
				

for (j in 1:nloops)
{
	cat("loop ",j,"\n")
	for (i in 1:ncol(genExp))
	{
		Gam<-Gibbs(Gam,genExp,y,i)
	}

	if (j>burnIn) 
	{
	matGam[k,]<-Gam 
	int<-rep(1,nrow(genExp))
	k<-k+1	
	}
	
}

k0<-k-1
matGamTemp<-matGam
prob<-0
selectModel<-matrix(rep(0,k0*(ncol(genExp)+1)),ncol=ncol(genExp)+1)
modelInd<-0
for (j in 1:k0)
{
	if (sum(matGamTemp[j,])!=9*ncol(genExp))
	{
		cat("j ",j,"\n")
		cat("to catch ",matGamTemp[j,],"\n")
		selectModel[modelInd,1:ncol(genExp)]<-matGamTemp[j,]
		hit<-1
		if (j<k0)
		{
			for (jj in (j+1):k0)
			{
				if ((sum(abs(matGamTemp[j,]-matGamTemp[jj,]))==0) && 
					sum(matGamTemp[jj,])!=9*ncol(genExp)) 
				{
					hit<-hit+1
					matGamTemp[jj,]<-rep(9,ncol(genExp))

				}
			}			
		}	
		modelInd<-modelInd+1
	
		prob<-prob+hit/k0
		selectModel[modelInd,ncol(genExp)+1]<-hit/k0
	}
	else
	{
		cat("j ",j," hit >1 \n")
	}
	if (prob==1) break

}
sortedSelect<-selectModel[order(as.matrix(selectModel[,ncol(genExp)+1]),decreasing=T),]
sortedSelect[1:10,]

sigma2<-rep(NA, nloops/2)
beta<-matrix(rep(NA, 10*nloops), nrow=nloops/2)
for (k in 1:10000)
{
  xGam<-cbind(int,genExp[,which(best==1)])
  sGam<-t(y)%*%y-c/(1+c)*t(t(xGam)%*%y)%*%ginv((t(as.matrix(xGam))%*%as.matrix(xGam)))%*%(t(xGam)%*%y)
  sigma2[k] <- rinvgamma(1, shape=nrow(genExp)/2, scale = sGam/2)
  beta <- rmvnorm(n=1000,mean=(c/(1+c))*ginv((t(as.matrix(xGam))%*%as.matrix(xGam)))%*%(t(xGam)%*%y), 
          sigma=(c/(1+c))*ginv((t(as.matrix(xGam))%*%as.matrix(xGam)))*sigma2[k])
}
quantile(sigma2, c(0.025, 0.5, 0.975))
> plot(sigma2, type="l")
> title("History plot of sigma ")

try<-c(500000,100000,80000, 50000, 10000, 5000, 1000, 500,100,50,10,5,1,0.5,0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001)
wholeprod=matGam[1:10000,]%*% try
unique(wholeprod)[75]
prod1<-rep(NA, 964)
for (gg in 1: 964)
{
  prod1[gg]<-length(which(wholeprod==unique(wholeprod)[gg]))
}

