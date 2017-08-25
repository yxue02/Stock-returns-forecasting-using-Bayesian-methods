library('mnormt');
library('MCMCpack');
rm(list=ls())
library("gplots");
setwd("D:\\project");
data = read.csv('sp100monthly.csv');

y = data.matrix(data[,4:ncol(data)]);
d = ncol(y);
y.pre=matrix(NaN,d,1) #Predict y
sigma.pre=matrix(NaN,d,1) #Predict sigma
bound.pre=matrix(NaN,d,2) #95% credible interval
T=1100;
month.last=60;
ymatrix.pre=matrix(NaN,d,nrow(y)-month.last);
sigmamatrix.pre=matrix(NaN,d,nrow(y)-month.last);
coverage=rep(NaN,nrow(y)-month.last);
mse=rep(NaN,nrow(y)-month.last);
q1matrix=matrix(NaN,d,nrow(y)-month.last);
q2matrix=matrix(NaN,d,nrow(y)-month.last);

for (mf in 1: (nrow(y)-month.last)){
y.true=y[month.last+mf,1];
	for (j in 1:d){
		y1=y[1:month.last,j];
		N = length(y1);
		#mu1, mu2 and phi1, phi2 are normal gamma distributed.
		mu0=0;
		kappa0=10;
		nu0=10;
		SS0=1;
		SS=sum((y1-mean(y1))^2);
		mu=matrix(NaN,T,2);
		phi=matrix(NaN,T,2);
		h.seq=matrix(NaN,T,1);
		phi[1,]=rgamma(2,nu0/2,SS0/2);
		mu[1,]=rnorm(2,0,sqrt(1/kappa0/phi[1,]));

		muphi=rbind(mu[1,],phi[1,]);
		muphi.new=muphi[,order(muphi[1,])];
		mu[1,]=muphi.new[1,];
		phi[1,]=muphi.new[2,];
		
			

		# q=[q1,q2] where q1~Beta(a1,b1) is the probability of 1 to 1. q2~Beta(a2,b2) is the probability of 2 to 1.
		a1=1;
		b1=1;
		a2=1;
		b2=1
		q=matrix(NaN,T,2);
		q[1,]=c(0.5,0.5);
		pi.l=matrix(NaN,N,2);
		pi=matrix(NaN,N,2);
		P.h=matrix(NaN,N,2);
		P.h1=matrix(NaN,N,2);
		P.l=matrix(NaN,N,4); #1 to 1, 1 to 2, 2 to 1, 2 to 2
		SS.hat=matrix(NaN,1,2);
		mu.hat=matrix(NaN,1,2);
		n.seq=matrix(NaN,T,2);
		#nt is the number of transitions
		nt=rbind(c(0,0),c(0,0));

		h= sample(c(1,2),N, replace = TRUE, prob = c(0.5,0.5));
		for (i in 2:T){

			pi.l[1,]=c(dnorm(y1[1],mu[i-1,1],sqrt(1/phi[i-1,1])),dnorm(y1[1],mu[i-1,2],sqrt(1/phi[i-1,2])));
			pi[1,]=pi.l[1,]/sum(pi.l[1,]);
			for (t in 2:N){
				for (k in 1:4){
					index1=ceiling(k/2);
					index2=k%%2;
					P.l[t,k]=pi[t-1,index1]*(q[i-1,index1])^index2*(1-q[i-1,index1])^(1-index2)*dnorm(y1[t],mu[i-1,index1],sqrt(1/phi[i-1,index1]));
					}
				pi.l[t,]=c(P.l[t,1]+P.l[t,3],P.l[t,2]+P.l[t,4]);
				pi[t,]=pi.l[t,]/sum(pi.l[t,]);
				}
			P.h[N,]=pi[N,]
			for (t in N:2){
				h[t]=sample(c(1,2),1,,prob=P.h[t,]);
				P.h1[t-1,]=c(P.l[t,h[t]],P.l[t,2+h[t]]);
				P.h[t-1,]=P.h1[t-1,]/sum(P.h1[t-1,]);
				}
			h[1]=sample(c(1,2),1,,prob=P.h[1,]);

			h.seq[i]=h[N];
			n=c(sum(h==1),sum(h==2));
			n.seq[i,]=n;
			nt=rbind(c(0,0),c(0,0));
			for (t in 2:N){
				nt=nt+rbind(c((h[t-1]==1)*(h[t]==1),(h[t-1]==1)*(h[t]==2)),c((h[t-1]==2)*(h[t]==1),(h[t-1]==2)*(h[t]==2)));
				}
			kappa.hat=kappa0+n;
			nu.hat=nu0+n;
			for (k in 1:2){
				if (n[k]!=0){
					SS.hat[k]=sum((y1[h==k]-mean(y1[h==k]))^2)+n[k]*kappa0/kappa.hat[k]*(mean(y1[h==k])-mu0)^2+SS0+SS;
					mu.hat[k]=(sum(y1[h==k])+kappa0*mu0)/kappa.hat[k];
				} else {
					SS.hat[k]=SS0+SS;;
					mu.hat[k]=mu0;
				}
			}
			phi[i,]=rgamma(2,nu.hat/2,SS.hat/2);
			mu[i,]=rnorm(2,mu.hat,sqrt(1/kappa.hat/phi[i,]));
			a1.hat=nt[1,1]+a1;
			b1.hat=nt[1,2]+b1;
			a2.hat=nt[2,1]+a2;
			b2.hat=nt[2,2]+b2;
			q[i,]=c(rbeta(1,a1.hat,b1.hat),rbeta(1,a2.hat,b2.hat));
			#When there is label switching q needs changes
			if (mu[i,1]>mu[i,2]){
				q[i,]=1-q[i,];
			}
			muphinq=rbind(mu[i,],phi[i,],n,q[i,]);	
			muphinq.new=muphinq[,order(muphinq[1,])];
			mu[i,]=muphinq.new[1,];
			phi[i,]=muphinq.new[2,];
			n=muphinq.new[3,];
			q[i,]=muphinq.new[4,];
		}
		burnin=100;
		mu.remain=mu[(burnin+1):T,];
		phi.remain=phi[(burnin+1):T,];
		sigma.remain=sqrt(1/phi.remain);
		h.remain=h.seq[(burnin+1):T,];
		q.remain=q[(burnin+1):T,];
		q1matrix[j,mf]=mean(q.remain[,1]);
		q2matrix[j,mf]=mean(q.remain[,2]);
		mu.plot=mu.remain[1:1000,];
		phi.plot=phi.remain[1:1000,];
		sigma.plot=sqrt(1/phi.plot);	
		plot(seq(1,1000,1),mu.plot[,1],type="l",col=1,ylim=c(min(mu.plot),max(mu.plot)));
		lines(seq(1,1000,1),mu.plot[,2],lty=3,col=2);
		#plot(seq(1,1000,1),sigma.plot[,1],type="l",col=1,ylim=c(min(sigma.plot),max(sigma.plot)));
		#lines(seq(1,1000,1),sigma.plot[,2],lty=3,col=2);
		#Predict next period;
		remain=T-burnin;
		y1.pre=matrix(NaN,remain,1);
		h.pre=matrix(NaN,remain,1);
			for (i in 1:remain){
			p.pre=c((h.remain[i]==1)*q.remain[i,1]+(h.remain[i]==2)*q.remain[i,2],(h.remain[i]==1)*(1-q.remain[i,1])+(h.remain[i]==2)*(1-q.remain[i,2]));
			h.pre[i]=sample(c(1,2),1,,prob=p.pre);
			y1.pre[i]=sum(c(h.pre[i]==1,h.pre[i]==2)*rnorm(2,mu.remain[i,],sigma.remain[i,]));
			}
		bound.pre[j,]=c(quantile(y1.pre,0.025),quantile(y1.pre,0.975));
		y.pre[j]=mean(y1.pre);
		sigma.pre[j]=sd(c(y1.pre));	
		if (j%%5==0){print(sprintf("mf = %d, j=%d.", mf,j))
		}
	}
	coverage[mf]=mean((bound.pre[,1]<y.true)*(bound.pre[,2]>y.true));
	mse[mf]=sd(c(y.pre-y.true));
	ymatrix.pre[,mf]=y.pre;
	sigmamatrix.pre[,mf]=sigma.pre;	
}
result=matrix(NaN,(4*d+2),nrow(y)-month.last);
result[1:d,]=t(ymatrix.pre);
result[(d+1):(2*d),]=t(sigmamatrix.pre);
result[(2*d+1):(3*d),]=t(q1matrix);
result[(3*d+1):(4*d),]=t(q2matrix);
result[(4*d+1),]=coverage;
result[(4*d+2),]=mse;
tresult=t(result);
write.csv(result, file = "sp100_FB", row.names = FALSE)
hist(q1matrix);
dev.copy2eps(file="q1FB.eps")
hist(q2matrix);
dev.copy2eps(file="q2FB.eps")
print(sprintf("The mean q1 is %.2f and the mean q2 is %.2f using FB method", mean(q1matrix), mean(q2matrix)));