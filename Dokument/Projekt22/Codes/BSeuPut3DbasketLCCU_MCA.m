function [UAC,UA]=BSeuPut3DbasketLCCU_MCA(param)
N=1e7; 
n=3;
sigma=0.2;
rho=0.25;
V=sigma^2*(rho*ones(n,n)+(1-rho)*eye(n));
A=sqrtm(V);
r=0.06;
T=1;
K=40;
S0=K*[ones(1,n)]/n;
tic;
X=A*randn(n,N);
Y1=max(K-S0*exp(r*T-T*repmat(diag(V),1,N)/2+sqrt(T)*X),0)*exp(-r*T);
Y2=max(K-S0*exp(r*T-T*repmat(diag(V),1,N)/2-sqrt(T)*X),0)*exp(-r*T);
Y1G=max(K-K*exp(r*T-T*mean(diag(V))/2+sqrt(T)*mean(X)),0)*exp(-r*T);
Y2G=max(K-K*exp(r*T-T*mean(diag(V))/2-sqrt(T)*mean(X)),0)*exp(-r*T);
[~,mPG]=BS(K,K,T,r,r-sigma^2/2,sigma*sqrt(rho*(1-1/n)+1/n));
Sxy=sum((Y1+Y2)/2.*((Y1G+Y2G)/2-mPG));
Sxx=sum(((Y1G+Y2G)/2-mPG).^2);
beta=Sxy/Sxx;
e0=std((Y1+Y2)/2-beta*((Y1G+Y2G)/2-mPG))/sqrt(N);
dt=toc;
nrep=ceil((2*e0/param)^2);
if nrep>1
 disp(['Estimated running time: '  num2str(dt*(nrep-1)) ' s.']);   
end
PA(1)=mean((Y1+Y2)/2);
PG(1)=mean((Y1G+Y2G)/2);
for k=2:nrep
 X=A*randn(n,N);
 Y1=max(K-S0*exp(r*T-T*repmat(diag(V),1,N)/2+sqrt(T)*X),0)*exp(-r*T);
 Y2=max(K-S0*exp(r*T-T*repmat(diag(V),1,N)/2-sqrt(T)*X),0)*exp(-r*T);
 Y1G=max(K-K*exp(r*T-T*mean(diag(V))/2+sqrt(T)*mean(X)),0)*exp(-r*T);
 Y2G=max(K-K*exp(r*T-T*mean(diag(V))/2-sqrt(T)*mean(X)),0)*exp(-r*T);
 Sxy=Sxy+sum((Y1+Y2)/2.*((Y1G+Y2G)/2-mPG));
 Sxx=Sxx+sum(((Y1G+Y2G)/2-mPG).^2);
 PA(k)=mean((Y1+Y2)/2);
 PG(k)=mean((Y1G+Y2G)/2);
end
beta=Sxy/Sxx;
UA=mean(PA);
UAC=mean(PA-beta*(PG-mPG));