function [UAC,UA]=BSeuPut4DbasketU_MCA(param)
N=1e7; 
sigma=0.3; 
n=4;
V=sigma^2*[1 0.3 0.4 0.5
     0.3 1 0.2 0.25
     0.4 0.2 1 0.3
     0.5 0.25 0.3 1];
A=sqrtm(V);
r=0.06;
T=1;
K=40;
S0=K*[ones(1,n)]/n;
tic
X=A*randn(n,N);
Y1=max(K-S0*exp(r*T-T*repmat(diag(V),1,N)/2+sqrt(T)*X),0)*exp(-r*T);
Y2=max(K-S0*exp(r*T-T*repmat(diag(V),1,N)/2-sqrt(T)*X),0)*exp(-r*T);
Y1G=max(K-K*exp(r*T-T*mean(diag(V))/2+sqrt(T)*mean(X)),0)*exp(-r*T);
Y2G=max(K-K*exp(r*T-T*mean(diag(V))/2-sqrt(T)*mean(X)),0)*exp(-r*T);
[~,mPG]=BS(K,K,T,r,r-mean(diag(V))/2,sqrt(sum(V(:))/n^2));
Sxy=sum((Y1+Y2)/2.*((Y1G+Y2G)/2-mPG));
Sxx=sum(((Y1G+Y2G)/2-mPG).^2);
beta=Sxy/Sxx;
e0=std((Y1+Y2)/2-beta*((Y1G+Y2G)/2-mPG))/sqrt(N);
dt=toc;
nrep=ceil((2*e0/param)^2);
if nrep>1
 disp(['Estimated running time: '  num2str(dt*(nrep-1)) ' s.']);   
end
PA(1)=mean((Y1+Y2)/2)
PG(1)=mean((Y1G+Y2G)/2)
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