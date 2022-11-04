function [DeltaAC,DeltaA]=BSeuPut5DbasketDelta_MCA(param)
K=1;   
r=0.05;
T=1;
N=1e7;
n=5;
s=[0.518; 0.648; 0.623; 0.570; 0.530];
w= [0.381; 0.065; 0.057; 0.270; 0.227];
rho=[1.00 0.79 0.82 0.91 0.84 
      0.79 1.00 0.73 0.80 0.76
      0.82 0.73 1.00 0.77 0.72
      0.91 0.80 0.77 1.00 0.90
      0.84 0.76 0.72 0.90 1.00];
V=diag(s)*rho*diag(s);
S0=K*w'; 
A=sqrtm(V);
tic
X=A*randn(n,N);
for j=1:5
 Y1(:,j)=-w(j)*exp(r*T-T*V(j,j)/2+sqrt(T)*X(j,:)).*(K-S0*exp(r*T-T*repmat(diag(V),1,N)/2+sqrt(T)*X)>0)*exp(-r*T);
 Y2(:,j)=-w(j)*exp(r*T-T*V(j,j)/2-sqrt(T)*X(j,:)).*(K-S0*exp(r*T-T*repmat(diag(V),1,N)/2-sqrt(T)*X)>0)*exp(-r*T);
end 
Y1G=-exp(r*T-T*w'*diag(V)/2+sqrt(T)*w'*X).*(K>K*exp(r*T-T*w'*diag(V)/2+sqrt(T)*w'*X))*exp(-r*T);
Y2G=-exp(r*T-T*w'*diag(V)/2-sqrt(T)*w'*X).*(K>K*exp(r*T-T*w'*diag(V)/2-sqrt(T)*w'*X))*exp(-r*T);
[mDG]=DeltaBSPut(K,K,T,r,r-w'*diag(V)/2,sqrt(w'*V*w));
Sxy=sum((Y1+Y2)/2.*repmat(((Y1G+Y2G)/2-mDG)',1,5));
Sxx=sum(((Y1G+Y2G)/2-mDG).^2);
beta=Sxy/Sxx;
e0=max(std((Y1+Y2)/2-repmat(beta,N,1).*repmat(((Y1G+Y2G)/2-mDG)',1,5)))/sqrt(N);
dt=toc;
nrep=ceil((2*e0/param)^2);
if nrep>1
 disp(['Estimated running time: '  num2str(dt*(nrep-1)) ' s.']);   
end    
PA(1,:)=mean((Y1+Y2)/2);
PG(1)=mean((Y1G+Y2G)/2);
for k=2:nrep
 X=A*randn(n,N);
 for j=1:5
  Y1(:,j)=-w(j)*exp(r*T-T*V(j,j)/2+sqrt(T)*X(j,:)).*(K-S0*exp(r*T-T*repmat(diag(V),1,N)/2+sqrt(T)*X)>0)*exp(-r*T);
  Y2(:,j)=-w(j)*exp(r*T-T*V(j,j)/2-sqrt(T)*X(j,:)).*(K-S0*exp(r*T-T*repmat(diag(V),1,N)/2-sqrt(T)*X)>0)*exp(-r*T);
 end
 Y1G=-exp(r*T-T*w'*diag(V)/2+sqrt(T)*w'*X).*(K>K*exp(r*T-T*w'*diag(V)/2+sqrt(T)*w'*X))*exp(-r*T);
 Y2G=-exp(r*T-T*w'*diag(V)/2-sqrt(T)*w'*X).*(K>K*exp(r*T-T*w'*diag(V)/2-sqrt(T)*w'*X))*exp(-r*T);
 Sxy=Sxy+sum((Y1+Y2)/2.*repmat(((Y1G+Y2G)/2-mDG)',1,5));
 Sxx=Sxx+sum(((Y1G+Y2G)/2-mDG).^2);
 PA(k,:)=mean((Y1+Y2)/2);
 PG(k)=mean((Y1G+Y2G)/2);
end
beta=Sxy/Sxx;
DeltaA=mean(PA,1);
DeltaAC=mean(PA,1)-beta*mean(PG-mDG);

function [Delta]=DeltaBSPut(S0,K,T,r,mu,sigma);
d1=1./(sigma.*sqrt(T)).*(log(S0./K)+mu.*T)+sigma.*sqrt(T);
Delta=exp(T.*(sigma.^2/2+mu-r)).*S0.*(normcdf(d1)-1);

