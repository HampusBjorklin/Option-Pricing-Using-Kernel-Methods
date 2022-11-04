function [VegaAC,VegaA]=BSeuPut5DbasketVega_MCA(param)
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
A=chol(V)';
tic
X=A*randn(n,N);
for j=1:5
 Y1(:,j)=-w(j)*(X(j,:)/s(j)*sqrt(T)-s(j)*T).*exp(r*T-T*V(j,j)/2+sqrt(T)*X(j,:)).*(K-S0*exp(r*T-T*repmat(diag(V),1,N)/2+sqrt(T)*X)>0)*exp(-r*T);
 Y2(:,j)=-w(j)*(-X(j,:)/s(j)*sqrt(T)-s(j)*T).*exp(r*T-T*V(j,j)/2-sqrt(T)*X(j,:)).*(K-S0*exp(r*T-T*repmat(diag(V),1,N)/2-sqrt(T)*X)>0)*exp(-r*T);
 Y1G(:,j)=-w(j)*(X(j,:)/s(j)*sqrt(T)-s(j)*T).*exp(r*T-T*w'*diag(V)/2+sqrt(T)*w'*X).*(K>K*exp(r*T-T*w'*diag(V)/2+sqrt(T)*w'*X))*exp(-r*T);
 Y2G(:,j)=-w(j)*(-X(j,:)/s(j)*sqrt(T)-s(j)*T).*exp(r*T-T*w'*diag(V)/2-sqrt(T)*w'*X).*(K>K*exp(r*T-T*w'*diag(V)/2-sqrt(T)*w'*X))*exp(-r*T);
end 
mVG=BSVegaPut(K,K,T,r,r-w'*diag(V)/2,sqrt(w'*V*w),s,w,rho);
Sxy=(Y1+Y2)'/2*((Y1G+Y2G)/2-repmat(mVG,N,1));
Sxx=((Y1G+Y2G)/2-repmat(mVG,N,1))'*((Y1G+Y2G)/2-repmat(mVG,N,1));
beta=Sxy\Sxx;
e0=max(std((Y1+Y2)/2-((Y1G+Y2G)/2-repmat(mVG,N,1))*beta))/sqrt(N);
dt=toc;
nrep=ceil((2*e0/param)^2);
if nrep>1
 disp(['Estimated running time: '  num2str(dt*(nrep-1)) ' s.']);   
end    
PA(1,:)=mean((Y1+Y2)/2);
PG(1,:)=mean((Y1G+Y2G)/2);
for k=2:nrep
 X=A*randn(n,N);
 for j=1:5
  Y1(:,j)=-w(j)*(X(j,:)/s(j)*sqrt(T)-s(j)*T).*exp(r*T-T*V(j,j)/2+sqrt(T)*X(j,:)).*(K-S0*exp(r*T-T*repmat(diag(V),1,N)/2+sqrt(T)*X)>0)*exp(-r*T);
  Y2(:,j)=-w(j)*(-X(j,:)/s(j)*sqrt(T)-s(j)*T).*exp(r*T-T*V(j,j)/2-sqrt(T)*X(j,:)).*(K-S0*exp(r*T-T*repmat(diag(V),1,N)/2-sqrt(T)*X)>0)*exp(-r*T);
  Y1G(:,j)=-w(j)*(X(j,:)/s(j)*sqrt(T)-s(j)*T).*exp(r*T-T*w'*diag(V)/2+sqrt(T)*w'*X).*(K>K*exp(r*T-T*w'*diag(V)/2+sqrt(T)*w'*X))*exp(-r*T);
  Y2G(:,j)=-w(j)*(-X(j,:)/s(j)*sqrt(T)-s(j)*T).*exp(r*T-T*w'*diag(V)/2-sqrt(T)*w'*X).*(K>K*exp(r*T-T*w'*diag(V)/2-sqrt(T)*w'*X))*exp(-r*T);
 end
 Sxy=Sxy+(Y1+Y2)'/2*((Y1G+Y2G)/2-repmat(mVG,N,1));
 Sxx=Sxx+((Y1G+Y2G)/2-repmat(mVG,N,1))'*((Y1G+Y2G)/2-repmat(mVG,N,1));
 PA(k,:)=mean((Y1+Y2)/2);
 PG(k,:)=mean((Y1G+Y2G)/2);
end
beta=Sxy/Sxx;
VegaA=mean(PA,1);
VegaAC=mean(PA,1)-(mean(PG,1)-mVG)*beta;

function [Vega]=BSVegaPut(S0,K,T,r,mu,sigma,s,w,rho)
d1=1./(sigma.*sqrt(T)).*(log(S0./K)+mu.*T)+sigma.*sqrt(T);
d2=d1-sigma.*sqrt(T);
Vega=exp(T.*(sigma.^2/2+mu-r)).*S0.*(normcdf(d1)-1)*T*(s'*diag(w)*rho*diag(w)-s'*diag(w))...
+exp(T.*(sigma.^2/2+mu-r)).*S0.*normpdf(d1)*...
(s'*diag(w)*rho*diag(w)/sigma*(-1./(sigma^2.*sqrt(T)).*(log(S0./K)+mu.*T)+sqrt(T))-s'*diag(w)*sqrt(T)/sigma)...
-exp(-r.*T).*K.*normpdf(d2)*(s'*diag(w)*rho*diag(w)/sigma*(-1./(sigma^2.*sqrt(T)).*(log(S0./K)+mu.*T))-s'*diag(w)*sqrt(T)/sigma);



