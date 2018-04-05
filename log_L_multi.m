function logr= log_L_multi(lam1,lam2,tau,cj,Nj_r,t,adj)% Compute the log likelihood of multiple individuals
%log L of multile drivers
% tau is a vector of all the change-points, cj is vector, t is a matrix with each row for a individual with the same change-point 
% Nj is a vector recording the # of events for each driver
% K is the # of drivers
N=sum(Nj_r);K=length(cj);%# of drivers
Nj_1=0;
for j=1:K
Nj_1=Nj_1+sum(t(j,1:Nj_r(j))<=tau(j));
end
logr=(lam2-lam1)*sum(tau)-lam2*sum(cj)+Nj_1*log(lam1)+(N-Nj_1)*log(lam2)+adj;
%r=exp(logr);
