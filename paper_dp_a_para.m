%variance lambda
% % implement non-parametric clustering method for change-point detection
% % when the number of clusters are fixed
% 
clear
clc
rng('shuffle')
tic
parpool(24)
% %***********************************************************
% %simulate data
u_v=[150 300];% two possible values of change-points
K_d=2;m=40; l1=0.25; l2=0.1; % true value of lam1, lam2
B=200;% generate B datasets 
%-------
%---------------
 cover_p=zeros(K_d+2,1);% whether the CI covers the true parameter, tau, l1,l2, l2
true_v=[u_v 1000*l1 1000*l2];
%-----------
estimate_all=zeros(K_d+4,B);%change-point of each cluster, lam1,lam2, # of clusters, % of correct assigned subjects,  
top3_all=zeros(2,2,B);
alpha0_m = zeros(1,B);
parfor s_b=1:B
 [tau_h,lam_1h,lam_2h,no_k_c,p_h,cover_p_f,alpha0_m(s_b)]=dp_MCMC_a_f(u_v,m,l1,l2,K_d,@latent_simu_f_lamj);   
 estimate_all(:,s_b)=[tau_h';1000*lam_1h;1000*lam_2h;no_k_c;p_h];
 cover_p=cover_p+cover_p_f;
end
% % %inference
 RMSE=zeros(K_d+2,1);Bias=zeros(K_d+2,1);
for i=1:K_d+2
RMSE(i)= sqrt((estimate_all(i,:)-true_v(i))*(estimate_all(i,:)-true_v(i))'/B);
Bias(i)=sum(estimate_all(i,:)-true_v(i))/(B*true_v(i))*100;
end
result=[true_v' mean(estimate_all(1:K_d+2,:),2) RMSE abs(Bias) 100*cover_p/B]
p_correct_no=100*sum(estimate_all(K_d+3,:)==K_d)/B;%percentage of correctly estimated number of clusters
p_correct_cluster=100*mean(estimate_all(K_d+4,:));%average percentage of correctly clustered subjects
p_correct = [p_correct_no p_correct_cluster]
toc
alpha0 = mean(alpha0_m);
t_time = 2.7
save dp_para_a.mat result p_correct alpha0 t_time;
%load dp_para_top3.mat;
% per_top = top3_all(2,top3_all(1,:,:)==K_d)%Percentages of the correct cluster no
% hist(100*per_top)
% xlabel('Percentages of the correct cluster no')
% ylabel('Frequency')