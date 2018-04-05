function [tau_h,lam_1h,lam_2h,no_k_c,p_h,cover_p,alpha0_m]=dp_MCMC_a_f(u_v,m,l1,l2,K_d,latent_simu_f_lamj)
 Bt=20000;%# of samples for Gibbs
Burn_in=10000;% burn in time
lag=5;%every 5th result is used
discard_CI=100;%only keep the values of +- 50 within the tau_h for each cluster to calculate CI; 
 cover_p=zeros(K_d+2,1);% whether the CI covers the true parameter, tau, l1,l2, l2
true_v=[u_v 1000*l1 1000*l2];
a1=0.5;b1=1;a2=0.1;b2=1;%prior of lam1,lam2,mean is a1/b1
l_tau=70;u_tau=350;ginv=u_tau-l_tau;% for rejection sampling of tau, the range of tau
x_tau=l_tau:2:u_tau;%tau is between l_tau and u_tau
alpha0=zeros(1,Bt);
 a0=0.0004;b0=2000;% 0.00045, P1 = 91.5%,  % prior of alpha_0, make sure alpha_0 is small
 a_alpha=0.0001;b_alpha=0.001;ginv_alpha=b_alpha-a_alpha;% for rejection sampling of alpha_0
 x_alpha=a_alpha:0.00005:b_alpha;%alpha_0 is between a_alpha and b_alpha
lam_1_gib=zeros(m,Bt);lam_2_gib=zeros(m,Bt);%log_k_gib=zeros(m_n,Bt);
lam_1_gib(:,1)=0.05; lam_2_gib(:,1)=0.15; 
    no_k=zeros(1,Bt);
    alpha0(1)=0.0004;
[z,Nj,C,tau_true_index]=latent_simu_f_lamj(u_v,m,l1,l2,K_d);%___________data simulation end_______________
 %________________________Gibbs start__________________________
%---------------------------
% initialize randomly
zj=unidrnd(K_d,m,1);tau_gib=u_v(1)*ones(m,Bt);
% zj is the index of tau_j in uk;
for k_t=2:K_d
tau_gib(zj==k_t,1)=u_v(k_t);% some drivers change at u1, others change at u2;
end
%_________initialize end__________________
    for t_g=2:Bt
        t_g   
% update tau
tau=tau_gib(:,t_g-1);%current tau
         for j=1:m   
             C_j=C(j);Nj_j=Nj(j);z_j=z(j,:);
              lam_1_gib_t = lam_1_gib(j,t_g-1);
              lam_2_gib_t = lam_2_gib(j,t_g-1);
          logf=logLj(lam_1_gib_t,lam_2_gib_t,x_tau,C_j,Nj_j,z_j,0);
           adj_c= -(max(logf));  
            % adjust by constant for computational purpose
           q0_adj=integral(@(x)Lj(lam_1_gib_t,lam_2_gib_t,x,C_j,Nj_j,z_j,adj_c),l_tau,u_tau,'AbsTol',1e-6,'RelTol',1e-3);
          tau_new=tau;tau_new(j)=[];% vector of tau except tau_j
          L_adj=Lj(lam_1_gib_t,lam_2_gib_t,tau_new,C_j,Nj_j,z_j,adj_c);sum_Ladj=sum(L_adj);
             rand_p=(alpha0(t_g-1)*q0_adj+sum_Ladj)*rand;% decide how to get tau_j              
          if rand_p<=sum_Ladj
                idx = sum(rand_p> cumsum(L_adj)) + 1;
                tau(j)=tau_new(idx);% sample from tau_(-j)
          % sample from H
          else
             tau(j)= ra_sampl_tau(l_tau,u_tau,ginv,adj_c,z_j,C_j,Nj_j,lam_1_gib_t,lam_2_gib_t);
          end
          %-------------------        
         end
         tau_gib(:,t_g)=tau;% current tau vector
         no_k(t_g)=length(unique(tau));
% % update lam_1
            for j=1:m
                Nj_1=sum(z(j,1:Nj(j))<=tau(j));            
            lam_1_gib(j,t_g)=gamrnd(a1+Nj_1,1/(b1+tau(j)),1);
% update lam_2
            lam_2_gib(j,t_g)=gamrnd(a2+Nj(j)-Nj_1,1/(b2+C(j)-tau(j)),1);
            end
            % update alpha_0
  k_a=length(unique(tau))
 alpha0(t_g)=a_west1995(alpha0(t_g-1),a0,b0,m,k_a);
    end 
% %________________________Gibbs end__________________________
% Decide the number of clusters 
no_k_c=floor(mean(no_k(Burn_in:lag:end)));% # of clusters
tau_burn=tau_gib(:,Burn_in:lag:end);
alpha0_m = mean(alpha0(Burn_in:lag:end));
%cluster assignment by confusion matrix
Confusion_v=zeros(1,m*(m-1)/2);i_0=1;%count the number of times each datapoint was assigned to the same compoenent with another datapoint
for i=1:m-1
    for j=i+1:m
       Confusion_v(i_0)=sum(tau_burn(i,:)==tau_burn(j,:));
      % Confusion_m(i,j)=Confusion_v(i_0);
    i_0=i_0+1;
    end
end
membership_h=zeros(m,1);
pppm=1-Confusion_v/(Bt-Burn_in);
Z_p = linkage(pppm,'complete');
membership = cluster(Z_p,'maxclust',K_d); 
   C_tau=zeros(1,K_d);
    for i=1:K_d
    C_tau(i)=mean2(tau_burn(membership==i,:));
    end
    [tau_h,I_d]=sort(C_tau);%cluster parameters ascending
    for i=1:K_d
     membership_h(membership==i)=I_d(i);% membership in order
    end
p_h=sum(membership_h==tau_true_index)/m;%percengate of correctly clustered subjects
 lam_1_use=lam_1_gib(:,Burn_in:lag:end);
  lam_2_use=lam_2_gib(:,Burn_in:lag:end);
  lam_1_use_v=lam_1_use(:);
  lam_2_use_v=lam_2_use(:);
 lam_1h=mean(lam_1_use_v);
lam_2h=mean(lam_2_use_v);
 tau_CI1=tau_burn(membership_h==1,:);
tau_CI2=tau_burn(membership_h==2,:);
tau_CI1_v=tau_CI1(tau_CI1<tau_h(1)+discard_CI );%make it a vector, discard some values for CI
tau_CI2_v=tau_CI2(tau_CI2>tau_h(2)-discard_CI);

CI=[ hpdr(tau_CI1_v)  hpdr(tau_CI2_v) 1000*quantile(lam_1_use_v,[0.025,0.975]) 1000*quantile(lam_2_use_v,[0.025,0.975])];
    for i_ci=1:K_d+2
        if true_v(i_ci)<=CI(2*i_ci) && true_v(i_ci)>=CI(2*i_ci-1)
      cover_p(i_ci)=1;    
        end
    end
%For each data set, list the top 3 number of clusters and the percentage;
no_k_h=no_k(Burn_in:lag:end);
unique_no=unique(no_k_h);
len_unique=length(unique_no);
per_unique=zeros(0,len_unique);
  for i=1:len_unique
    per_unique(i)=sum(no_k_h==unique_no(i))/length(no_k_h);  
  end
  [s_per,s_index] = sort(per_unique,'descend');
  top3=[unique_no(s_index(1:2));s_per(1:2)];
  %top3 end
end
