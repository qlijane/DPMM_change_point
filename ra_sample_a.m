function alpha0  = ra_sample_a(a0,b0,tau,m,b_alpha, a_alpha )
%rejection sampleing of alpha0
%   Detailed explanation goes here
ginv_alpha=b_alpha-a_alpha;
x_alpha=a_alpha:0.1:b_alpha;%alpha_0 is between a_alpha and b_alpha
   i=1;
                logf=zeros(1,length(x_alpha));
            for u=x_alpha
                logf(i)=logfa(a0,b0,u,tau,m);
                i=i+1;
            end
            adju_c= log(1/ginv_alpha)-(max(logf));
          %  logf_a=logf+adju_c;% adjust by constant for computational purpose
           % c_u=exp(max(logf_a))*ginv_alpha;% how to choose c? 
           c_u=1;
                  accept = false;
                while accept == false
                  u = rand();
                  v = a_alpha + (b_alpha -a_alpha )*rand();              
                    f1=logfa(a0,b0,v,tau,m)+adju_c;
                    if c_u*u <= exp(f1)*ginv_alpha
                        alpha0 = v;
                        accept = true;
                    end
                end