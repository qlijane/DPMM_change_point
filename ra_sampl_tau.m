function u_s=ra_sampl_tau(a,b,ginv,adj_c,t,Cj,Nj,lam1,lam2)

%t is the event time for the jth driver;Nj is the # of events for the jth
%driver
% rejection and acceptance sampling from H

 adju_ch= log(1/ginv)+adj_c;
            c_u=1;% how to choose c? 
                  accept = false;
                while accept == false
                  u = rand();
                  v = a + (b-a)*rand();
                 if v>=a&&v<=b%V need to be between a and b
                     f1=Lj(lam1,lam2,v,Cj,Nj,t,adju_ch);
                     if c_u*u <= f1*ginv
                        u_s = v;
                        accept = true;
                    end
                 end 
                end