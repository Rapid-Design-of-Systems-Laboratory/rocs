K                 = 200;  
option.Mt         = cat(3 , [0.6 ; 6] , [1 ; -10] ,[10 ; -1] , [ 0 ; 10] , [1 ; -3] , [5 ; 5]);                                       %(d x 1 x L)
option.St         = corr2cov(cat(3 , [1 0.9 ; 0.9 1] , [1 -0.9 ; -0.9 1] , [2 0 ; 0 2] , [2 0 ; 0 2] , [2 0 ; 0 2] , [2 0 ; 0 2]  )); %(d x d x L)
option.pt         = cat(3 , [0.1] , [0.1]  , [0.2] , [0.2] , [0.2] , [0.2] );                                                         %(1 x 1 x L)

Z                 = mvgmmrnd(option.Mt , option.St , option.pt , K);
option.L          = size(option.Mt , 3);
option.Mmin       = [];
option.Mmax       = [];
option.Rvarmin    = 0.8;
option.Rvarmax    = [];
option.Rcorrmin   = -0.95;
option.Rcorrmax   = 0.95;
option.pmin       = 0;
option.pmax       = 1;
option.N          = 100;
option.rho        = 20*10e-3;
option.alpha      = 0.7;
option.beta       = 0.95;
option.q          = 9;
option.epsi       = 10e-2;
option.h          = 10;
option.cons       = 5;
option.T_max      = 2500;
option.verbose    = 1;

[M_opt , S_opt , p_opt , loglt , T] = ce_mvgm(Z , option);
logl_true                           = -loglikelihood(Z , option.Mt , option.St , option.pt);        
[xt , yt]                           = ndellipse(option.Mt , option.St);  
[x , y]                             = ndellipse(M_opt , S_opt); 
figure(1)
plot(Z(1 , :) , Z(2 , :) , '+' ,  xt , yt , 'r' , x , y , 'g' , 'linewidth' , 2 , 'markersize', 2);
grid on
figure(2)
plot((1:T) , loglt , (1:T) , logl_true(: , ones(1 , T)) , 'g' , 'linewidth' , 2)
grid on 
axis([1 , T , 0.8*logl_true , 1.2*max(loglt) ])
legend(['Logl(t)'] , ['Logl_{true}'],1)