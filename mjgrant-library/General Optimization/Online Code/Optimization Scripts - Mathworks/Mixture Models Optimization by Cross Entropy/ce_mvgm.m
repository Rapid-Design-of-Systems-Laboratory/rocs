function [M_opt , S_opt , p_opt , loglt , T] = ce_mvgm(Z , option);

% Stochastic Maximum Likelihood optimization with a Cross Entropy (CE) algorithm for
% gaussian mixture estimation application
%
% Usage
% ------
%
%  [M_opt , S_opt , p_opt , loglt , T] = gauss_ce(Z , option);
%
% Inputs
% -------
%
% Z                       Data (d x K)
% option                  (optional) structure for the option of the CE algorithm
%
%   option.L              Number of mixture (default = 1)
%   option.Mmin           Minimum' means value (default = [] = no minimum)
%   option.Mmax           Maximum' means value (default = [] = no maximum)
%   option.Rvarmin        Minimum' variance value (default = 0.8 ,  [] = no minimum)
%   option.Rvarmax        Maximum' variance value (default = [] = no maximum)
%   option.Rcorrmin       Minimum' correlation value (default = -1)
%   option.Rcorrmax       Maximum' correlation value (default = 1)
%   option.pmin           Minimum' weight value (default = 0)
%   option.pmax           Maximum' weight value (default = 1)
%   option.N              Number of gaussian mixture to generate for each iteration (default 150)
%   option.rho            Threshold for selectionning the best rho*N mixture (default 15*10e-3)
%   option.alpha          Weight for updating Gaussian's mean of parameters (default 0.8)
%   option.beta           Starting weight value for updating variance (default 0.95)
%   option.q              Cooling temperature parameter for updating beta;
%   option.epsi           Minimum loglikelihood difference for extra variance injection (default 10e-2)
%   option.h              Extra variance injection parameter (default 5)
%   option.cons           Indicate how many CE iterations to wait until declaring a solution is found (default 5)
%   option.T_max          Maximum number of iterations (default 2500)
%   option.verbose        1 for real time plot if d = 2 (default 0);
%   option.Mt             True Mean of the data (d x 1 x M). option.verbose = 1 is requiered;  
%   option.St             True variance (d x d x M). option.verbose = 1 is requiered;  
%   option.pt             True weight (1 x 1 x M). option.verbose = 1 is requiered;  
%
%
% Outputs
% --------
% 
% M_opt                   Optimal mean. a (d x 1 x L) matrix.
% S_opt                   Optimal variance. a (d x d x L) matrix.
% p_opt                   Optimal path found. a (1 x 1 x L) ) matrix.
% loglt                   LogLikelihood of the global solution versus iteration
% T                       Number of CE iteration ( < option.T_max))
%
%
% Example
%---------
%
%  K                 = 200;  
%  option.Mt         = cat(3 , [0.6 ; 6] , [1 ; -10] ,[10 ; -1] , [ 0 ; 10] , [1 ; -3] , [5 ; 5]);                                       %(d x 1 x L)
%  option.St         = corr2cov(cat(3 , [1 0.9 ; 0.9 1] , [1 -0.9 ; -0.9 1] , [2 0 ; 0 2] , [2 0 ; 0 2] , [2 0 ; 0 2] , [2 0 ; 0 2]  )); %(d x d x L)
%  option.pt         = cat(3 , [0.1] , [0.1]  , [0.2] , [0.2] , [0.2] , [0.2] );                                                         %(1 x 1 x L)
%  Z                 = mvgmmrnd(option.Mt , option.St , option.pt , K);
%  option.L          = size(option.Mt , 3);
%  option.Mmin       = [];
%  option.Mmax       = [];
%  option.Rvarmin    = 0.8;
%  option.Rvarmax    = [];
%  option.Rcorrmin   = -0.95;
%  option.Rcorrmax   = 0.95;
%  option.pmin       = 0;
%  option.pmax       = 1;
%  option.N          = 100;
%  option.rho        = 20*10e-3;
%  option.alpha      = 0.7;
%  option.beta       = 0.95;
%  option.q          = 9;
%  option.epsi       = 10e-2;
%  option.h          = 10;
%  option.cons       = 5;
%  option.T_max      = 2500;
%  option.verbose    = 1;
%  [M_opt , S_opt , p_opt , loglt , T] = ce_mvgm(Z , option);
%  logl_true                           = -loglikelihood(Z , option.Mt , option.St , option.pt);        
%  [xt , yt]                           = ndellipse(option.Mt , option.St);  
%  [x , y]                             = ndellipse(M_opt , S_opt); 
%  figure(1)
%  plot(Z(1 , :) , Z(2 , :) , '+' ,  xt , yt , 'r' , x , y , 'g' , 'linewidth' , 2 , 'markersize', 2);
%  grid on
%  figure(2)
%  plot((1:T) , loglt , (1:T) , logl_true(: , ones(1 , T)) , 'g' , 'linewidth' , 2)
%  grid on 
%  axis([1 , T , 0.8*logl_true , 1.2*max(loglt) ])
%  legend(['Logl(t)'] , ['Logl_{true}'],1)
%

%
% Author : Sébastien PARIS (sebastien.paris@lsis.org)
%
% Ver 1.0 (03/04/2005)
%
% Ref    : "Global Likelihood optimization via the cross-entropy method with an application to mixture models",
% Zdravko Boetev, Dirk. P. Kroese, Proceeding of the 2004 Winter simulation conference.
%
% see    http://iew3.technion.ac.il/CE/about.php




[d , K] = size(Z);

if ( (nargin < 2) || isempty(option))
    
    option.L          = 1;
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
    option.h          = 3;
    option.cons       = 5;
    option.T_max      = 2500;
    option.verbose    = 0;
    
end

L       = option.L;

alpha   = option.alpha;


if (option.verbose == 1)
    
    if(d ~= 2)
        
        disp('Verbose only for d = 2');
        
        option.verbose == 0;
        
    end
    
    if (~any(strcmp(fieldnames(option) , 'Mt')))
        
        disp('option.Mt must be provided tp verbose, switching verbose = 0');
        
        option.verbose == 0;
        
    else
        
        fig1        = figure(1);
        
        set(fig1 , 'renderer' , 'zbuffer');
        
        set(fig1 , 'doublebuffer' , 'on');
        
        set(gca , 'nextplot' , 'replace');
        
        H                = [1 0 ; 0 1];
        
        [xt , yt]        = ndellipse(option.Mt , option.St , H);  
        
        maxi             = 1.6*(max(Z , [] , 2));
        
        mini             = 1.6*(min(Z , [] , 2));
        
    end
    
end


Nrho                      = round(option.rho*option.N);

cteNrho                   = 1/Nrho;

cte1Nrho                  = 1/(Nrho - 1);

Od                        = ones(d , 1);

OL                        = ones(1 , L);

ONrho                     = ones(1 , Nrho);

loglt                     = zeros(1 , option.T_max);


ind_Nrho                  = (1 : Nrho);


ind_y                     = (1 : d)';

i1                        = (1 : d+1 :d*d)';

i2                        = (0:L - 1)*d*d;

indice1                    = reshape((i1(: , OL , :) + i2(Od , :)) , d*L , 1);

indice2                   = (1:d*d*L);

indice2(indice1)          = [];


%%==== Initial draw  =====%%

%a) Initial Mean %  M ~ N(Ma , Mb , Mmin , Mmax) with Ma = E[Z], Mb =%(4*4)*Cov[Z]

tpM                     = sum(Z , 2)/K;

Ma                      = tpM(: , : , OL);


res                     = (Z - tpM(: , ones(1 , K)));

tpR                     = cov2corr(res*res'/(K - 1));


temp                    = (4*4)*diag(tpR);

Mb                      = temp(: , : , OL);

if (isempty(option.Mmin))
    
    Mmin     = reshape(option.Mmin , 0 , 0 , L);
    
else
    
    Mmin     = option.Mmin(Od , : , OL);
    
end

if (isempty(option.Mmax))
    
    Mmax     = reshape(option.Mmax , 0 , 0 , L);
    
else
    
    Mmax     = option.Mmax(Od , : , OL);
    
end


%b) Initial Covariance %   R ~ N(Ra , Rb , Rmin , Rmax)


Ra            = tpR(: , : , OL)*0.5;

Rb            = 3*abs(Ra)*0.5;


Rmin          = zeros(d , d , L);

Rmax          = zeros(d , d , L);


if(isempty(option.Rvarmin))
    
    Rmin(indice1) = Ra(indice1)/9 - 2*Ra(indice1)/9;
    
else
    
    Rmin(indice1) = option.Rvarmin;
    
end

if (isempty(option.Rcorrmin))
    
    Rmin(indice2) = -1;
    
else
    
    
    Rmin(indice2) = option.Rcorrmin;
    
end


if(isempty(option.Rvarmax))
    
    Rmax(indice1) = Ra(indice1)/9 + 2*Ra(indice1)/9;
    
else
    
    Rmax(indice1) = option.Rvarmax;
    
end



if (isempty(option.Rcorrmax))
    
    Rmax(indice2) = 1;
    
else
    
    Rmax(indice2) = option.Rcorrmax;
    
end



pa         = ones(1 , 1 , L)/L;    % p ~ N(pa , pb² , pmin , pmax) with pa centered on 1/L and pb variance of a uniform quantization of [0 , 1] by L

pb         = (1/(12*L))*ones(1 , 1 , L);


if (isempty(option.pmax))
    
    pmax     = reshape(option.pmax , 0 , 0 , L);
    
else
    
    
    pmax     = option.pmax(: , : , OL);
    
end


if (isempty(option.pmin))
    
    pmin     = reshape(option.pmin , 0 , 0 , L);
    
else
    
    
    pmin     = option.pmin(: , : , OL);
    
end



t        = 1;

logl_opt = +inf;

cons     = 0;

while( (t < option.T_max) & (cons < option.cons) )
    
    [M , R , p]               = sample_gaussian_mixture(Ma , Mb , Mmin , Mmax , Ra , Rb , Rmin , Rmax , pa , pb , pmin , pmax ,  option.N);
    
    S                         = corr2cov(R);                                       
    
    [logl , indice]           = sort(-loglikelihood(Z , M , S , p));
    
    logl_new                  = logl(1);
    
    
    if (max(Rb(indice1)) < option.epsi);
        
        Rb(indice1) = Rb(indice1) + abs(logl_new - logl_old)*option.h;
        
        cons        = cons + 1;
        
    end
    
    I                         = indice(ind_Nrho);
    
    tp_M                      = M(: , : , : , I);
    
    tp_R                      = R(: , : , : , I);
    
    tp_p                      = p(: , : , : , I);
    
    
    
    beta                      = option.beta - option.beta*(1 - 1/t)^(option.q);
    
    
    mean_M                    = sum(tp_M , 4)*cteNrho;
    
    res_M                     = (tp_M - mean_M(: , : , : , ONrho));
    
    Ma                        = alpha*mean_M + (1 - alpha)*Ma;
    
    Mb                        = beta*(sum(res_M.*res_M , 4)*cte1Nrho) + (1 - beta)*Mb;
    
    
    
    mean_R                    = sum(tp_R , 4)*cteNrho;
    
    res_R                     = (tp_R - mean_R(: , : , : , ONrho));
    
    Ra                        = alpha*mean_R + (1 - alpha)*Ra;
    
    Rb                        = beta*(sum(res_R.*res_R , 4)*cte1Nrho) + (1 - beta)*Rb;
    
    
    mean_p                    = sum(tp_p , 4)*cteNrho;
    
    res_p                     = (tp_p - mean_p(: , : , : , ONrho));
    
    pa                        = alpha*mean_p + (1 - alpha)*pa;
    
    pb                        = beta*(sum(res_p.*res_p , 4)*cte1Nrho) + (1 - beta)*pb;
    
    if (logl_new < logl_opt)
        
        M_opt                     = M(: , : , : , I(1));
        
        R_opt                     = R(: , : , : , I(1));
        
        S_opt                     = S(: , : , : , I(1));
        
        p_opt                     = p(: , : , : , I(1));
        
        logl_opt                  = logl_new;
        
    end
    
    logl_old                  = logl_new;
    
    loglt(t)                  = logl_opt;
        
    
    t                         = t + 1;
    
    
    %%===========  Display ===========%%
    
    if (option.verbose)
        
        fig1             = figure(1);
        
        [x , y]          = ndellipse(M_opt , S_opt , H );   
        
        plot(Z(1 , :) , Z(2 , :) , '+' ,  xt , yt , 'r' , x , y , 'g' , 'linewidth' , 2 , 'markersize', 2);
        
        title(sprintf('t = %d, current sol = %4.4f, global sol = %4.4f, d = %d' , t , logl(1) , logl_opt , cons))
        
        axis([mini(1) , maxi(1) , mini(2) , maxi(2)])
        
        drawnow
        
    else
        
        disp(sprintf('t = %d, current sol = %4.4f, global sol = %4.4f, d = %d' , t , logl(1) , logl_opt , cons));
        
    end  
    
end

loglt(t:option.T_max)  = [];

T                      = t - 1;
