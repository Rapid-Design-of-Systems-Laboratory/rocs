function [mu , sigma , p] = sample_gaussian_mixture(mua , mub , mumin , mumax , ...
                                           sigmaa , sigmab , sigmamin , sigmamax , ...
                                           pa , pb , pmin , pmax ,  N);
 
% [mu , sigma , p] = sample_gaussian_mixture(mua , mub , mumin , mumax , ...
%                                            sigmaa , sigmab , sigmamin , sigmamax , ...
%                                            pa , pb , pmin , pmax ,  N);
%                                       
%
% d         = 3;
% M         = 4;
% N         = 1000;
% mua       = randn(d , 1 , M);
% mub       = rand(d , 1 , M);
% mumin     = [];
% mumax     = [];
%
% sigmaa    = randn(d , d , M);
% sigmab    = rand(d , d , M);
% sigmamin  = zeros(d , d , M);
% sigmamax  = zeros(d , d , M);
% i1        = (1:d+1:d*d)';
% i2        = (0:M-1)*d*d;
% index     = reshape((i1(: , ones(1 , M) , :) + i2(ones(d , 1) , :)) , d*M , 1);
% indice    = (1:d*d*M);
% indice(index) = [];
%
% sigmamin(indice) = -0.95;
% sigmamin(index)  = 0.7;
% sigmamax(indice) = +0.95;
% sigmamax(index)  = 10;
%
%
% pa       = rand(1 , 1 , M);
% pb       = rand(1 , 1 , M);
% pmin     = zeros(1 , 1 , M);
% pmax     = ones(1 , 1 , M);
% 
%


[d , v , M]  = size(mua);

ON           = ones(1 , N);

%%% Sample mu %%% N(mua , mub , mumin , mumax)

mu          = randnt(mua(: , : , : , ON) , mub(: , : , : , ON) , mumin(: , : , : , ON) , mumax(: , : , : , ON));

%%% Sample sigma %%% N(mua , mub , mumin , mumax)

sigma       = randnt(sigmaa(: , : , : , ON) , sigmab(: , : , : , ON) , sigmamin(: , : , : , ON) , sigmamax(: , : , : , ON));

%%% Sample p %%% N(pa , pb , pmin , pmax)

p          = randnt(pa(: , : , : , ON) , pb(: , : , : , ON) , pmin(: , : , : , ON) , pmax(: , : , : , ON));

%%%% Normalization %%%%

sump       = sum(p , 3);

p          = p./sump(: , : , ones(1 , M) , :);  %(1 x 1 x M x N);