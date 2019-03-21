% demopsobehavior.m
% demo of the pso.m function
% the pso tries to find the minimum of the f6 function, a standard
% benchmark
%
% on the plots, blue is current position, green is Pbest, and red is Gbest

% Brian Birge
% Rev 2.0
% 10/30/03

clear all
close all
clc
help demopsobehavior
warning off

%aviobj = avifile(['demopsobehavior.avi']);
%aviobj.Quality=100;
%aviobj.Compression='Cinepak';
%aviobj.fps=15;

%xrng=input('Input search range for X, e.g. [-10,10] ? ');
%yrng=input('Input search range for Y ? ');
xrng=[-20,20];
yrng=[-20,20];

% if =0 then we look for minimum, =1 then max
 % minmax=input('Search for min (0) or max (1) ?');
  minmax=1;
  mvden = input('Max velocity divisor ? '); 
  ps    = input('How many particles? ');
  modl  = input('Choose PSO model? ');
 % note: if errgoal=NaN then unconstrained min or max is performed
  if minmax==1
    %  errgoal=0.97643183; % max for f6 function (close enough for termination)
      errgoal=NaN;
  else
     % errgoal=0; % min
      errgoal=NaN;
  end
  minx = xrng(1);
  maxx = xrng(2);
  miny = yrng(1);
  maxy = yrng(2);

%--------------------------------------------------------------------------
  dims=2;
  varrange=[];
  mv=[];
  for i=1:dims
      varrange=[varrange;minx maxx];
      mv=[mv;(varrange(i,2)-varrange(i,1))/mvden];
  end
  
  ac      = [2.1,2.1];% acceleration constants, only used for modl=0
  Iwt     = [0.9,0.6];  % intertia weights, only used for modl=0
  shw     = 100;   % how often to update display
  epoch   = 20000; % max iterations
  wt_end  = 1500; % iterations it takes to go from Iwt(1) to Iwt(2), only for modl=0
  errgrad = 1e-99;   % lowest error gradient tolerance
  errgraditer=5000; % max # of epochs without error change >= errgrad
  PSOseed = 0;    % if=1 then can input particle starting positions, if= 0 then all random
  % starting particle positions (first 20 at zero, just for an example)
   PSOseedValue = repmat([0],ps-10,1);
  
  psoparams=...
   [shw epoch ps ac(1) ac(2) Iwt(1) Iwt(2) ...
    wt_end errgrad errgraditer errgoal modl PSOseed];

 % run pso
 % vectorized version (normally instead of goplotpso4demo use goplotpso or your own)
  pso_out=pso_Trelea_vectorized('f6', dims,...
      mv, varrange, minmax, psoparams,'goplotpso4demo',PSOseedValue);
  %% vectorized version (normally instead of goplotpso4demo use goplotpso or your own)
  %[pso_out,tr,te]=pso_Trelea_vectorized('f6', dims,...
  %    mv, varrange, minmax, psoparams,'goplotpso',PSOseedValue);
  
  %% non-vectorized (much slower but some cost functions can't be vectorized)
 %  pso_out=pso_Trelea('f6', dims,...
 %      mv, varrange, minmax, psoparams,'goplotpso4demo');  
%-------------------------------------------------------------------------- 
disp('Best fit parameters: cost = f6(x,y)');
disp(['    x = ',num2str(pso_out(1))]);
disp(['    y = ',num2str(pso_out(2))]);
disp([' cost = ',num2str(pso_out(3))]);
disp(['mean(te) = ',num2str(mean(te))]);
%aviobj = close(aviobj);

indx=find(bestpos(:,3)<1e-4);
figure
plot(gbest_pred(:,1),gbest_pred(:,2),'g.','Markersize',5)
hold on
plot(bestpos(:,1),bestpos(:,2),'b.','markersize',5)
plot(bestpos(indx,1),bestpos(indx,2),'r.','markersize',8)
%clear all,for i=1:5000, offset(i)=linear_dyn(0.5);, pause(0.01),end,x=offset;,y=x;
%clear all,for i=1:2000, [x(i),y(i)]=spiral_dyn(1,.1);, pause(0.01),end
%plot(x,y,'k-'), grid on, axis equal
%title('red=good match gbest, blue=gbest, black is known')
hold off