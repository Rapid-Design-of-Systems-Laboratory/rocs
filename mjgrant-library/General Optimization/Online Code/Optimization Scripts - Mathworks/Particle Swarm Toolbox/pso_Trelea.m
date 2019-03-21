function [OUT,varargout]=pso_Trelea(functname,D,varargin)
% pso_Trelea.m
% a generic particle swarm optimizer
% to find the minimum or maximum of any 
% MISO matlab function
%
% Implements both of Trelea's types as well as generic inertia
%
% Usage:
%    [optOUT]=PSO(functname,D)
%   or:
%    [optOUT,tr,te]=PSO(functname,D,mv,VarRange,minmax,PSOparams,plotfcn)
%
% Inputs:
%    functname - string of matlab function to optimize
%    D - # of inputs to the function (dimension of problem)
%    
% Optional Inputs:
%    mv - max particle velocity, either a scalar or a vector of length D
%           (this allows each component to have it's own max velocity), 
%           default = 4, set if not input or input as NaN
%
%    VarRange - matrix of ranges for each input variable, 
%      default -100 to 100, of form:
%       [ min1 max1 
%         min2 max2
%            ...
%         minD maxD ]
%
%    minmax = 0, funct minimized (default)
%           = 1, funct maximized
%           = 2, funct is targeted to P(12) (minimizes distance to errgoal)
%
%    PSOparams - PSO parameters
%      P(1) - Epochs between updating display, default = 25. if 0, no display
%      P(2) - Maximum number of iterations (epochs) to train, default = 2000.
%      P(3) - population size, default = 20
%
%      P(4) - acceleration const 1 (local best influence), default = 2
%      P(5) - acceleration const 2 (global best influence), default = 2
%      P(6) - Initial inertia weight, default = 0.9
%      P(7) - Final inertia weight, default = 0.4
%      P(8) - Epoch when inertial weight at final value, default = 1500
%      P(9)- minimum global error gradient, 
%                 if abs(Gbest(i+1)-Gbest(i)) < gradient over 
%                 certain length of epochs, terminate run, default = 1e-9
%      P(10)- epochs before error gradient criterion terminates run, 
%                 default = 50, if the SSE does not change over 50 epochs
%                               then exit
%      P(11)- error goal, if NaN then unconstrained min or max, default=NaN
%      P(12)- type flag, 1,2, for Trelea types, 0 for common inertia type
%                 default = 0
%      
%    plotfcn - optional name of plotting function, default 'goplotpso',
%              make your own and put here
%
% Outputs:
%    optOUT - optimal inputs and associated min/max output of function, of form:
%        [ bestin1
%          bestin2
%            ...
%          bestinD
%          bestOUT ]
%
% Optional Outputs:
%    tr    - Gbest at every iteration, traces flight of swarm
%    te    - epochs to train, returned as a vector 1:endepoch
%
% Example:  out=pso_Trelea('f6',2)
%
% See Also: pso_Trelea_vectorized

% Brian Birge
% Rev 2.5
% 4/29/04

rand('state',sum(100*clock));
if nargin < 2
   error('Not enough arguments.');
end

% PSO PARAMETERS
if nargin == 2    % only specified functname and D
   VRmin=ones(D,1)*-100; 
   VRmax=ones(D,1)*100;    
   VR=[VRmin,VRmax];
   minmax = 0;
   P = [];
   mv = 4;
   plotfcn='goplotpso';   
elseif nargin == 3  % specified functname, D, and mv
   VRmin=ones(D,1)*-100; 
   VRmax=ones(D,1)*100;    
   VR=[VRmin,VRmax];
   minmax = 0;
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end
   P = [];
   plotfcn='goplotpso';   
elseif nargin == 4  % specified functname, D, mv, Varrange
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end
   VR=varargin{2}; 
   minmax = 0;
   P = [];
   plotfcn='goplotpso';   
elseif nargin == 5 % Functname, D, mv, Varrange, and minmax
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end    
   VR=varargin{2};
   minmax=varargin{3};
   P = [];
   plotfcn='goplotpso';
elseif nargin == 6 % Functname, D, mv, Varrange, minmax, and psoparams
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end    
   VR=varargin{2};
   minmax=varargin{3};
   P = varargin{4};  % psoparams
   plotfcn='goplotpso';   
elseif nargin == 7 % Functname, D, mv, Varrange, minmax, and psoparams, plotfcn
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end    
   VR=varargin{2};
   minmax=varargin{3};
   P = varargin{4};  % psoparams
   plotfcn = varargin{5};   
else    
   error('Wrong # of input arguments.');
end

% sets up default pso params
Pdef=[25 2000 20 2 2 0.9 0.4 1500 1e-9 50 NaN 0];
Plen=length(P);
P=[P,Pdef(Plen+1:end)];

df  = P(1);
me  = P(2);
ps  = P(3);
ac1 = P(4);
ac2 = P(5);
iw1 = P(6);
iw2 = P(7);
iwe = P(8);
ergrd=P(9);
ergrdep=P(10);
errgoal=P(11);
trelea=P(12);

if ((minmax==2) & isnan(errgoal))
    error('errgoal = NaN, minmax = 2, change these in input');
end

if (P(1))~=0
  plotflg=1;
else
  plotflg=0;
end

% take care of setting max velocity param here
if length(mv)==1
 velmaskmin = -mv*ones(1,D);
 velmaskmax = mv*ones(1,D);
elseif length(mv)==D
 velmaskmin = forcerow(-mv);
 velmaskmax = forcerow(mv);
else
 error('Max vel must be either a scalar or same length as prob dimension D');
end

% PLOTTING
 message = sprintf('PSO: %%g/%g iterations, GBest = %%g.\n',me);
 
% INITIALIZE INITIALIZE INITIALIZE INITIALIZE INITIALIZE INITIALIZE

% initialize population of particles and their velocities at time zero,
% format of pos= (particle#, dimension, epoch)
 % construct random population positions bounded by VR
  pos(1:ps,1:D,1)=normalize(rand([ps,D]),VR',1); 
 % construct initial random velocities between -mv,mv
  vel(1:ps,1:D,1)=normalize(rand([ps,D]),...
      [forcecol(-mv),forcecol(mv)]',1);

% initial pbest positions vals
 pbest=pos;
 
 evstrg1=char(ones(ps,1)*['feval(''',functname,'''']);
 evstrg3=char(ones(ps,1)*[')']);
 
 for j=1:ps  % start particle loop
    numin=sprintf(',%20.50f',pos(j,1:D));    
    evstrg=strcat('feval(''',functname,'''',char(numin(1:end)),')');
    out(j)=eval(evstrg);     % evaluate desired function with particle j      
 end
 
 pbestval=out;   % initially, pbest is same as pos

% assign initial gbest here also (gbest and gbestval)
 if minmax==1
   % this picks gbestval when we want to maximize the function  
    [gbestval,idx1]=max(pbestval);
 elseif minmax==0
   % this works for straight minimization  
    [gbestval,idx1]=min(pbestval);
 elseif minmax==2
   % this works when you know target but not direction you need to go
   % good for a cost function that returns distance to target that can be either
   % negative or positive (direction info)
    [temp,idx1]=min((pbestval-ones(size(pbestval))*errgoal).^2);
    gbestval=pbestval(idx1);
 end
     
 gbest=pbest(idx1,:);  % this is gbest position
 gbesthist=gbest;
 tr(1)=gbestval;       % save for output

% INITIALIZE END INITIALIZE END INITIALIZE END INITIALIZE END

% start PSO iterative procedures
cnt=0; % counter used for updating display according to df in the options
cnt2=0; % counter used for the stopping subroutine based on error convergence

for i=1:me  % start epoch loop (iterations)

   for j=1:ps  % start particle loop
     numin=sprintf(',%20.50f',pos(j,1:D));    
     evstrg=strcat('feval(''',functname,'''',char(numin(1:end)),')');    
     out(j)=eval(evstrg); % evaluate desired function with particle j  
     e(j) = out(j); % use to minimize or maximize function to unknown values     
     
    % update pbest to reflect whether searching for max or min of function
     if minmax==0
       if pbestval(j)>=e(j);
          pbestval(j)=e(j);
          pbest(j,:)=pos(j,:);
       end
     elseif minmax==1
       if pbestval(j)<=e(j);
           pbestval(j)=e(j);
           pbest(j,:)=pos(j,:);
       end
     elseif minmax==2
       if (e(j)-errgoal)^2 <= (pbestval(j)-errgoal)^2
           pbestval(j)=e(j);
           pbest(j,:)=pos(j,:);
       end
     end
           
    % assign gbest by finding minimum of all particle pbests 
     if minmax==1
      % this picks gbestval when we want to maximize the function   
       [iterbestval,idx1]=max(pbestval);  
       if gbestval<=iterbestval
          gbestval=iterbestval;
          gbest=pbest(idx1,:);
  %        gbesthist=[gbesthist;gbest];  % this line will slow algo down a lot
       end       
     elseif minmax==0 
      % this works for straight minimization and for minimizing error 
      % to target   
       [iterbestval,idx1]=min(pbestval);
       if gbestval>=iterbestval
          gbestval=iterbestval;
          gbest=pbest(idx1,:);
  %        gbesthist=[gbesthist;gbest];  % this line will slow algo down a lot
       end       
     elseif minmax==2
       [temp,idx1]=min((pbestval-ones(size(pbestval))*errgoal).^2);
       iterbestval=pbestval(idx1);
       if (iterbestval-errgoal)^2 <= (gbestval-errgoal)^2
          gbestval=iterbestval;
          gbest=pbest(idx1,:);
  %        gbesthist=[gbesthist;gbest];  % this line will slow algo down a lot      
       end
     end    
    
     tr(i+1)=gbestval; % keep track of global best val
     te=i;             % returns epoch number to calling program when done
     
     assignin('base','temp_pso_out',[gbest';gbestval]);
     assignin('base','temp_te',[0:te]);
     assignin('base','temp_tr',[tr]);     
%PSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSO

    % get new velocities, positions (this is the heart of the PSO algorithm)     
    % each epoch get new set of random numbers
     rannum1=rand([1,D]);
     rannum2=rand([1,D]);
     ac11=rannum1.*ac1;
     ac22=rannum2.*ac2;
        
     if trelea==2    
      % from Trelea's paper, parameter set 2
       vel(j,1:D)=0.729.*vel(j,1:D)...                        % prev vel
                  +1.494*rannum1.*(pbest(j,1:D)-pos(j,1:D))...% independent
                  +1.494*rannum2.*(gbest(1,1:D)-pos(j,1:D));  % social  
     elseif trelea==1
      % from Trelea's paper, parameter set 1            
       vel(j,1:D)=0.6.*vel(j,1:D)...                          % prev vel
                  +1.7*rannum1.*(pbest(j,1:D)-pos(j,1:D))...  % independent
                  +1.7*rannum2.*(gbest(1,1:D)-pos(j,1:D));    % social              
     else
      % common PSO algo with inertia wt 
      % get inertia weight, just a linear funct w.r.t. epoch parameter iwe
       if i<=iwe
          iwt(i)=((iw2-iw1)/(iwe-1))*(i-1)+iw1; 
       else
          iwt(i)=iw2;
       end
       
       vel(j,1:D)=iwt(i).*vel(j,1:D)...                % prev vel
                  +ac11.*(pbest(j,1:D)-pos(j,1:D))...  % independent
                  +ac22.*(gbest(1,1:D)-pos(j,1:D));     % social                  
     end
            
    % limit velocities here using masking
     minvelmask_throwaway= vel(j,:) <= velmaskmin;       
     minvelmask_keep= vel(j,:) >  velmaskmin;
      
     newvelA=vel(j,:).*minvelmask_keep; % keeps good vals, zeros bad
     newvelB=velmaskmin.*minvelmask_throwaway; % 
     vel(j,:)=newvelA+newvelB;  % takes care of vals < -maxvel
       
     maxvelmask_throwaway= vel(j,:) >= velmaskmax;
     maxvelmask_keep= vel(j,:) <  velmaskmax;
       
     newvelA=vel(j,:).*maxvelmask_keep;
     newvelB=velmaskmax.*maxvelmask_throwaway;
     vel(j,:)=newvelA+newvelB; % takes care of vals > maxvel      
        
    % update new position (PSO algo)    
     pos(j,:)=pos(j,:)+vel(j,:);
     
    % position masking, limits positions to desired search space
    % method: 0) no position limiting, 1) saturation at limit,
    %         2) wraparound at limit , 3) bounce off limit
     posmaskmeth=3;     
     
     posmaskmin=VR(:,1)';
     posmaskmax=VR(:,2)';
     
     minposmask_throwaway = pos(j,:) <=posmaskmin;
     minposmask_keep      = pos(j,:) > posmaskmin;     
     maxposmask_throwaway = pos(j,:) >=posmaskmax;
     maxposmask_keep      = pos(j,:) < posmaskmax;
     
    if posmaskmeth==1
     % this is the saturation method (Eberhart inertia method)
      newposA=pos(j,:).*minposmask_keep;
      newposB=posmaskmin.*minposmask_throwaway;
      pos(j,:)=newposA+newposB;
      
      newposA=pos(j,:).*maxposmask_keep;
      newposB=posmaskmax.*maxposmask_throwaway;
      pos(j,:)=newposA+newposB;
      
      %vel(j,:)=zeros(size(vel(j,:)));
     elseif posmaskmeth==2
     % this is the wraparound method    
      newposA=pos(j,:).*minposmask_keep;
      newposB=posmaskmin.*maxposmask_throwaway;
      pos(j,:)=newposA+newposB;
      
      newposA=pos(j,:).*maxposmask_keep;
      newposB=posmaskmax.*minposmask_throwaway;
      pos(j,:)=newposA+newposB;      
     elseif posmaskmeth==3
     % this is the bounce method, particles bounce off the boundaries with -vel      
      newposA=pos(j,:).*minposmask_keep;
      newposB=posmaskmin.*minposmask_throwaway;
      pos(j,:)=newposA+newposB;
      vel(j,:)=vel(j,:).*minposmask_keep + -vel(j,:).*minposmask_throwaway;
      
      newposA=pos(j,:).*maxposmask_keep;
      newposB=posmaskmax.*maxposmask_throwaway;
      pos(j,:)=newposA+newposB;
      vel(j,:)=vel(j,:).*maxposmask_keep + -vel(j,:).*maxposmask_throwaway;
    end
    
%PSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSOPSO     
       
   end         % end particle loop
  %------------------------------------------------------------------------      
   % check for stopping criterion based on speed of convergence to desired 
   % error   
    tmp1=abs(tr(i)-gbestval);
    if tmp1>ergrd
       cnt2=0;
    elseif tmp1<=ergrd
       cnt2=cnt2+1;
       if cnt2>=ergrdep
         if plotflg==1
          fprintf(message,i,gbestval);           
          disp(' ');
          disp(['--> Solution likely, GBest hasn''t changed for ',...
                  num2str(cnt2),' epochs.']);
          eval(plotfcn);
         end       
         break
       end
    end
   % this stops if using constrained optimization and goal is reached
   if ~isnan(errgoal)
    if ((gbestval<=errgoal) & (minmax==0)) | ((gbestval>=errgoal) & (minmax==1))

        if plotflg==1
            fprintf(message,i,gbestval);
          disp(' ');            
            disp(['--> Error Goal reached, successful termination!']);            
            eval(plotfcn);
        end
        break
    end
   % this is stopping criterion for constrained from both sides    
    if minmax==2
      if ((tr(i)<errgoal) & (gbestval>=errgoal)) | ((tr(i)>errgoal) ...
              & (gbestval<=errgoal))        
        if plotflg==1
            fprintf(message,i,gbestval);
          disp(' ');            
            disp(['--> Error Goal reached, successful termination!']);            
            eval(plotfcn);
        end
        break              
      end
    end % end if minmax==2
   end  % end ~isnan if
   
% this section does the plots during iterations   
if plotflg==1      
  if (rem(i,df) == 0 ) | (i==me) | (i==1) 
     fprintf(message,i,gbestval);
     cnt=cnt+1; % count how many times we display (useful for movies)
       
       eval(plotfcn); % defined at top of script
       
  end  % end update display every df if statement    
 end % end plotflg if statement
end  % end epoch loop

% clear temp outputs
 evalin('base','clear temp_pso_out temp_te temp_tr;');
% hold off
% outputs
 OUT=[gbest';gbestval];
 varargout{1}=[0:te];
 varargout{2}=[tr];