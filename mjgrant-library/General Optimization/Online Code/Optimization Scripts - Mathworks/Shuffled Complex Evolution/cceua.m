%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [snew,fnew,icall]=cceua(s,sf,bl,bu,icall,maxn)
%  This is the subroutine for generating a new point in a simplex
%
%   s(.,.) = the sorted simplex in order of increasing function values
%   s(.) = function values in increasing order
%
% LIST OF LOCAL VARIABLES
%   sb(.) = the best point of the simplex
%   sw(.) = the worst point of the simplex
%   w2(.) = the second worst point of the simplex
%   fw = function value of the worst point
%   ce(.) = the centroid of the simplex excluding wo
%   snew(.) = new point generated from the simplex
%   iviol = flag indicating if constraints are violated
%         = 1 , yes
%         = 0 , no

[nps,nopt]=size(s);
n = nps;
m = nopt;
alpha = 1.0;
beta = 0.5;

% Assign the best and worst points:
sb=s(1,:); fb=sf(1);
sw=s(n,:); fw=sf(n);

% Compute the centroid of the simplex excluding the worst point:
ce=mean(s(1:n-1,:));

% Attempt a reflection point
snew = ce + alpha*(ce-sw);

% Check if is outside the bounds:
ibound=0;
s1=snew-bl; idx=find(s1<0); if ~isempty(idx); ibound=1; end;
s1=bu-snew; idx=find(s1<0); if ~isempty(idx); ibound=2; end;

if ibound >=1; 
    snew = bl + rand(1,nopt).*(bu-bl);
end;
fnew = functn(nopt,snew);
icall = icall + 1;

% Reflection failed; now attempt a contraction point:
if fnew > fw;
    snew = sw + beta*(ce-sw);
    fnew = functn(nopt,snew);
    icall = icall + 1;

% Both reflection and contraction have failed, attempt a random point;
    if fnew > fw;
        snew = bl + rand(1,nopt).*(bu-bl);
        fnew = functn(nopt,snew);
        icall = icall + 1;
    end;
end;

% END OF CCE
return;
