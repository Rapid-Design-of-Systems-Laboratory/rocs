%----------------------------------------------%
% Begin Function:  minimumTimeToClimbCompute.m %
%----------------------------------------------%
function [CD,CL,eta]=minimumTimeToClimbCompute(Mach,CONSTANTS)

CDdat  = CONSTANTS.CDdat;
CLdat  = CONSTANTS.CLdat;
etadat = CONSTANTS.etadat;

%ii     = find(Mach>=0.8);
%jj     = find(Mach<0.8);

ii     = Mach>=0.8;
jj     = Mach<0.8;

mpoly  = Mach(ii);
CD  = zeros(length(Mach),1);
CL = zeros(length(Mach),1);
eta = zeros(length(Mach),1);
if any(ii)
    CD(ii)  = ppval(CDdat,mpoly);
    CL(ii)  = ppval(CLdat,mpoly);
    eta(ii) = ppval(etadat,mpoly);
end

if any(jj)
    CD(jj)  = 0.013;
    CL(jj)  = 3.44;
    eta(jj) = 0.54;
end

%----------------------------------------------%
% End Function:  minimumTimeToClimbCompute.m   %
%----------------------------------------------%
