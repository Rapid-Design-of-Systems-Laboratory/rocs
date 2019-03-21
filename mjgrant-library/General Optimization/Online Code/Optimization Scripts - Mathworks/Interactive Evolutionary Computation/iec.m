function [popu] = iec(funname,popuin,nindin,ndimin);
%IEC - Interactive Evolutionary Computation (interactive ES)
%
%  [popu] = iec(funname,popuin,nind,ndim)
%    popu:     resulted population variable (structure)
%    funname:  user's function (string)
%    popuin:   initial population (structure), if popu=[] -> default init.
%    nind:     number of individuals (integer) if popu=[] else not used
%    ndim:     number of design var.s (integer) if popu=[] else not used
%
%    popu.LB:     low boundaries (real 1 x ndim)
%    popu.UB:     upper boundaries (real 1 x ndim)
%    popu.chrom:  chromosomes (real nind x ndim)
%    popu.devi:   strategy variables (real nind x ndim)
%    popu.niter:  number of iterations (integer)
%    popu.nfeval: number of fun. evaluations (integer)
%
% EXAMPLE:
% Calling the iec.m function:
%
%   popusize = 7; %size of population
%   ndim = 3;     %number of search space dimensions
%   popu = iec('testuserfun',[],7,3);
%
% The testuserfun.m function:
%
%   function [numplts] = testuserfun(x,plts);
%   % x:    chromosome (it contains the 3 variables)
%   % plts: vector of subplots of interactive figure (the 2 subplots)
%   %Return the number of plots for interactive figure, if x = []
%   if isempty(x),
%     numplts = 2; %we need two plots for every potential solution
%     return;
%   end
%   %Else chromosome (x) evaluation and draw plots for the user
%   % (it's only a simple example whithout any meaning)
%   % 1. plot
%   subplot(plts(1));
%   t = [0.01:0.01:1];
%   y = sin(x(1)*t+x(2))+x(3)*sin(x(1)*t);
%   plot(t,y,'r');
%   % 2. plot
%   subplot(plts(2));
%   t = [0.01:0.01:1];
%   y = cos(x(3)*t+x(1))+x(2)*sin(x(1)*t);
%   plot(t,y,'b');
%


% (c) Janos Madar (Madár János), University of Veszprem, Hungary, 2005

%---------------------------- Variables ----------------------------

%Population variable
popu = popuin;

if isempty(popu),
  nind = nindin;
  ndim = ndimin;
  popu.LB = zeros(1,ndim);
  popu.UB = ones(1,ndim);
  popu.chrom = rand(nind,ndim);
  popu.devi = rand(nind,ndim)*1/sqrt(ndim);
  popu.niter = 0;
  popu.nfeval = 0;
else
  nind = size(popu.chrom,1);
  ndim = size(popu.chrom,2);
end

%ES parameters
tau0 = 1/sqrt(2*ndim);
taui = 1/sqrt(2*sqrt(ndim));

%------------------------------ IEC-ES ------------------------------

ifig = [];
nplts = feval(funname,[]); %get number of subplots
eltv = [];
while 1,

  %Evaluation (interactive)
  exite = 3;
  while exite == 3,
    [exite,popu,ifig,parv,delv] = iec_iface(popu,funname,ifig,nplts,length(eltv));
    eltv = setdiff([1:nind],delv);
  end
  if exite == 1,
    break;
  end

  %Create new generation
  chrom = popu.chrom;
  devi = popu.devi;
  %reordered and insert elites
  npar = length(parv);
  dummy = [];
  for k = 1:npar,
    if ~isempty(find(eltv==parv(k))),
      dummy = [dummy parv(k)];
      eltv(find(eltv==parv(k))) = [];
    end
  end
  dummy = [dummy eltv];
  eltv = dummy;
  k = length(eltv);
  chrom(1:k,:) = popu.chrom(eltv,:);
  %reorder subplots
  ifig = iec_plotord(ifig,eltv);
  %new offsprings
  while k < nind,
    k = k+1;
    %choose parents
    ix1 = parv(floor(rand*npar)+1);
    ix2 = parv(floor(rand*npar)+1);
    %discrete recombination of design variables
    chrom(k,:) = popu.chrom(ix1,:);
    dummy = find(round(rand(1,ndim))==0);
    chrom(k,dummy) = popu.chrom(ix2,dummy);
    %intermediate recombination of strategy variables
    devi(k,:) = mean([popu.devi(ix1,:);popu.devi(ix2,:)]);
    %mutation of strategy variable
    z0 = tau0*randn;
    zi = taui*randn(1,ndim);
    devi(k,:) = devi(k,:).*exp(z0+zi);
    %mutation of design variables
    chrom(k,:) = chrom(k,:) + randn(1,ndim).*devi(k,:);
    %constraints
    if (~isempty(popu.LB)),
      chrom(k,:) = max(chrom(k,:),popu.LB);
    end
    if (~isempty(popu.UB)),
      chrom(k,:) = min(chrom(k,:),popu.UB);
    end
  end
  %the new population variable
  popu.chrom = chrom;
  popu.devi = devi;

  %Next iteration
  popu.niter = popu.niter + 1;

end



%======================================================================
%======================================================================

function [exite,popu,ifig,parv,delv] = iec_iface(popuin,fcnname,ifigin,nplts,nelt);
% Interavtive ES optimization - interface
%

%------------------------------ Begin ------------------------------

%Inic output and variables
exite = 0;
popu = popuin;
ifig = ifigin;
parv = [];
delv = [];
nind = size(popu.chrom,1);

%------------------------------ Figure ------------------------------

global buttonclick

%If the figure does not exist, create it
if isempty(ifig),
  %figure
  ifig.handle = figure('Color',[0.8 0.8 0.8]);
  set(ifig.handle,'DefaultAxesFontName','times');
  set(ifig.handle,'DefaultTextFontName','times');
  set(ifig.handle,'DefaultAxesFontSize',8);
  set(ifig.handle,'DefaultTextFontSize',8);
  %subplots
  for i = 1:nind,
    for j = 1:nplts,
      ifig.axs(i,j) = axes('Position',[2 2 0.2 0.2]);
    end
  end
  %header with line and labels
  ifig.axtop = axes('Position',[0 0.85 1 0.1]);
  axis off
  ifig.linetop = line([0;1],[0;0]);
  set(ifig.linetop,'Color',[0 0 0],'LineWidth',2);
  ifig.label1 = text(0,1.5,'Iteration: 0');
  ifig.label2 = text(0.2,1.5,'---');
  %selection buttons
  for i = 1:nind,
    ifig.btn1(i) = uicontrol('Units','normalized','Position',[2,2,0.1,0.1]);
    set(ifig.btn1(i),'String','Select');
    ss = sprintf('global buttonclick; buttonclick=%i; uiresume;',i+100);
    set(ifig.btn1(i),'CallBack',ss);
  end
  %delete buttons
  for i = 1:nind,
    ifig.btn2(i) = uicontrol('Units','normalized','Position',[2,2,0.1,0.1]);
    set(ifig.btn2(i),'String','Delete');
    ss = sprintf('global buttonclick; buttonclick=%i; uiresume;',i+200);
    set(ifig.btn2(i),'CallBack',ss);
  end
  %main control buttons
  ifig.btncnc = uicontrol('Units','normalized','Position',[0.5,0.975,0.15,0.03]);
  set(ifig.btncnc,'String','Cancel');
  set(ifig.btncnc,'CallBack','global buttonclick; buttonclick=1; uiresume;');
  ifig.btnok = uicontrol('Units','normalized','Position',[0.5,0.94,0.15,0.03]);
  set(ifig.btnok,'String','OK (next generation)');
  set(ifig.btnok,'CallBack','global buttonclick; buttonclick=2; uiresume;');
  ifig.btnexit = uicontrol('Units','normalized','Position',[0.85,0.96,0.1,0.03]);
  set(ifig.btnexit,'String','Exit');
  set(ifig.btnexit,'CallBack','global buttonclick; buttonclick=3; uiresume;');
  ifig.btnld = uicontrol('Units','normalized','Position',[0.7,0.96,0.1,0.03]);
  set(ifig.btnld,'String','Load');
  set(ifig.btnld,'CallBack','global buttonclick; buttonclick=4; uiresume;');
  ifig.btnsv = uicontrol('Units','normalized','Position',[0.7,0.92,0.1,0.03]);
  set(ifig.btnsv,'String','Save');
  set(ifig.btnsv,'CallBack','global buttonclick; buttonclick=5; uiresume;');
end

%Arrange select/delete buttons
for i = 1:nind,
  posw = 0.97/nind;
  posx = posw*(i-1)+0.1*posw+0.03;
  posw = posw*0.4;
  posy = 0.88;
  posh = 0.03;
  set(ifig.btn1(i),'Position',[posx,posy,posw,posh]);
  set(ifig.btn2(i),'Position',[posx+posw,posy,posw,posh]);
end
%arrange and reset subplots
for i = 1:nind,
  for j = 1:nplts,
    posw = 0.97/nind;
    posh = 0.9/nplts;
    posx = posw*(i-1)+0.1*posw+0.03;
    posy = 0.9-posh*j+0.1*posh;
    posw = posw*0.8;
    posh = posh*0.8;
    set(ifig.axs(i,j),'Position',[posx posy posw posh]);
    if i > nelt, 
      cla(ifig.axs(i,j));
    end
  end
end

%------------------------ Evaluation & User -------------------------

%Evaluate individuals
for i = 1:nind,
  if i > nelt,
    feval(fcnname,popu.chrom(i,:),ifig.axs(i,:));
    popu.nfeval = popu.nfeval + 1;
  end
end

%Some info
ss = sprintf('Iter:%i, feval:%i',popu.niter,popu.nfeval);
set(ifig.label1,'String',ss);

%Event handler
exite = 0;
parv = [];
delv = [];
while exite == 0;
  %label2 - selected parent pairs, and deleted individuals
  ss = 'Parents: ';
  for k = 1:length(parv),
    ss = strcat(ss,sprintf('%i, ',parv(k)));
  end
  ss = strcat(ss,'Delete: ');
  for k = 1:length(delv),
    ss = strcat(ss,sprintf('%i, ',delv(k)));
  end
  set(ifig.label2,'String',ss);
  %wait for button onclick (or window onclose) event
  buttonclick = 0;
  set(ifig.handle, 'vis', 'on', 'waitstatus', 'waiting');
  waitfor(ifig.handle, 'waitstatus', 'inactive');
  %if closed or <exit>
  if buttonclick == 0 | buttonclick == 3,
    exite = 1;
    break;
  end
  %<cancel>
  if buttonclick == 1,
    parv = [];
    delv = [];
    continue;
  end
  %<ok>
  if buttonclick == 2,
    if length(parv) > 0 & length(delv) > 0,
      exite = 2;
      break;
    end
  end
  %<load> 
  if buttonclick == 4,
    clear popu
    load popusave
    parv = [];
    delv = [1:nind];
    exite = 3;
    break;
  end
  %<save> 
  if buttonclick == 5,
    save popusave popu
    continue;
  end
  %i-th <select>
  if buttonclick > 100 & buttonclick < 200,
    if isempty(find(parv==buttonclick-100)),
      parv = [parv, buttonclick-100];
      parv = sort(parv);
    end
    continue;
  end
  %i-th <delete>
  if buttonclick > 200 & buttonclick < 300,
    if isempty(find(delv==buttonclick-200)),
      delv = [delv, buttonclick-200];
      delv = sort(delv);
    end
    continue;
  end
end

%======================================================================
%======================================================================
function [ifig] = iec_plotord(ifigin,ixv);
%Interactive figure subplot reorder

ifig = ifigin;
n = size(ifig.axs,1);
m = length(ixv);
ixv2 = setdiff([1:n],ixv);
ifig.axs(1:m,:) = ifigin.axs(ixv,:);
ifig.axs(m+1:n,:) = ifigin.axs(ixv2,:);
