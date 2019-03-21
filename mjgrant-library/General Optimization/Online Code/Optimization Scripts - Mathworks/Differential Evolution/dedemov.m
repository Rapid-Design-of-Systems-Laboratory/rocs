function dedemov(action);
% DE function minimization demonstration.
% dedemov is called with no parameters.
%
% Differential Evolution for MATLAB
% Copyright (C) 1996 R. Storn
% International Computer Science Institute (ICSI)
% 1947 Center Street, Suite 600
% Berkeley, CA 94704
% E-mail: storn@icsi.berkeley.edu
% WWW:    http://http.icsi.berkeley.edu/~storn
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 1, or (at your option)
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. A copy of the GNU 
% General Public License can be obtained from the 
% Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 	

if nargin<1,   %---number of arguments less than one ?
    action='initialize';
end;

if strcmp(action,'initialize'),
%---Header of Window--------------------------------
    figNumber=figure( ...
	'Name','Differential Evolution Demo', ...
        'NumberTitle','off', ...
	'Visible','off');

%---define area for graphics------------------------
    axes( ...
	'Units','normalized', ...
	'Visible','off', ...
        'Position',[0.1 0.35 0.60 0.75]);

%---Set up the Comment Window-----------------------
    top1=0.25;
    left1=0.05;
    right1=0.75;
    bottom1=0.05;
    labelHt1=0.05;
    spacing1=0.005;
    promptStr=str2mat(' ',...
               ' This demo shows the track of the best', ...
	       ' population member for each generation.');

%---First, the MiniCommand Window frame------------------------
    frmBorder=0.02;
    frmPos1=[left1-frmBorder bottom1-frmBorder ...
	(right1-left1)+2*frmBorder (top1-bottom1)+2*frmBorder];
    uicontrol( ...
        'Style','frame', ...
        'Units','normalized', ...
        'Position',frmPos1, ...
	'BackgroundColor',[0.5 0.5 0.5]);
%---Then the text label----------------------------------------
    labelPos1=[left1 top1-labelHt1 (right1-left1) labelHt1];
    uicontrol( ...
	'Style','text', ...
        'Units','normalized', ...
        'Position',labelPos1, ...
        'BackgroundColor',[0.5 0.5 0.5], ...
	'ForegroundColor',[1 1 1], ...
        'String','Comment Window');
%---Then the editable text field-------------------------------
    txtPos1=[left1 bottom1 (right1-left1) top1-bottom1-labelHt1-spacing1];
    txtHndl=uicontrol( ...
   	'Style','edit', ...
        'Units','normalized', ...
        'Max',10, ...
        'BackgroundColor',[1 1 1], ...
        'Position',txtPos1, ...
   	'Callback','dedemov(''eval'')', ...
        'String',promptStr);

%---Save the text handle for future use------------------------
    set(gcf,'UserData',txtHndl)


%==============================================================
%---Information for all buttons (and menus)--------------------
    labelColor=[0.8 0.8 0.8];
    yInitPos=0.90;
    menutop=0.95;
    btnTop = 0.6;
    top=0.75;
    left=0.785;
    btnWid=0.175;
    btnHt=0.06;
    textHeight = 0.05;
    textWidth = 0.06;
    %---Spacing between the button and the next command's label
    spacing=0.019;

%---The CONSOLE frame------------------------------------------
    frmBorder=0.019; frmBottom=0.04; 
    frmHeight = 0.92; frmWidth = btnWid;
    yPos=frmBottom-frmBorder;
    frmPos=[left-frmBorder yPos frmWidth+2*frmBorder frmHeight+2*frmBorder];
    h=uicontrol( ...
        'Style','frame', ...
        'Units','normalized', ...
        'Position',frmPos, ...
	'BackgroundColor',[0.5 0.5 0.5]);

%===Define the individual selection items======================

%---The DE strategy Selection Menu-----------------------------
    btnNumber=1;
    yPos=menutop-(btnNumber-1)*(btnHt+spacing);
    btnPos=[left yPos-btnHt btnWid btnHt];
    labelStr='rand/1|best/1|best/2|rtb/1|rand/2';
    callbackStr='dedemov(''changemethod'');';
    MethodHndl=uicontrol( ...
        'Style','popupmenu', ...
        'Units','normalized', ...
        'Position',btnPos, ...
        'String',labelStr, ...
        'Interruptible','yes', ...
        'Callback',callbackStr);

    
%---Size of Population NP--------------------------------------
    btnNumber=1;
    yPos=menutop-(btnNumber-1)*(btnHt+spacing);
    top = yPos - btnHt - spacing;
    labelWidth = frmWidth-textWidth-.01;
    labelBottom=top-textHeight;
    labelLeft = left;
    labelPos = [labelLeft labelBottom labelWidth textHeight];
    h = uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
	'Position',labelPos, ...
        'Horiz','left', ...
	'String','NP', ...
        'Interruptible','no', ...
        'BackgroundColor',[0.5 0.5 0.5], ...
        'ForegroundColor','white');


%---Weighting factor F-----------------------------------------
    btnNumber=2;
    yPos=menutop-(btnNumber-1)*(btnHt+spacing);
    top = yPos - btnHt - spacing;
    labelWidth = frmWidth-textWidth-.01;
    labelBottom=top-textHeight;
    labelLeft = left;
    labelPos = [labelLeft labelBottom labelWidth textHeight];
    h = uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
        'Position',labelPos, ...
        'Horiz','left', ...
	'String','F', ...
        'Interruptible','no', ...
	'BackgroundColor',[0.5 0.5 0.5], ...
	'ForegroundColor','white');


%---Crossover probability CR-----------------------------------
    btnNumber=3;
    yPos=menutop-(btnNumber-1)*(btnHt+spacing);
    top = yPos - btnHt - spacing;
    labelWidth = frmWidth-textWidth-.01;
    labelBottom=top-textHeight;
    labelLeft = left;
    labelPos = [labelLeft labelBottom labelWidth textHeight];
    h = uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
	'Position',labelPos, ...
        'Horiz','left', ...
	'String','CR', ...
        'Interruptible','no', ...
        'BackgroundColor',[0.5 0.5 0.5], ...
	'ForegroundColor','white');

%---Maximum number of iterations itermax------------------------
    btnNumber=4;
    yPos=menutop-(btnNumber-1)*(btnHt+spacing);
    top = yPos - btnHt - spacing;
    labelWidth = frmWidth-textWidth-.01;
    labelBottom=top-textHeight;
    labelLeft = left;
    labelPos = [labelLeft labelBottom labelWidth textHeight];
    h = uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
	'Position',labelPos, ...
        'Horiz','left', ...
	'String','itermax', ...
        'Interruptible','no', ...
        'BackgroundColor',[0.5 0.5 0.5], ...
	'ForegroundColor','white');

%---Text field---------------------------------------------------
    textPos = [labelLeft+labelWidth-.015 labelBottom textWidth+.025 7*.85*textHeight];
    callbackStr = 'dedemov(''setContr'')';

%---Enter the default values-------------------------------------
    str = sprintf('15\n\n0.9\n\n0.9\n\n200');
    mat = [15; 0.9; 0.9; 200]; % default values
    CvarsHndl = uicontrol( ...
        'Style','edit', ...
        'Units','normalized', ...
        'Position',textPos, ...
        'Max',2, ... % makes this a multiline edit object
        'Horiz','right', ...
        'Background','white', ...
        'Foreground','black', ...
        'String',str,'Userdata',mat, ...
        'callback',callbackStr);

%---The INFO button-----------------------------------------------
    uicontrol( ...
	'Style','push', ...
        'Units','normalized', ...
        'Position',[left bottom1+btnHt+spacing btnWid btnHt], ...
        'String','Info', ...
        'Callback','dedemov(''info'')');

%---The CLOSE button--------------------------------------------
    uicontrol( ...
	'Style','push', ...
        'Units','normalized', ...
        'Position',[left bottom1 btnWid btnHt], ...
        'String','Close', ...
        'Callback','close(gcf)');

%---Get Handles--------------------------------------------------
    chndlList=[MethodHndl CvarsHndl txtHndl];

    %---Now uncover the figure-----------------------------------
    set(figNumber, ...
	'Visible','on', ...
	'UserData',chndlList);

    dedemov('demobutton') %---> Go and optimize
    return

%===Different actions============================================

%---Modified control variables-----------------------------------
elseif strcmp(action,'setContr'),
    v = get(gco,'Userdata');
    str = get(gco,'String');
    ind = find(abs(str)<32);
    str(ind) = 32*ones(size(ind));
    str = str';
    vv = eval(['[' str(:)' ']'],'-1')';
    if vv == -1
        vv = v;
    elseif  length(vv)~=4  %must be 4 items as specified above
        vv = v;
    end

    str = sprintf('%g\n\n%g\n\n%g\n\n%g\n\n%g\n\n%g',vv(1),vv(2),vv(3),vv(4));%,...
          %vv(5), vv(6));
    set(gco,'Userdata',vv,'String',str)  %take new data

    dedemov('demobutton') %---> Go and optimize
    return

%---Modified method----------------------------------------------
elseif strcmp(action,'changemethod'),
    v = get(gco,'value');  % 1 = DE/rand/1, 2 = DE/best/1, 3 = DE/best/2,
			   % 4 = DE/rand-to-best/1, 5 = DE/rand/2
    chndlList = get(gcf,'UserData');

    dedemov('demobutton') %---> Go and optimize
    return

%-----Start optimization-----------------------------------
elseif strcmp(action,'demobutton'),

%---Plotting the banana function---------------------------
    cla reset;
    %axis off; 
    x=-2:.2:2;
    y=-1:.2:3;
    [xx,yy]=meshgrid(x,y);
    zz=100*(yy-xx.^2).^2+(1-xx).^2;

    %---Set up the appropriate colormap--------------------------------------
    %---In this case, the colormap has been chosen to give the surface plot--
    %---a nice healthy banana color.-----------------------------------------
    hsv2=hsv;
    hsv3=[hsv2(11:64,:); hsv2(1:10,:)];

    %---draw the surface plot-----------------------
    surfHndl=surface(x,y,zz,'EdgeColor',[.8 .8 .8]);
    axis([-2 2 -2 2]);
    text(-1,-3,'x(1)'); %---substitute for xlabel
    ylabel('x(2)');
    %axis off;
    view(10,55);
    colormap(hsv3);
    hold on;
    [c,contHndl]=contour3(x,y,zz+50,[100 500],'k');
    set(contHndl,'Color',[.8 .8 .8]);
    drawnow
    plot3(1,1,0,'ko', ...
    	'MarkerSize',15, ...
    	'LineWidth',2, ...
    	'EraseMode','none');
    text(0.8,1.4,0,'   goal', ...
    	'Color',[0 0 0], ...
    	'EraseMode','none');

%---Get handles----------------------------------------

    chndlList=get(gcf,'Userdata');
    MethodHndl = chndlList(1);
    CvarsHndl  = chndlList(2);
    txtHndl    = chndlList(3);

%---initialize variables-------------------------------
    cvars = get(CvarsHndl,'UserData');
    NP = cvars(1);                   % population size
    F  = cvars(2);                   % weighting factor
    CR = cvars(3);                   % crossover probability
    itermax = cvars(4);              % maximum number of iterations
    strategy = get(MethodHndl,'value');
    D = 2;

    %---as a reminder labelStr='rand/1|best/1|best/2|rtb/1';

%---Now run the DE demo eventually---------------------

    [bestmem,nfeval]=devec(NP,D,F,CR,itermax,strategy);
    str1=sprintf(' x(1) = %d   x(2) = %d',bestmem(1),bestmem(2));
    str2=sprintf(' Number of function evaluations: %d', nfeval);
    %str1=get(txtHndl,'String');
    strout=str2mat(str1,' ',str2);
    set(txtHndl,'String',strout);

%---Info button----------------------------------------
elseif strcmp(action,'info'),
    ttlStr=get(gcf,'Name');
    hlpStr= ...                                               
        ['                                                 '  
         ' This demonstration shows the minimization       '  
	 ' of Rosenbrock''s "banana function":              '
         '                                                 '  
         '     f(x)= 100*(x(2)-x(1)^2)^2+(1-x(1))^2        '  
         '                                                 '  
         ' by Differential Evolution (DE). DE is a sto-    '  
         ' chastic minimization procedure for continuous   '  
         ' space functions that may be non-differentiable, '  
         ' nonlinear and multimodal. DE requires just 3    '  
         ' control variables:                              '  
         '                                                 '  
         '   Number of population members NP               '  
         '   Difference vector weight     F   ex [0, 2]    '  
         '   Crossover probability        CR  ex [0, 1]    '  
	 '                                                 '
	 ' As a first guess NP=10*(number of parameters)   '
	 ' is usually a good choice.                       '
	 '                                                 '
         ' File name: dedemov.m                            '];
    helpfun(ttlStr,hlpStr);                                   

end;    % if strcmp(action, ...
