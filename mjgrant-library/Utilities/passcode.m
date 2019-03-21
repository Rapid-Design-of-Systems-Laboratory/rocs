function CurrWord = passcode(action);
%PASSCODE  password input dialog box.
%  answer = PASSCODE creates a modal dialog box that returns user 
%  password input. Given characters are substituted with '*'-Signs 
%  like in usual Windows dialogs.
%
%  Allowed input characters can be predefined within the 
%  main function (see 'allowed_charakter'-subfunction at the end of the script).
%
%  See also PCODE.

% Version: v1.2 (22-Jul-2005)
% Author:  Elmar.Tarajan@Mathworks.de
% Modified: Michael J. Grant / SSDL GATech

if nargin == 0  % LAUNCH GUI
   ScreenSize = get(0,'ScreenSize');
   h = figure('Menubar','none', ...
              'Units','Pixels', ...
              'Resize','off', ...
              'NumberTitle','off', ...
              'CloseRequestFcn','set(findobj(''Tag'',''password''),''UserData'',[]);uiresume', ...
              'Name',['password required'], ...
              'KeyPressFcn', 'passcode(''KeyPress_Callback'');', ...       
              'Position',[ (ScreenSize(3:4)-[300 75])/2 300 75 ], ...
              'Color',[0.8 0.8 0.8], ...
              'WindowStyle','modal');
   uicontrol( 'Parent',h, ...
              'Style','Edit', ...
              'Enable','inactive', ...
              'Units','Pixels','Position',[49 28 202 22], ...
              'FontSize',15, ...
              'BackGroundColor',[0.7 0.7 0.7]);
   uicontrol( 'Parent',h, ...
              'Style','Text', ...
              'Tag','password', ...
              'Units','Pixels','Position',[51 30 198 18], ...
              'FontSize',15, ...
              'BackGroundColor',[1 1 1]);
   uicontrol( 'Parent',h, ...
              'Style','Text', ...
              'Tag','error', ...
              'Units','Pixels','Position',[50 2 200 20], ...
              'FontSize',8, ...
              'String','character not allowed',...
              'Visible','off',...
              'ForeGroundColor',[1 0 0], ...              
              'BackGroundColor',[0.8 0.8 0.8]);
   uiwait
   CurrWord = get(findobj('Tag','password'),'UserData');      
   delete(h)
else
   switch action
   case 'KeyPress_Callback'
      CurrChar = get(gcf,'CurrentCharacter');
      CurrWord = get(findobj('Tag','password'),'UserData');
      %
      if int8(CurrChar) == 8
          CurrWord = CurrWord(1:end-1);
      elseif int8(CurrChar) == 13
         uiresume
         return
      elseif ~isempty(CurrChar)
         if any(allowed_characters==CurrChar)
            CurrWord = [CurrWord CurrChar];
          else
            set(findobj('Tag','error'),'Visible','on')
            pause(0.5)
            set(findobj('Tag','error'),'Visible','off')   
         end% if
      end% if
      set(findobj('Tag','password'),'UserData',CurrWord, ...
          'String',char('*'*ones(1,length(CurrWord))))      
   end% switch
end% if
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tmp = allowed_characters
tmp = ['ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789' ...
       '<>[]{}()@!?*#=~-+_.,;:�$%&/|\'];   


