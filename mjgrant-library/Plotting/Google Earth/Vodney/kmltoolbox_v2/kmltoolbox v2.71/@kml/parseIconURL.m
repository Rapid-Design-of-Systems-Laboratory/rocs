function  iconURL = parseIconURL(iconAlias)
%KML.PARSEICONURL(iconAlias) is an internal function to convert an icon alias to
%  the convention used by Google Earth
%  
%  A list of all icons can be found <a href="http://jyotirmaya.blogspot.com/2008/03/google-map-files-kml-icon.html">here</a>.
%
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $
    
    
    persistent allURLs alias
    
    if isempty(iconAlias)
        iconURL = 'http://maps.google.com/mapfiles/kml/pal2/icon18.png';
        return
    end
    
    if ~iscell(iconAlias)
        iconAlias = {iconAlias};
        cellOutput = false;
    else
        cellOutput = true;
    end
    
        
    if isempty(allURLs)

        %List of icons taken from http://jyotirmaya.blogspot.com/2008/03/google-map-files-kml-icon.html
        pal = @(p,n) cellfun(@(x)sprintf(['p' num2str(p) 'i%i'],x),num2cell(1:n),'UniformOutput',false);
        padleIcons = {'wht-stars', 'wht-stars-lv', 'wht-blank', 'wht-blank-lv', 'wht-circle', 'wht-circle-lv', 'wht-circle_maps', 'wht-blank_maps', ...
                      'ylw-stars', 'ylw-stars-lv', 'ylw-blank', 'ylw-blank-lv', 'ylw-circle', 'ylw-circle-lv', 'ylw-circle_maps', 'ylw-blank_maps', ...
                      'red-stars', 'red-stars-lv', 'red-blank', 'red-blank-lv', 'red-circle', 'red-circle-lv', 'red-circle_maps', 'red-blank_maps', ...
                      'grn-stars', 'grn-stars-lv', 'grn-blank', 'grn-blank-lv', 'grn-circle', 'grn-circle-lv', 'grn-circle_maps', 'grn-blank_maps', ...
                      'blu-stars', 'blu-stars-lv', 'blu-blank', 'blu-blank-lv', 'blu-circle', 'blu-circle-lv', 'blu-circle_maps', 'blu-blank_maps', ...
                      'pink-stars', 'pink-stars-lv', 'pink-blank', 'pink-blank-lv', 'pink-circle', 'pink-circle-lv', 'pink-circle_maps', 'pink-blank_maps', ...
                      'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', ...
                      'A-lv', 'B-lv', 'C-lv', 'D-lv', 'E-lv', 'F-lv', 'G-lv', 'H-lv', 'I-lv', 'J-lv', 'K-lv', 'L-lv', 'M-lv', 'N-lv', 'O-lv', 'P-lv', 'Q-lv', 'R-lv', 'S-lv', 'T-lv', 'U-lv', 'V-lv', 'W-lv', 'X-lv', 'Y-lv', 'Z-lv', ...
                      'A_maps', 'B_maps', 'C_maps', 'D_maps', 'E_maps', 'F_maps', 'G_maps', 'H_maps', 'I_maps', 'J_maps', 'K_maps', 'L_maps', 'M_maps', 'N_maps', 'O_maps', 'P_maps', 'Q_maps', 'R_maps', 'S_maps', 'T_maps', 'U_maps', 'V_maps', 'W_maps', 'X_maps', 'Y_maps', 'Z_maps', ...
                      '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', ...
                      '1-lv', '2-lv', '3-lv', '4-lv', '5-lv', '6-lv', '7-lv', '8-lv', '9-lv', '10-lv', ...
                      '1_maps', '2_maps', '3_maps', '4_maps', '5_maps', '6_maps', '7_maps', '8_maps', '9_maps', '10_maps'};
        alias = cat(2,pal(2,63),pal(3,63),pal(4,63),pal(5,63),padleIcons);

        allURLs = { ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon1.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon2.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon3.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon4.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon5.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon6.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon7.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon8.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon9.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon10.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon11.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon12.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon13.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon14.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon15.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon16.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon17.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon18.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon19.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon20.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon21.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon22.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon23.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon24.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon25.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon26.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon27.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon28.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon29.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon30.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon31.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon32.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon33.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon34.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon35.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon36.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon37.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon38.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon39.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon40.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon41.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon42.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon43.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon44.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon45.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon46.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon47.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon48.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon49.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon50.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon51.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon52.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon53.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon54.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon55.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon56.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon57.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon58.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon59.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon60.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon61.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon62.png', ...
                    'http://maps.google.com/mapfiles/kml/pal2/icon63.png', ...
                    ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon1.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon2.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon3.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon4.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon5.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon6.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon7.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon8.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon9.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon10.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon11.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon12.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon13.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon14.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon15.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon16.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon17.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon18.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon19.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon20.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon21.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon22.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon23.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon24.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon25.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon26.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon27.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon28.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon29.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon30.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon31.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon32.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon33.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon34.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon35.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon36.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon37.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon38.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon39.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon40.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon41.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon42.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon43.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon44.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon45.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon46.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon47.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon48.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon49.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon50.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon51.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon52.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon53.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon54.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon55.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon56.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon57.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon58.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon59.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon60.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon61.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon62.png', ...
                    'http://maps.google.com/mapfiles/kml/pal3/icon63.png', ...
                    ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon1.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon2.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon3.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon4.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon5.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon6.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon7.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon8.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon9.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon10.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon11.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon12.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon13.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon14.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon15.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon16.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon17.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon18.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon19.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon20.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon21.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon22.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon23.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon24.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon25.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon26.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon27.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon28.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon29.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon30.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon31.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon32.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon33.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon34.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon35.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon36.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon37.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon38.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon39.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon40.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon41.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon42.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon43.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon44.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon45.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon46.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon47.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon48.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon49.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon50.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon51.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon52.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon53.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon54.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon55.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon56.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon57.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon58.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon59.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon60.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon61.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon62.png', ...
                    'http://maps.google.com/mapfiles/kml/pal4/icon63.png', ...
                    ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon1.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon2.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon3.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon4.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon5.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon6.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon7.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon8.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon9.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon10.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon11.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon12.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon13.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon14.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon15.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon16.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon17.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon18.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon19.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon20.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon21.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon22.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon23.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon24.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon25.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon26.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon27.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon28.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon29.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon30.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon31.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon32.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon33.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon34.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon35.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon36.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon37.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon38.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon39.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon40.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon41.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon42.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon43.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon44.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon45.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon46.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon47.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon48.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon49.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon50.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon51.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon52.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon53.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon54.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon55.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon56.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon57.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon58.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon59.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon60.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon61.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon62.png', ...
                    'http://maps.google.com/mapfiles/kml/pal5/icon63.png', ...
                    ...
                    'http://maps.google.com/mapfiles/kml/paddle/wht-stars.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/wht-stars-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/wht-blank.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/wht-blank-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/wht-circle.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/wht-circle-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/wht-circle_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/wht-blank_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/ylw-stars.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/ylw-stars-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/ylw-blank.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/ylw-blank-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/ylw-circle.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/ylw-circle-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/ylw-circle_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/ylw-blank_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/red-stars.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/red-stars-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/red-blank.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/red-blank-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/red-circle.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/red-circle-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/red-circle_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/red-blank_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/grn-stars.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/grn-stars-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/grn-blank.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/grn-blank-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/grn-circle.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/grn-circle-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/grn-circle_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/grn-blank_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/blu-stars.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/blu-stars-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/blu-blank.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/blu-blank-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/blu-circle.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/blu-circle-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/blu-circle_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/blu-blank_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/pink-stars.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/pink-stars-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/pink-blank.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/pink-blank-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/pink-circle.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/pink-circle-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/pink-circle_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/pink-blank_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/A.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/B.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/C.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/D.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/E.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/F.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/G.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/H.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/I.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/J.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/K.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/L.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/M.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/N.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/O.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/P.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/Q.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/R.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/S.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/T.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/U.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/V.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/W.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/X.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/Y.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/Z.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/A-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/B-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/C-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/D-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/E-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/F-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/G-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/H-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/I-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/J-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/K-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/L-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/M-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/N-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/O-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/P-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/Q-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/R-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/S-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/T-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/U-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/V-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/W-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/X-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/Y-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/Z-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/A_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/B_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/C_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/D_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/E_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/F_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/G_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/H_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/I_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/J_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/K_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/L_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/M_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/N_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/O_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/P_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/Q_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/R_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/S_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/T_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/U_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/V_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/W_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/X_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/Y_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/Z_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/1.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/2.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/3.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/4.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/5.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/6.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/7.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/8.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/9.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/10.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/1-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/2-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/3-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/4-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/5-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/6-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/7-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/8-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/9-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/10-lv.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/1_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/2_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/3_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/4_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/5_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/6_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/7_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/8_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/9_maps.png', ...
                    'http://maps.google.com/mapfiles/kml/paddle/10_maps.png'};

    end
    
    [ism,loc] = ismember(iconAlias,alias);
    loc(~ism) = 18; %pal2 icon18
    isURLalready = cellfun(@(x)~isempty(x),strfind(iconAlias,'http://'));
    isLocal      = cellfun(@(x)exist(x,'file'),iconAlias);

    toCopy = ~(isURLalready | isLocal);
    
    iconURL = iconAlias;
    iconURL(toCopy) = allURLs(loc(toCopy));
    
    if ~cellOutput
        iconURL = iconURL{1};
    end
end