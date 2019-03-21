function [] = REV4googleEarth(in)
    
    % iterate through number of in.googleEarth(i)
    for i = 1:size(in.googleEarth,2)
        
       % load result file related to i^th data set in in.googleEarth
       load(in.googleEarth(i).resultsFile);
       trajdata_in = out.setCONT(in.googleEarth(i).contIndex).CONT;
        
        % create KML file
        if isfield(in.googleEarth(i), 'kmlFile') && ~isempty(in.googleEarth(i).kmlFile)
            plotKML = kml([in.googleEarth(i).kmlFile ' plot ' num2str(i)]);
        else
            plotKML = kml(['plot ' num2str(i)]);
        end
        
        if ~isfield(in.googleEarth(i), 'offTrajectory') || isempty(in.googleEarth(i).offTrajectory)
            % check for skip value
            if ~isfield(in.googleEarth(i), 'skip') || isempty(in.googleEarth(i).skip)
                in.googleEarth(i).skip = 1;
            end
            
            % collect trajectory data
            trajNum = 0;
            for j = 1:(in.googleEarth(i).skip):length(trajdata_in)
                % determines index j for in.googleEarth(i).trajectory(j)
                trajNum = trajNum + 1;

                % assign trajectory data
                in.googleEarth(i).trajectory(trajNum).lon = trajdata_in(j).sol.y(2,:) * 180/pi;
                in.googleEarth(i).trajectory(trajNum).lat = trajdata_in(j).sol.y(3,:) * 180/pi;
                in.googleEarth(i).trajectory(trajNum).alt = trajdata_in(j).sol.y(1,:) - trajdata_in(j).in.const.re{1};

            end

            if isfield(in.googleEarth(i), 'trajectory') && ~isempty(in.googleEarth(i).trajectory)
                % create KML trajectories folders
                folderTrajectories = plotKML.createFolder('trajectories');
            
                % plot trajectories
                for k = 1:length(in.googleEarth(i).trajectory)
                    % initialize structure if nonexistent
                    if ~isfield(in.googleEarth(i).trajectory(k), 'altTrack')
                        in.googleEarth(i).trajectory(k).altTrack = [];
                    end

                    % plot altTrack
                    if in.googleEarth(i).trajectory(k).togAltTr
                        % determine opacity of altTrack
                        if isfield(in.googleEarth(i).trajectory(k).altTrack, 'opacity') && ~isempty(in.googleEarth(i).trajectory(k).altTrack.opacity)
                            altOpac = in.googleEarth(i).trajectory(k).altTrack.opacity;
                        else
                            altOpac = 'FF';
                        end
                        
                        if isfield(in.googleEarth(i).trajectory(k).altTrack, 'color') && isnumeric(in.googleEarth(i).trajectory(k).altTrack.color) && ~isempty(in.googleEarth(i).trajectory(k).altTrack.color)
                            % format RGB color data
                            altColor = fliplr(dec2base(in.googleEarth(i).trajectory(k).altTrack.color, 16)');

                            % combine color and opacity hex strings
                            altColor = [altOpac altColor(:)'];
                        else
                            % determine color of altTrack
                            if isfield(in.googleEarth(i).trajectory(k).altTrack, 'color') && ~isempty(in.googleEarth(i).trajectory(k).altTrack.color)
                                altColor = in.googleEarth(i).trajectory(k).altTrack.color;
                            else
                                altColor = 'g';
                            end

                            % determine saturation of altTrack
                            if isfield(in.googleEarth(i).trajectory(k).altTrack, 'saturation') && ~isempty(in.googleEarth(i).trajectory(k).altTrack.saturation)
                                altSat = in.googleEarth(i).trajectory(k).altTrack.saturation;
                            else
                                altSat = 'FF';
                            end

                            % convert color to hex values
                            altColor = color2Hex(altColor, altSat, altOpac);
                        end

                        % determine altWidth of altTrack
                        if isfield(in.googleEarth(i).trajectory(k).altTrack, 'altWidth') && ~isempty(in.googleEarth(i).trajectory(k).altTrack.altWidth)
                            altWidth = in.googleEarth(i).trajectory(k).altTrack.altWidth;
                        else
                            altWidth = 2;
                        end

                        % plot altTrack KML
                        folderTrajectories.plot3(in.googleEarth(i).trajectory(k).lon,...
                            in.googleEarth(i).trajectory(k).lat,...
                            in.googleEarth(i).trajectory(k).alt, 'name',...
                            ['altTrack ' num2str(k)], 'lineWidth', altWidth,...
                            'lineColor', altColor, 'altitudeMode','absolute'); 
                    end

                    % initialize structure if nonexistent
                    if ~isfield(in.googleEarth(i).trajectory(k), 'groundTrack')
                        in.googleEarth(i).trajectory(k).groundTrack = [];
                    end

                    % plot groundTrack
                    if in.googleEarth(i).trajectory(k).togGrndTr
                        % determine opacity of groundTrack
                        if isfield(in.googleEarth(i).trajectory(k).groundTrack, 'opacity') && ~isempty(in.googleEarth(i).trajectory(k).groundTrack.opacity)
                            groundOpac = in.googleEarth(i).trajectory(k).groundTrack.opacity;
                        else
                            groundOpac = 'FF';
                        end
                        
                        if isfield(in.googleEarth(i).trajectory(k).groundTrack, 'color') && isnumeric(in.googleEarth(i).trajectory(k).groundTrack.color) && ~isempty(in.googleEarth(i).trajectory(k).groundTrack.color)
                            % format RGB color data
                            groundColor = fliplr(dec2base(in.googleEarth(i).trajectory(k).groundTrack.color, 16)');

                            % combine color and opacity hex strings
                            groundColor = [groundOpac groundColor(:)'];
                        else
                            % determine color of groundTrack
                            if isfield(in.googleEarth(i).trajectory(k).groundTrack, 'color') && ~isempty(in.googleEarth(i).trajectory(k).groundTrack.color)
                                groundColor = in.googleEarth(i).trajectory(k).groundTrack.color;
                            else
                                groundColor = 'r';
                            end

                            % determine saturation of groundTrack
                            if isfield(in.googleEarth(i).trajectory(k).groundTrack, 'saturation') && ~isempty(in.googleEarth(i).trajectory(k).groundTrack.saturation)
                                groundSat = in.googleEarth(i).trajectory(k).groundTrack.saturation;
                            else
                                groundSat = 'FF';
                            end

                            % convert color to hex values
                            groundColor = color2Hex(groundColor, groundSat, groundOpac);
                        end

                        % determine groundWidth of groundTrack
                        if isfield(in.googleEarth(i).trajectory(k).groundTrack, 'groundWidth') && ~isempty(in.googleEarth(i).trajectory(k).groundTrack.groundWidth)
                            groundWidth = in.googleEarth(i).trajectory(k).groundTrack.groundWidth;
                        else
                            groundWidth = 3;
                        end

                        % plot groundTrack KML
                        folderTrajectories.plot(in.googleEarth(i).trajectory(k).lon,...
                            in.googleEarth(i).trajectory(k).lat,...
                            'name', ['groundTrack ' num2str(k)], 'lineWidth',...
                            groundWidth,'lineColor', groundColor); 
                    end

                    % create movie
                    if isfield(in.googleEarth(i).trajectory(k), 'movie') && ~isempty(in.googleEarth(i).trajectory(k).movie)
                        % check for model scale
                        if ~isfield(in.googleEarth(i), 'modelScale') || isempty(in.googleEarth(i).modelScale)
                            in.googleEarth(i).modelScale = 1;
                        end

                        % check for camera mode
                        if ~isfield(in.googleEarth(i).trajectory(k), 'cameraMode') || isempty(in.googleEarth(i).trajectory(k).cameraMode)
                            in.googleEarth(i).trajectory(k).cameraMode = 'behind';
                        end

                        % check for camera distance
                        if ~isfield(in.googleEarth(i).trajectory(k), 'cameraDist') || isempty(in.googleEarth(i).trajectory(k).cameraDist)
                            in.googleEarth(i).trajectory(k).cameraDist = 4;
                        end

                        % write movie to KML file
                        folderTrajectories.modelTour(in.googleEarth(i).trajectory(k).movieT,...
                            in.googleEarth(i).trajectory(k).lon,...
                            in.googleEarth(i).trajectory(k).lat,...
                            in.googleEarth(i).trajectory(k).alt,...
                            in.googleEarth(i).trajectory(k).yaw,...
                            in.googleEarth(i).trajectory(k).pitch,...
                            in.googleEarth(i).trajectory(k).roll,...
                            'model', 'missile.dae','scale',...
                            in.googleEarth(i).modelScale, 'cameraMode',...
                            in.googleEarth(i).trajectory(k).cameraMode,...
                            'cameraDistance', in.googleEarth(i).trajectory(k).cameraDist)
                    end
                end
            end
        end
        
        if ~isfield(in.googleEarth(i), 'offModels') || isempty(in.googleEarth(i).offModels)
            if isfield(in.googleEarth(i), 'modelTraj') && ~isempty(in.googleEarth(i).modelTraj)
                % create KML models folder
                folderModels = plotKML.createFolder('models');

                % place models
                for l = 1:length(in.googleEarth(i).modelTraj)
                    % determine trajectory to place models on
                    trajIndex = in.googleEarth(i).modelTraj(l).trajIndex;

                    % check for model scale
                    if ~isfield(in.googleEarth(i), 'modelScale')
                        in.googleEarth(i).modelScale = 2000000;
                    end

                    % group models in folders by trajectory number
                    folderModTraj_L = folderModels.createFolder(sprintf('trajectory %i', trajIndex));

                    for m = 1:length(in.googleEarth(i).modelTraj(l).index)
                        modelName = ['model ' num2str(m)];

                        % determine point on trajectory to place model
                        modIndex = in.googleEarth(i).modelTraj(l).index(m);

                        % check for specified orientation values
                        if ~isfield(in.googleEarth(i).modelTraj(l), 'yaw')
                            modLon = in.googleEarth(i).trajectory(trajIndex).lon;
                            modLat = in.googleEarth(i).trajectory(trajIndex).lat;
                            modAlt = in.googleEarth(i).trajectory(trajIndex).alt;

                            % approximate tangential orientations
                            if modIndex < length(modLon)
                                dLon = (modLon(modIndex + 1) - modLon(modIndex)) * (pi / 180) * 6731000;
                                dLat = (modLat(modIndex + 1) - modLat(modIndex)) * (pi / 180) * 6731000;
                                dAlt = modAlt(modIndex + 1) - modAlt(modIndex);
                            else
                                dLon = (modLon(modIndex) - modLon(modIndex - 1)) * (pi / 180) * 6731000;
                                dLat = (modLat(modIndex) - modLat(modIndex - 1)) * (pi / 180) * 6731000;
                                dAlt = modAlt(modIndex) - modAlt(modIndex - 1);
                            end

                            in.googleEarth(i).model(l).yaw(m) = -atand(dLat / dLon);
                            in.googleEarth(i).model(l).roll(m) = 0;
                            in.googleEarth(i).model(l).pitch(m) = -atand(dAlt / sqrt(dLon^2 + dLat^2));
                        end     

                        % place model KML
                        folderModTraj_L.model(in.googleEarth(i).trajectory(trajIndex).lon(modIndex),...
                            in.googleEarth(i).trajectory(trajIndex).lat(modIndex),...
                            in.googleEarth(i).trajectory(trajIndex).alt(modIndex),...
                            in.googleEarth(i).model(l).yaw(m),...
                            in.googleEarth(i).model(l).roll(m),...
                            in.googleEarth(i).model(l).pitch(m),...
                            'model', 'missile.dae', 'name', modelName, 'scale',...
                            in.googleEarth(i).modelScale, 'altitudeMode', 'absolute');
                    end          
                end
            end
        end
        
        if ~isfield(in.googleEarth(i), 'offBorders') || isempty(in.googleEarth(i).offBorders)
            if isfield(in.googleEarth(i), 'bordCountry') && ~isempty(in.googleEarth(i).bordCountry)
                % create country border KML folder
                folderCountries = plotKML.createFolder('country constraints');
            
                % create country border polygons
                for n = 1:length(in.googleEarth(i).bordCountry)
                    % create individual country border KML folder
                    folderCountry = folderCountries.createFolder(sprintf('country %i', n));

                    % determine altitude of polygons
                    if isfield(in.googleEarth(i).bordCountry(n), 'alt') && ~isempty(in.googleEarth(i).bordCountry(n).alt)
                        altConst = in.googleEarth(i).bordCountry(n).alt;
                    else
                            altConst = 500;
                    end

                    if isfield(in.googleEarth(i).bordCountry(n), 'togWall') && in.googleEarth(i).bordCountry(n).togWall
                        % determine opacity of wall
                        if isfield(in.googleEarth(i).bordCountry(n), 'wallOpac') && ~isempty(in.googleEarth(i).bordCountry(n).wallOpac)
                            wallOpac = in.googleEarth(i).bordCountry(n).wallOpac;
                        else
                            wallOpac = 'AA';
                        end
                        
                        if isfield(in.googleEarth(i).bordCountry(n), 'wallColor') && isnumeric(in.googleEarth(i).bordCountry(n).wallColor) && ~isempty(in.googleEarth(i).bordCountry(n).wallColor)
                            % format RGB color data
                            wallColor = fliplr(dec2base(in.googleEarth(i).bordCountry(n).wallColor, 16)');

                            % combine color and opacity hex strings
                            wallColor = [wallOpac wallColor(:)'];
                        else
                            % determine wall color
                            if isfield(in.googleEarth(i).bordCountry(n), 'wallColor') && ~isempty(in.googleEarth(i).bordCountry(n).wallColor)
                                wallColor = in.googleEarth(i).bordCountry(n).wallColor;
                            else
                                wallColor = 'k';
                            end

                            % determine saturation of wall
                            if isfield(in.googleEarth(i).bordCountry(n), 'wallSat') && ~isempty(in.googleEarth(i).bordCountry(n).wallSat)
                                wallSat = in.googleEarth(i).bordCountry(n).wallSat;
                            else
                                wallSat = 'FF';
                            end

                            % convert color to hex values
                            wallColor = color2Hex(wallColor, wallSat, wallOpac);
                        end

                        % plot polygon wall around country border
                        bordLon = in.googleEarth(i).bordCountry(n).bordLon;
                        bordLat = in.googleEarth(i).bordCountry(n).bordLat;
                        for o = 1:numel(bordLon)-1;
                            folderCountry.poly3([bordLon(o) bordLon(o+1) bordLon(o+1) bordLon(o) bordLon(o)],...
                                [bordLat(o) bordLat(o+1) bordLat(o+1) bordLat(o) bordLat(o)],...
                                [0 0 1 1 0] * altConst, 'polyColor', wallColor,...
                                'extrude', false, 'lineWidth', 0);
                        end
                    end

                    if isfield(in.googleEarth(i).bordCountry(n), 'togCap') && in.googleEarth(i).bordCountry(n).togCap
                        % determine opacity of cap
                        if isfield(in.googleEarth(i).bordCountry(n), 'capOpac') && ~isempty(in.googleEarth(i).bordCountry(n).capOpac)
                            capOpac = in.googleEarth(i).bordCountry(n).capOpac;
                        else
                            capOpac = 'CC';
                        end
                        
                        if isfield(in.googleEarth(i).bordCountry(n), 'capColor') && isnumeric(in.googleEarth(i).bordCountry(n).capColor) && ~isempty(in.googleEarth(i).bordCountry(n).capColor)
                            % format RGB color data
                            capColor = fliplr(dec2base(in.googleEarth(i).bordCountry(n).capColor, 16)');

                            % combine color and opacity hex strings
                            capColor = [capOpac capColor(:)'];
                        else
                            % determine cap color
                            if isfield(in.googleEarth(i).bordCountry(n), 'capColor') && ~isempty(in.googleEarth(i).bordCountry(n).capColor)
                                capColor = in.googleEarth(i).bordCountry(n).capColor;
                            else
                                capColor = 'y';
                            end

                            % determine saturation of cap
                            if isfield(in.googleEarth(i).bordCountry(n), 'capSat') && ~isempty(in.googleEarth(i).bordCountry(n).capSat)
                                capSat = in.googleEarth(i).bordCountry(n).capSat;
                            else
                                capSat = 'FF';
                            end

                            % convert color to hex values
                            capColor = color2Hex(capColor, capSat, capOpac);
                        end
                        
                        % plot floating country-shaped polygon
                        bordLon = in.googleEarth(i).bordCountry(n).bordLon;
                        bordLat = in.googleEarth(i).bordCountry(n).bordLat;
                        folderCountry.poly3(bordLon, bordLat,...
                            linspace(1,1,length(bordLon)) * altConst,...
                            'polyColor', capColor, 'extrude', false, 'lineWidth', 0);
                    end
                end
            end
        end
        
        if ~isfield(in.googleEarth(i), 'offImpact') || isempty(in.googleEarth(i).offImpact)
            if isfield(in.googleEarth(i), 'impactKE') && ~isempty(in.googleEarth(i).impactKE)
                % create impact energy KML folder
                folderImpactKE = plotKML.createFolder('impact points');
            
                % create energy gradient polygons
                for p = 1:length(in.googleEarth(i).impactKE)
                    % create folder for each gradient
                    folderGradient = folderImpactKE.createFolder(sprintf('impact %i', p));

                    % set latitude and longitude variables for impact energy data
                    coordLon = in.googleEarth(i).impactKE(p).coordLon;
                    coordLat = in.googleEarth(i).impactKE(p).coordLat;

                    % determine lower and upper bounds of impact data
                    minKE = min(in.googleEarth(i).impactKE(p).energy);
                    maxKE = max(in.googleEarth(i).impactKE(p).energy);
                    
                    % check for RGB values
                    colorKE = in.googleEarth(i).impactKE(p).colors;
                    if isfield(in.googleEarth(i).impactKE(p), 'colors') && isnumeric(in.googleEarth(i).impactKE(p).colors) && ~isempty(in.googleEarth(i).impactKE(p).colors)
                        colorNumbers = fliplr(in.googleEarth(i).impactKE(p).colors);
                        in.googleEarth(i).impactKE(p).saturation = '01';
                    elseif isfield(in.googleEarth(i).impactKE(p), 'colors') && ~isempty(in.googleEarth(i).impactKE(p).colors)
                        % translate colors to numbers
                        for q = 1:size(colorKE, 1)
                            if (colorKE(q) == 'y')
                                colorNumbers(q, :) = [0 1 1];
                            elseif (colorKE(q) == 'm')
                                colorNumbers(q, :) = [1 0 1];
                            elseif (colorKE(q) == 'c')
                                colorNumbers(q, :) = [1 1 0];
                            elseif (colorKE(q) == 'r')
                                colorNumbers(q, :) = [0 0 1];
                            elseif (colorKE(q) == 'g')
                                colorNumbers(q, :) = [0 1 0];
                            elseif (colorKE(q) == 'b')
                                colorNumbers(q, :) = [1 0 0];
                            elseif (colorKE(q) == 'w')
                                colorNumbers(q, :) = [1 1 1];
                            else
                                colorNumbers(q, :) = [0 0 0];
                            end 
                        end
                    else
                        colorNumbers = [0 0 1; 0 1 1; 0 1 0];
                    end

                    % determine and plot color values
                    sampleDelta = in.googleEarth(i).impactKE(p).sampleDelta;
                    for r = 1:length(in.googleEarth(i).impactKE(p).energy);
                        valEnergy = in.googleEarth(i).impactKE(p).energy(r);

                        % choose algorithm based on number of colors
                        if size(colorKE, 1) == 1
                            % single color gradient
                            colorPercent = (valEnergy - minKE) / (maxKE - minKE);
                            colorOutput = dec2base(round(colorNumbers(1,:) * colorPercent *...
                                base2dec(in.googleEarth(i).impactKE(p).saturation, 16)), 16, 2)';
                        elseif size(colorKE, 1) == 2
                            % two color gradient
                            colorPercent = (valEnergy - minKE) / (maxKE - minKE);
                            colorOutput = dec2base(round((colorNumbers(1,:) * colorPercent + colorNumbers(2,:) * (1 - colorPercent)) *...
                                base2dec(in.googleEarth(i).impactKE(p).saturation, 16)), 16, 2)';
                        else
                            % three color gradient
                            transVal = minKE + 0.5 * (maxKE - minKE);
                            if (valEnergy < transVal)
                                colorPercent = ((valEnergy - minKE) / (transVal - minKE));
                                colorOutput = dec2base(round((colorNumbers(2,:) * colorPercent + colorNumbers(3,:) * (1 - colorPercent)) *...
                                base2dec(in.googleEarth(i).impactKE(p).saturation, 16)), 16, 2)';
                            else
                                colorPercent = ((valEnergy - transVal) / (maxKE - transVal));
                                colorOutput = dec2base(round((colorNumbers(1,:) * colorPercent + colorNumbers(2,:) * (1 - colorPercent)) *...
                                base2dec(in.googleEarth(i).impactKE(p).saturation, 16)), 16, 2)';
                            end
                        end
                        
                        % finish formatting hex string
                        colorOutput = [in.googleEarth(i).impactKE(p).opacity colorOutput(:)'];

                        % plot polygon individual polygons
                        folderGradient.poly3([(coordLon(r) - sampleDelta) (coordLon(r) + sampleDelta) (coordLon(r) + sampleDelta) (coordLon(r) - sampleDelta) (coordLon(r) - sampleDelta)],...
                            [(coordLat(r) - sampleDelta) (coordLat(r) - sampleDelta) (coordLat(r) + sampleDelta) (coordLat(r) + sampleDelta) (coordLat(r) - sampleDelta)],...
                            [1 1 1 1 1] * in.googleEarth(i).impactKE(p).alt,...
                            'polyColor', colorOutput, 'extrude',false,'lineWidth', 0);
                    end
                end
            end
        end
        
        % send KML file to Google Earth
        plotKML.run;
    end
end