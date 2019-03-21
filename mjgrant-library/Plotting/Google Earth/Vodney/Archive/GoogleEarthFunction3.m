function PlotGE(varargin)
% Trajectory is a two-dimensional matrix with latitudes, longitudes, and
% altitudes as three rows
% Orientation is a two-dimensional matrix with roll, pitch, and yaw as
% three rows (optional, set to zero if unused)


if mod(length(varargin),3) ~= 0
    
    fprintf(1,'\ninput error\n\n');
    
    return
end

for N = 1:((length(varargin)) / 3)
    
    % create kml object
    Path = kml(['trajectory ' num2str(N)]);
    Missile = kml(['missile ' num2str(N)]);

    % separate inputs
    Trajectory = cell2mat(varargin(1 + 3 * (N - 1)));
    modPlace = cell2mat(varargin(2 + 3 * (N - 1)));
    Orientation = cell2mat(varargin(3 + 3 * (N - 1)));


    % define trajectory matricies

    % hardcoded test case ///////////////////////////////////
    Lat = linspace(-60,40,1000);
    Lng = linspace(-200,200,1000);
    Alt = N * 500000 - 00 * Lng.^2;
    % ///////////////////////////////////////////////////////

%     Lat = Trajectory(1,:);
%     Lng = Trajectory(2,:);
%     Alt = Trajectory(3,:);

    % plot trajectory
    Path.plot3(Lng, Lat, Alt, 'name', 'trajectory','lineWidth',1, ...
        'lineColor','FF98FB98','altitudeMode','absolute');

    Path.plot(Lng, Lat,'name','ground track','lineWidth',2, ...
        'linecolor','fb9898FF');

    Path.run;

    % define orientation matricies
    if Orientation == 0

        dLng = [];
        dLat = [];
        dAlt = [];

        for I = 1:length(modPlace)

            K = modPlace(I);

            if K < length(Lat)

                dLng = [dLng (Lng(K + 1) - Lng(K)) * (pi / 180) * 6731000];
                dLat = [dLat (Lat(K + 1) - Lat(K)) * (pi / 180) * 6731000];
                dAlt = [dAlt Alt(K + 1) - Alt(K)];

            else

                dLng = [dLng (Lng(K) - Lng(K - 1)) * (pi / 180) * 6731000];
                dLat = [dLat (Lat(K) - Lat(K - 1)) * (pi / 180) * 6731000];
                dAlt = [dAlt Alt(K) - Alt(K - 1)];

            end

            Yaw = -[atand(dLat ./ dLng) 0];
            Roll = linspace(0,0,length(modPlace));
            Pitch = -[atand(dAlt ./ sqrt(dLng.^2 + dLat.^2)) 0];

        end


    else

        % use input 
        Pitch = Orientation(1,:);
        Roll = Orientation(2,:);
        Yaw = Orientation(3,:);

    end


    for I = 1:length(modPlace)

        
        K = modPlace(I);
        modelFolder = Missile.createFolder(sprintf('location %i', I));
        modelName = ['model ' num2str(I)];
        
        % set the scale of the model
        modelSize = 2000000;

        % place model (lng, lat, alt, yaw, roll, pitch...)
        modelFolder.model(Lng(K),Lat(K),Alt(K),Yaw(I),Roll(I),Pitch(I), ...
            'model','missile.dae','name',modelName,'scale',modelSize, ...
            'altitudeMode','absolute');

        % set distance and altitude from camera to the model
        camDist = 5;
        
        % place camera (lng, lat, alt, yaw, roll, pitch...)
        degLng_1 = sind(Yaw(I)) * camDist * (6731000 / (Alt(I) + 6731000));
        degLat_1 = cosd(Yaw(I)) * camDist * (6731000 / (Alt(I) + 6731000));
        degLng_2 = cosd(Yaw(I)) * camDist * (6731000 / (Alt(I) + 6731000));
        degLat_2 = -sind(Yaw(I)) * camDist * (6731000 / (Alt(I) + 6731000));
        degAlt = 5000 * (abs(Alt(K) - 6731000) / 6731)^(1/2);

        delLng = degLng_1 * (pi / 180) * 6731000; 
        delLat = degLat_1 * (pi / 180) * 6731000;
        camPitch = -abs(atand(degAlt / sqrt(delLng^2 + delLat^2))) + 90;

%         % camera placement visualization
%         modelFolder.model((Lng(K) + degLng_1),(Lat(K) + degLat_1),(Alt(K) + ...
%             degAlt),Yaw(I)-90,Roll(I),-camPitch + 90,'model','missile.dae', ...
%             'name','right arrow','scale',1000000,'altitudeMode','absolute');
% 
%         modelFolder.model((Lng(K) - degLng_1),(Lat(K) - degLat_1),(Alt(K) + ...
%             degAlt),Yaw(I)+90,Roll(I),-camPitch + 90,'model','missile.dae', ...
%             'name','left arrow','scale',1000000,'altitudeMode','absolute');
% 
%         modelFolder.model((Lng(K) + degLng_2),(Lat(K) + degLat_2),(Alt(K) + ...
%             degAlt),Yaw(I),Roll(I),-camPitch + 90,'model','missile.dae', ...
%             'name','back arrow','scale',1000000,'altitudeMode','absolute');
% 
%         modelFolder.model((Lng(K) - degLng_2),(Lat(K) - degLat_2),(Alt(K) + ...
%             degAlt),Yaw(I) + 180,Roll(I),-camPitch + 90,'model','missile.dae', ...
%             'name','front arrow','scale',1000000,'altitudeMode','absolute');

        % camera placement
        modelFolder.camera((Lng(K) + degLng_1),(Lat(K) + degLat_1), ...
            (Alt(K) + degAlt),Yaw(I) + 180,camPitch,0,'name','right view');

        modelFolder.camera((Lng(K) - degLng_1),(Lat(K) - degLat_1), ...
            (Alt(K) + degAlt),Yaw(I),camPitch,0,'name','left view');  

        modelFolder.camera((Lng(K) + degLng_2),(Lat(K) + degLat_2), ...
            (Alt(K) + degAlt),Yaw(I) -90,camPitch,0,'name','back view');

        modelFolder.camera((Lng(K) - degLng_2),(Lat(K) - degLat_2), ...
            (Alt(K) + degAlt),Yaw(I) + 90,camPitch,0,'name','front view');     

    end

    % send to Google Earth
    Missile.run;

end
end