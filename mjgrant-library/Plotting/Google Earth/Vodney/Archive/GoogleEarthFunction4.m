function PlotGE(varargin)
%% Inputs
% varargin expects inputs in sets of three matricies:
%   - trajectory matrix: latitudes in row 1, longitudes in row 2, and
%     altitude (absolute) in row 3
%   - model placement: one dimensional array containing the column index
%     of the trajectory matrix, and those columns (lat, lng, alt) are where
%     the models are placed (NOTE: length of the matrix does not matter, 
%     input 0 for model placement if models are not wanted)
%   - orientation matrix: pitch values in row 1, roll in row 2, yaw in
%     row 3 (NOTE: has to match dimensions of trajectory matrix, otherwise
%     input 0 and models will be placed tangent to trajectory by default)
%
% for plotting multiple trajectories, input another set of the three
% matricies (for example: varargin = [trajectory1 modPlace1 Orientation1
% Trajectory2 modPlace2 Orientation2 ... etc])
%
% the function will output an error message if the number of inputs is not
% a multiple of three

if mod(length(varargin),3) ~= 0
    
    fprintf(1,'\ninput error\n\n');
    
    return
end

%% Trajectory Plotting

for N = 1:((length(varargin)) / 3)
    
    % create kml file
    Path = kml(['trajectory ' num2str(N)]);
    Missile = kml(['missile ' num2str(N)]);

    % separate inputs
    Trajectory = cell2mat(varargin(1 + 3 * (N - 1)));
    modPlace = cell2mat(varargin(2 + 3 * (N - 1)));
    Orientation = cell2mat(varargin(3 + 3 * (N - 1)));


    % define trajectory matricies
    
    % hardcoded test case ///////////////////////////////////
%     Lat = linspace(-60,40,1000);
%     Lng = linspace(-200,200,1000);
%     Alt = N * 500000 - 00 * Lng.^2;
    % ///////////////////////////////////////////////////////
    
    Lat = Trajectory(1,:);
    Lng = Trajectory(2,:);
    Alt = Trajectory(3,:);

    % plot trajectory
    Path.plot3(Lng, Lat, Alt, 'name', 'trajectory','lineWidth',1, ...
        'lineColor','FF98FB98','altitudeMode','absolute');

    % plot groundtrack
    Path.plot(Lng, Lat,'name','ground track','lineWidth',2, ...
        'linecolor','fb9898FF');

    % send to Google Earth
    Path.run;
    
%% Orientation Preparation for Model Placement

    % define orientation matricies
    if Orientation == 0

        dLng = [];
        dLat = [];
        dAlt = [];

        % calculate tangential orientations
        for I = 1:length(modPlace)

            K = modPlace(I);

            % find the change in value between indexed trajectory point and
            % previous trajectory point
            if K < length(Lat)

                dLng = [dLng (Lng(K + 1) - Lng(K)) * (pi / 180) * 6731000];
                dLat = [dLat (Lat(K + 1) - Lat(K)) * (pi / 180) * 6731000];
                dAlt = [dAlt Alt(K + 1) - Alt(K)];

            else

                dLng = [dLng (Lng(K) - Lng(K - 1)) * (pi / 180) * 6731000];
                dLat = [dLat (Lat(K) - Lat(K - 1)) * (pi / 180) * 6731000];
                dAlt = [dAlt Alt(K) - Alt(K - 1)];

            end

            % define orientation values from calculated deltas
            Yaw = -[atand(dLat ./ dLng) 0];
            Roll = linspace(0,0,length(modPlace));
            Pitch = -[atand(dAlt ./ sqrt(dLng.^2 + dLat.^2)) 0];

        end


    else

        % use input values for orientation
        Pitch = Orientation(1,:);
        Roll = Orientation(2,:);
        Yaw = Orientation(3,:);

    end


    % iterate through given indicies for model placement
    for I = 1:length(modPlace)

        % create folder for specified location
        K = modPlace(I);
        modelFolder = Missile.createFolder(sprintf('location %i', I));
        modelName = ['model ' num2str(I)];
        
        % set the scale of the model
        % (NOTE: manually set this value if different result is desired)
        modelSize = 2000000;

        % place model (lng, lat, alt, yaw, roll, pitch...)
        modelFolder.model(Lng(K),Lat(K),Alt(K),Yaw(I),Roll(I),Pitch(I), ...
            'model','missile.dae','name',modelName,'scale',modelSize, ...
            'altitudeMode','absolute');

        % set distance and altitude from camera to the model
        % (NOTE: manually set this value if different result is desired)
        camDist = 5;
        
        % calculate camera orientations based on missile orientations
        degLng_1 = sind(Yaw(I)) * camDist * (6731000 / (Alt(I) + 6731000));
        degLat_1 = cosd(Yaw(I)) * camDist * (6731000 / (Alt(I) + 6731000));
        degLng_2 = cosd(Yaw(I)) * camDist * (6731000 / (Alt(I) + 6731000));
        degLat_2 = -sind(Yaw(I)) * camDist * (6731000 / (Alt(I) + 6731000));
        degAlt = 5000 * (abs(Alt(K) - 6731000) / 6731)^(1/2);

        % define camera orientation values from calculated deltas
        delLng = degLng_1 * (pi / 180) * 6731000; 
        delLat = degLat_1 * (pi / 180) * 6731000;
        camPitch = -abs(atand(degAlt / sqrt(delLng^2 + delLat^2))) + 90;

        % place camera (lng, lat, alt, yaw, roll, pitch...)
        modelFolder.camera((Lng(K) + degLng_1),(Lat(K) + degLat_1), ...
            (Alt(K) + degAlt),Yaw(I) + 180,camPitch,0,'name','left view');

        modelFolder.camera((Lng(K) - degLng_1),(Lat(K) - degLat_1), ...
            (Alt(K) + degAlt),Yaw(I),camPitch,0,'name','right view');  

        modelFolder.camera((Lng(K) + degLng_2),(Lat(K) + degLat_2), ...
            (Alt(K) + degAlt),Yaw(I) -90,camPitch,0,'name','front view');

        modelFolder.camera((Lng(K) - degLng_2),(Lat(K) - degLat_2), ...
            (Alt(K) + degAlt),Yaw(I) + 90,camPitch,0,'name','back view');     

    end

    % send to Google Earth
    Missile.run;

end
end