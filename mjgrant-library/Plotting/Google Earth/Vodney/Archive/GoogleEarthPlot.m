function [] = GoogleEarthPlot(results_file, cont_index, skip, model, orient)
%% Inputs
%  
%  results_file     : path to the results MAT file
%  cont_index       : index of the continuation case to plot
%  skip             : number of trajectories to skip
%  model (optional) : Input passed through to the Google earth function
%  orient(optional) : Input passed through to the Google earth function
%
    load(results_file);
    trajdata_in = out.setCONT(cont_index).CONT;
    
    ctr = 0;
    arguments = cell(3*length(trajdata_in)/skip,1);
    for i=1:skip:length(trajdata_in)
        lon = trajdata_in(i).sol.y(2,:)*180/pi;
        lat = trajdata_in(i).sol.y(3,:)*180/pi;
        alt = trajdata_in(i).sol.y(1,:) - trajdata_in(i).in.const.re{1};
        
        traj = [lat; lon; alt];
        if nargin < 5
            orient = zeros(size(traj));
        end
        if nargin < 4
            if i == 1
                model = 1;
            else
                model = [];
            end
        end
        arguments{3*ctr+1} = traj;
        arguments{3*ctr+2} = model;
        arguments{3*ctr+3} = orient;
        ctr = ctr + 1;
    end
    GoogleEarthFunction4(arguments{:});
end

