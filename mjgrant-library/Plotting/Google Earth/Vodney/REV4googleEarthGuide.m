% Google Earth in.googleEarth structure guide
%
% NOTE: '(opt)' means it is optional to leave variable empty or not create
%       the structure
%
%       '(1/0)' means boolean input of 1 or 0
%
%       '#xN#' matrix means that number of columns has to match any other
%       matrix with the same 'S#'
%
%       available colors use same string input as Matlab default for
%       plotting
%
%       colors can also be input as RGB values of 0-255 for more color options
%       (for impact energy, put all RGB colors in one matrix with each row
%       being a separate color)
%
% [in.]
% in.googleEarth(i)   : each googleEarth(i) has everything needed to plot
%                       multiple trajectories from one results file
%
% [in.googleEarth(i).]
% kmlFile (opt)       : names KML file so function can run multiple times
%                       without overwriting
% resultsFile         : sets results file to reference for i^th
%                       googleEarth(i) iteration
% trajectory(k)       : contains all data on trajectories
% contIndex           : index of the continuation case to plot
% skip (opt)          : number of trajectories to skip
% modelTraj(l) (opt)  : contains model data for each specified trajectory
% modelScale (opt)    : determines visual size of models
% bordCountry(n)      : contains all data for country borders
% impactKE(p)         : contains all data for impact point gradients
%
% // only use the following four variables to isolate/avoid an error //
% offTrajectory (opt) : any value turns off this section of the function
% offModels (opt)     : any value turns off this section of the function
% offBorders (opt)    : any value turns off this section of the function
% offImpact (opt)     : any value turns off this section of the function
%
% [in.googleEarth(i).trajectory(k).]
% lon                 : 1xN1 matrix of longitudes
% lat                 : 1xN1 matrix of latitudes
% alt                 : 1xN1 matrix of altitudes
% togAltTr (1/0)      : value of 1 turns on altTrack
% togGrndTr (1/0)     : value of 1 turns on groundTrack
% altTrack            : contains all altTrack data
% groundTrack         : contains all groundTrack data
% movie (opt)         : any value turns on this section of the function
% movieT (opt)        : (not optional if movie structure exists) this is a
%                       1xS1 matrix related to the time for the animation
%                       (not really sure what the values refer to, except
%                       that the example uses position and orientation as
%                       functions of this matrix)
% roll (opt)          : 1xN1 matrix of roll values along trajectory (not optional if creating movie)
% pitch (opt)         : 1xN1 matrix of pitch values along trajectory (not optional if creating movie)
% yaw (opt)           : 1xN1 matrix of yaw values along trajectory (not optional if creating movie)
% cameraMode (opt)    : either 'above' or 'behind'
% cameraDist (opt)    : distance from camera to center of model during movie
%                        
% [in.googleEarth(i).trajectory(k).altTrack.]
% color (opt)         : sets color of trajectory
% saturation (opt)    : sets color saturation (hex string 00 to FF)
% opacity (opt)       : sets color opacity (hex string 00 to FF)
% altWidth (opt)      : sets width of line
%
% [in.googleEarth(i).trajectory(k).groundTrack.]
% color (opt)         : sets color of trajectory
% saturation (opt)    : sets color saturation (hex string 00 to FF)
% opacity (opt)       : sets color opacity (hex string 00 to FF)
% groundWidth (opt)   : sets width of line
%
% [in.googleEarth(i).modelTraj(l).] (opt)
% trajIndex           : index of trajectory to place models on
% index               : 1xN2 matrix of indicies to place models along
%                       trajectory
% yaw (opt)           : 1xN2 matrix of yaw for requested models
% roll (opt)          : 1xN2 matrix of roll for requested models
% pitch (opt)         : 1xN2 matrix of pitch for requested models
%
% [in.googleEarth(i).bordCountry(n).] (opt)
% bordLat             : 1xN3 matrix of latitude points on border
% bordLon             : 1xN3 matrix of longitude points on border
%                       (border points should be counterclockwise, or else
%                       color and opacity are messed up)
%                       as long both are the same direction)
% alt (opt)           : height of polygons
% togWall (1/0)       : value of 1 turns on wall polygons
% togCap (1/0)        : value of 1 turns on cap polygons
% wallColor (opt)     : sets wall color (hex string 00 to FF)
% wallSat (opt)       : sets wall color saturation (hex string 00 to FF)
% wallOpac (opt)      : sets wall color opacity (hex string 00 to FF)
% capColor (opt)      : sets cap color (hex string 00 to FF)
% capSat (opt)        : sets cap color saturation (hex string 00 to FF)
% capOpac (opt)       : sets cap color opacity (hex string 00 to FF)
%
% [in.googleEarth(i).impactKE(p).] (opt)
% coordLat            : 1xN4 matrix of 
% coordLon            : 1xN4 matrix of longitude points with impact energy data
% energy              : 1xN4 matrix of energy value corresponding to coordinates
% colors              : single column matrix of up to three rows, with the
%                       gradient colors
% opacity             : opacity of energy gradient (hex string 00 to FF)
% saturation          : saturation of energy gradient colors (hex string 00 to FF)
% alt                 : altitude of gradient polygons
% sampleDelta         : distance between data points (IMPORTANT: makes adjacent polygons touch)