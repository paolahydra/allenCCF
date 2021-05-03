function fwireframe = plotWireFrame(T_roi, varargin)
% plot points from a single brain on fwireframe
%
% suggested colors over black brain:
% BrainColors = .75*[1.3 1.3 1.3; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 .9; .65 .4 .25; .7 .95 .3; .7 0 0; .6 0 .7; 1 .6 0]; 
   

% sort out inputs
if nargin == 1
    braincolor = .75*[1.3 1.3 1.3];
    black_brain = true;
    fwireframe = [];
elseif nargin == 2
    braincolor = varargin{1};
    black_brain = true;
    fwireframe = [];
elseif nargin == 3
    braincolor = varargin{1};
    black_brain = varargin{2};
    fwireframe = [];
elseif nargin >= 4
    braincolor = varargin{1};
    black_brain = varargin{2};
    fwireframe = varargin{3};
end
bregma = allenCCFbregma(); % bregma position in reference data space
atlas_resolution = 0.010; % mm


% transform coordinates
ap_pixel = bregma(1) - T_roi.AP_location./atlas_resolution; %OK
ml_pixel = bregma(3) + T_roi.ML_location./atlas_resolution; %OK
dv_pixel = bregma(2) + T_roi.DV_location./atlas_resolution; %OK


% plot in wireframe
if isempty(fwireframe)
    fwireframe = plotBrainGrid([], [], fwireframe, black_brain);
    hold on; 
    fwireframe.InvertHardcopy = 'off';
else
    figure(fwireframe)
end

plot3(ap_pixel, ml_pixel, dv_pixel, '.','linewidth',2, 'color',braincolor,'markers',10);

end