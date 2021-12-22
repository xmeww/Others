function posDeg = pospix2deg(posPix,viewDist,res)
% Converts positions in pixels to visual angle given the viewing distance and resolution. 
% 
% USAGE: 
%   posDeg = pospix2deg(posPix,viewDist,res) 
%
% INPUT:
%   posPix: position in pixels, numeric . 
%   viewDist: The viewing distance in mm, i.e the distance of the 
%       observer's eyes from the monitor (see DETAILS for assumptions),
%       scalar. 
%   res: screen resolution in pixels/mm, scalar. 
% 
% OUTPUT: 
%   posDeg: the position in degrees. 
% 
% DETAILS:
%   It is assumed that the observer is point-like. 

%% Checking Matlab version
if verLessThan('matlab', '8.1.0.604')
    error('In order to run the function, 2013a or newer version of Matlab is needed. ');
end

%% Parsing input
p = inputParser;

addRequired(p,'posPix',@isnumeric);
addRequired(p,'viewDist',@isscalar);
addRequired(p,'res',@isscalar);

parse(p,posPix,viewDist,res);

posPix = p.Results.posPix;
viewDist = p.Results.viewDist;
res = p.Results.res;

%%
posDeg = atand(posPix/(viewDist*res));

end

