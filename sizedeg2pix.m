function sizePix = sizedeg2pix(sizeDeg,posDeg,viewDist,res)
% Converts sizes in visual angle to pixels given the position of the object, the viewing distance and the resolution. 
% 
% DETAILS:
%   It is assumed that the observer is point-like (cyclops).  
% 
% USAGE: 
%   sizePix = posdeg2pix(sizeDeg,posDeg,viewDist,res) 
%
% INPUT:
%   sizeDeg: the size in degrees visual angle (either azimuth or 
%       elevation).
%   posDeg: the position in degrees visual angle (either azimuth or 
%       elevation). 
%   viewDist: The viewing distance in any unit, i.e the distance of the 
%       observer'seyes from the monitor. (see DETAILS for assumptions)
%   res: screen resolution in pixels/ the unit of viewDist. 
% 
% OUTPUT: 
%   posPix: the position in pixels. 

%% Checking Matlab version
if verLessThan('matlab', '8.1.0.604')
    error('In order to run the function, 2013a or newer version of Matlab is needed. ');
end

%% Parsing input
p = inputParser;

checkPosDeg = @(x) all(size(x) == size(sizeDeg));
addRequired(p,'sizeDeg',@isvector);
addRequired(p,'posDeg',checkPosDeg);
addRequired(p,'viewDist',@isscalar);
addRequired(p,'res',@isscalar);

parse(p,sizeDeg,posDeg,viewDist,res);

sizeDeg = p.Results.sizeDeg;
posDeg = p.Results.posDeg;
viewDist = p.Results.viewDist;
res = p.Results.res;

%%

nearPart = (viewDist.*sind(sizeDeg/2))./(cosd(posDeg).*sind(repmat(90,size(posDeg))+posDeg-(sizeDeg/2)));
farPart = (viewDist.*sind(sizeDeg/2))./(cosd(posDeg).*sind(repmat(90,size(posDeg))-posDeg-(sizeDeg/2)));

sizePix = (nearPart+farPart)*res;

end

