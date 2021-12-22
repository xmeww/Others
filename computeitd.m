function delay = computeitd(azimuth,varargin)
% Computes the Interaural Time Difference given the location of the sound source. 
% 
% DETAILS: 
%   Computes the amount of time with which the LEFT channel has to be
%   delayed with respect to the right in order to achieve the required
%   lateralization effect. 
%   The computation is based on Figure 7.3 from the book 'An Introduction 
%   to the Psychology of Hearing' from Brian C. Moore. 
%   Assumptions: 
%       1) distant sound sources
%       2) the head is regarded as 'two holes separated by a spherical 
%          obstacle'. The 'spherical obstacle' has a radius given in the
%          headRadius argument. 
% 
% USAGE: 
%   delay = computeitd(azimuth)
%   delay = computeitd(azimuth,Name,Value)
% 
% INPUT: 
%   azimuth: vector of the locations of the sound sources with respect to 
%       the observer in the azimuthal plane in degrees (azimuth). 
%       Straight ahead corresponds to 0 degrees, sources from the left
%       have negative azimuth values, sources from the right have
%       positive azimuth values. 
%       The value has to be between -90 and 90 degrees (sounds sources
%       ahead of, not behind of the observer). 
%   'Name' - Value arguments: 
%       headRadius: the radius of the head (assumed to be spherical). The
%           default is 0.09 metres (9 cm) according to the book. 
%       soundSpeed: the speed of the sound in metres per second. The 
%           default is 343.2 metres per second corresponding to the speed 
%           of sound in dry air at 20 ?C. 
%   
% OUTPUT: 
%   delay: vector containging the amount of time in seconds with which 
%       the LEFT channel has to be delayed in order to achieve the required
%       lateralization effect corresponding to the specified azimuth values. 
%       For sound sources to the left of the observer the delay is a 
%       negative number (as well as the location of the sound source in 
%       the azimuthal plane) and corresponds to advancing the left
%       channel with the given amount (this is however usually done by
%       delaying the right channel). 

 
%% Checking Matlab version
if verLessThan('matlab', '8.1.0.604')
    error('In order to run the function, 2013a or newer version of Matlab is needed. ');
end

%% Parsing input
p = inputParser;

checkAzimuth = @(x) all(x >= -90) && all(x <= 90);
checkPosScalar = @(x) isscalar(x) && x > 0;

addRequired(p,'azimuth',checkAzimuth);
addParamValue(p,'headRadius',0.09,checkPosScalar);
addParamValue(p,'soundSpeed',343.2,checkPosScalar);


parse(p,azimuth,varargin{:});

azimuth = p.Results.azimuth;
headRadius = p.Results.headRadius;
soundSpeed = p.Results.soundSpeed;

%% 
% d = r(sin(theta) + theat)
delay = headRadius * (sind(azimuth) + (azimuth * pi / 180)) / soundSpeed;

end

