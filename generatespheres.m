function spheres = generatespheres(xstd,ystd,num,sphereDiam,fig)
% Generate a random 2D gaussian cloud of spheres with the desired specifications. 
% 
% DETAILS: 
%   The 2D gaussian cloud is centered at (0,0). All coordinates are
%   in azimuth. New spheres are resampled till all spheres meet the
%   following criteria: 
%       1) std of the spheres must be within the 90-110% of the desired
%           std. 
%       2) all spheres must be within the 95% confidence interval of the 2D
%           gaussian defined by xstd and ystd. 
%       3) the spheres must not overlap. 
% 
% USAGE: 
%   spheres = generatespheres(xstd,ystd,num,sphereDiam,fig)
% 
% INPUT: 
%   xstd: horizontal std of the 2D gaussian in azimuth
%   ystd: vertical std of the 2D gaussian 
%   num: number of spheres in the gaussian cloud
%   sphereDiam: diameter of each sphere 
%   fig: whether to plot the generated gaussian cloud (nonzero = yes)
% 
% OUTPUT: 
%   spheres: num by 2 matrix containing the coordinates of the centers of 
%       the spheres. col 1: x coordinates, col 2: y coordinates. 

 
%% Checking Matlab version
if verLessThan('matlab', '8.1.0.604')
    error('In order to run the function, 2013a or newer version of Matlab is needed. ');
end
%% Parsing input.
p = inputParser;

addRequired(p,'xstd',@isscalar);
addRequired(p,'ystd',@isscalar);
addRequired(p,'num',@isscalar);
addRequired(p,'sphereDiam',@isscalar);
addRequired(p,'fig',@isscalar);

parse(p,xstd,ystd,num,sphereDiam,fig);

xstd = p.Results.xstd;
ystd = p.Results.ystd;
num = p.Results.num;
sphereDiam = p.Results.sphereDiam;
fig = p.Results.fig;

%% Main part of the function

% Just a rough criterion to prevent the selection from going into an
% infinite loop.
if sphereDiam >= min(xstd,ystd)/2
   error(['Sphere diameter is too large for the given stds! ' ...
        'High chance of going into an infinite loop with the selection. ']);
end

% Default randomization for 2D Gaussian
mu = [0 0];
sigma = [xstd^2 0; 0 ystd^2];
spheres = mvnrnd(mu,sigma,num);

% Check whether the spheres are meeting the criteria: 
% 1) std of the spheres must be within the 90-110% of the desired std.
% 2) all spheres must be within the 95% confidence interval of the 2D
%    gaussian defined by xstd and ystd. 
% 3) the spheres must not overlap. - not yet implemented. 
[isStdBad,isOut,isOverlap] = checkspheres(spheres,xstd,ystd,sphereDiam);

while isStdBad || any(isOut) || any(isOverlap)
    if isStdBad
        % If the std is bad, obviously we have to generate a completely new
        % set of spheres. 
        spheres = mvnrnd(mu,sigma,num);
        % Recheck the new set of spheres. 
        [isStdBad,isOut,isOverlap] = checkspheres(spheres,xstd,ystd,sphereDiam);
    else
        % If there are spheres outside of the 95% confidence interval of
        % the 2D gaussian, or any of those are overlapping, we just 
        % resample those affected. 
        while any(isOut) || any(isOverlap)
            spheres(isOut | isOverlap,:) = mvnrnd(mu,sigma,sum(isOut | isOverlap));
            % Recheck the new spheres.
            [isStdBad,isOut,isOverlap] = checkspheres(spheres,xstd,ystd,sphereDiam);
        end
    end  
end

% Plot if needed
if fig
    figure;
    N = 50;
    theta = 0:1/N:2*pi;
    x = cos(theta) * 2 * xstd; 
    y = sin(theta) * 2 * ystd;
    plot(x,y,'r');
    hold on
    viscircles(spheres,repmat(sphereDiam/2,num,1),'EdgeColor','b');
    axis equal;
end

end


function [isStdBad,isOut,isOverlap] = checkspheres(spheres,xstd,ystd,sphereDiam)

% Check whether the std of the generated spheres is outside 90-110% of the
% desired std.
isStdBad = std(spheres(:,1)) > xstd*1.1 || std(spheres(:,1)) < xstd*0.9;

% Check and save the ids of spheres which are outside the 95% confidence
% interval of the 2D gaussian. This is based on the equation of an ellipse
% with diameters 2*xstd and 2*ystd minus the radius of the spheres.
isOut = (spheres(:,1)./((2*xstd)-(sphereDiam/2))).^2 + (spheres(:,2)./((2*ystd)-(sphereDiam/2))).^2 > 1;

% Check the Eucledian distance of the spheres pairwise and find those pairs
% which distance is less than the sphere's diameter (basically
% if they overlap). Return just one member of each pair (for
% resampling). 
% find the indexes of one member of overlapping pairs. 
[~,ids] = find(tril(squareform(pdist(spheres) <= sphereDiam)));
% Logical vector of length number of spheres indicating the above mentioned
% spheres. 
isOverlap = zeros(size(spheres,1),1);
isOverlap(unique(ids)) = 1;

end