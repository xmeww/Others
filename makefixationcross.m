function [fixImage,fixRect] = makefixationcross(fixWidth,lineWidth,backgroundColor,lineColor)
% Makes a fixation cross. 
% 
% USAGE: 
%   [image,rect] = makefixationcross(fixWidth,lineWidth,ppd,backgroundColor,lineColor)
% DETAILS: 
%   Returns a fixation cross image and the rectangle enclosing the 
%   image (in pixels). 
% INPUT: 
%   fixWidth: width of the fixation cross in pixels. 
%   lineWidth: width of the line of the cross in pixels. 
%   backgroundColor: color of the background (grayscale or [r g b] triplet)
%   lineColor: color of the line (grayscale or [r g b] triplet)
% OUTPUT:
%   fixImage: the fixation cross image. 
%   fixRect: coordinates of the rectangle enclosing the image (in pixels), top
%         left corner is (0,0). 

%% Parsing input
p = inputParser;

checkpositivescalar = @(x) isscalar(x) && x > 0;
checkcolor = @(x) (all(size(x) == [1 1]) || all(size(x) == [1 3])) && all(x >= 0 & x <= 255);

addRequired(p,'fixWidth',checkpositivescalar);
addRequired(p,'lineWidth',checkpositivescalar);
addRequired(p,'backgroundColor',checkcolor);
addRequired(p,'lineColor',checkcolor);

parse(p,fixWidth,lineWidth,backgroundColor,lineColor);

fixWidth = p.Results.fixWidth;
lineWidth = p.Results.lineWidth;
backgroundColor = p.Results.backgroundColor;
lineColor = p.Results.lineColor;

%% Main body of the function

% If the width of the line is greater than the fixationcorss's width swap
% them. 
if fixWidth < lineWidth
   temp = fixWidth;
   fixWidth = lineWidth;
   lineWidth = temp; 
end

% The width of the line must be at least one pixel.
if lineWidth < 1
    lineWidth = 1;
end

% If grayscale colors are given, expand them to [r g b] vectors. 
if all(size(backgroundColor) == [1 1]);
    backgroundColor = repmat(backgroundColor,1,3);
end

if all(size(lineColor) == [1 1]);
    lineColor = repmat(lineColor,1,3);
end

% This makes sure that both widths are even or both are odd.
if mod(lineWidth,2) ~= mod(fixWidth,2)
    fixWidth = fixWidth + 1;
end

% Create an image of fixation cross
fixImage = ones(fixWidth,fixWidth,3);
for i = 1:size(backgroundColor,2)
   fixImage(:,:,i) = fixImage(:,:,i) * backgroundColor(i);
end

for i = 1:size(lineColor,2)
    fixImage(((fixWidth-lineWidth)/2)+1:((fixWidth-lineWidth)/2)+lineWidth,:,i) = lineColor(i);
    fixImage(:,((fixWidth-lineWidth)/2)+1:((fixWidth-lineWidth)/2)+lineWidth,i) = lineColor(i);
end

% The coordinates of the enclosing rectangle.
fixRect = [0 0 fixWidth fixWidth];

end

