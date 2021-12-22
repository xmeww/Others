function [grating gabor] = gaborPatch(imSize, sigma, lamda)
% imSize (i.e. 200) = image size: n X n
% sigma (i.e. 30) = gaussian standard deviation in pixels
% lamda (i.e. 50) = wavelength (number of pixels per cycle)


theta = 45;                             % grating orientation

phase = 1;                              % phase (0 -> 1)
trim = .005;    

X = 1:imSize;                           % X is a vector from 1 to imageSize
X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5

freq = imSize/lamda;                    % compute frequency from wavelength
phaseRad = (phase * 2* pi);             % convert to radians: 0 -> 2*pi

[Xm Ym] = meshgrid(X0, X0);             % 2D matrices

thetaRad = (theta / 360) * 2*pi;        % convert theta (orientation) to radians
Xt = Xm * cos(thetaRad);                % compute proportion of Xm for given orientation
Yt = Ym * sin(thetaRad);                % compute proportion of Ym for given orientation
XYt =  Xt + Yt ;                        % sum X and Y components
XYf = XYt * freq * 2*pi;                % convert to radians and scale by frequency
grating = sin( XYf + phaseRad);         % make 2D sinewave

s = sigma / imSize;                     % gaussian width as fraction of imageSize

gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian

gauss(gauss < trim) = 0;                 % trim around edges (for 8-bit colour displays)
gabor = grating .* gauss;                % use .* dot-product
