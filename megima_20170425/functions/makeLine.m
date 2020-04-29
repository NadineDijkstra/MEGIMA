function img = makeLine(iw, ih, DiamPix, inrad, outrad, pixperdeg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin<6, pixperdeg = 1; end
if nargin<5, outrad = inf; end
if nargin<4, inrad = 0; end
if nargin<3, DiamPix = 1; end
inrad = inrad * pixperdeg; 
outrad = outrad * pixperdeg;


iw = roundEven(iw); %These must be even to be able to precisely center the image on the screen (which also has an even number of pixels)
ih = roundEven(ih);
DiamPix = roundEven(DiamPix);

range_x = linspace(iw/-2, iw/2, iw);
range_y = linspace(ih/-2, ih/2, ih);
[x, y] = meshgrid(range_x, range_y);
img = (y > (0-DiamPix/2) & y < (DiamPix/2));

eccen = sqrt(x.^2 + y.^2);
img(or(eccen > outrad, eccen < inrad)) = 0;



end

