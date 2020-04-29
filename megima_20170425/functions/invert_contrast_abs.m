function imgOut = invert_contrast_abs(imgIn)
%imgOut = invert_contrast_abs(imgIn)
%
% Inverts the contrast of a greyscale image.
%

imgOut = abs(imgIn-1);
