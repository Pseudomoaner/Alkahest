function orients = findImageOrients(img,SD)
%FINDIMAGEORIENTS uses the tensor method to find the local orientation of
%the input image.
%
%   INPUTS:
%       -img: greyscale input image. Matrix of floats (doubles).
%       -SD: Standard deviation of smoothing kernal. Increase to increase
%       the spatial scale over which the orientations are calculated.
%
%   OUTPUTS:
%       -orients: output image indicating local orientation at each image
%       point in radians
 
[gX,gY] = imgradientxy(img);
Ixx = gX.*gX;
Iyy = gY.*gY;
Ixy = gX.*gY;

tIxx = imgaussfilt(Ixx,SD,'FilterSize',2*ceil(SD*8)+1); %Extra filter size should suffice to ensure smoothness of image gradients (i.e. no zero values)
tIyy = imgaussfilt(Iyy,SD,'FilterSize',2*ceil(SD*8)+1);
tIxy = imgaussfilt(Ixy,SD,'FilterSize',2*ceil(SD*8)+1);

orients = 0.5 * atan2(2*tIxy,tIyy-tIxx);