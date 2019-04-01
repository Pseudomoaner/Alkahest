function outImg = paintEllipse(inImg,xs,ys,majors,minors,phis,intensities,pxSize)
%PAINTELLIPSE paints the ellipses of the given dimensions and positions 
%into the image inImg. Paints in grayscale, based on the intensities given 
%in the intensities variable.
%
%   INPUTS:
%       -inImg: a matrix of initial image intensities. Should be mostly
%       ones, but can contain other values, indicating e.g. confinement 
%       geometries.
%       -xs: X-coordinates of ellipses to be painted. Nx1 vector
%       -ys: Y-coordinates of ellipses to be painted. Nx1 vector
%       -majors: Major axes of ellipses. Nx1 vector
%       -minors: Minor axes of ellipses. Nx1 vector
%       -phis: Orientations of ellipses, in degress. Nx1 vector
%       -intensities: Value that should be subtracted from inImg over mask
%       indicated by each ellipse. Nx1 vector
%       -pxSize: Scaling factor that relates pixel size in inImg to
%       physical units used in xs, ys, majors and minors.
%
%   OUTPUTS:
%       -outImg: Result of subtracting each ellipse from inImg
%
%   Author: Oliver J. Meacock

xPxs = ceil(xs/pxSize);
yPxs = ceil(ys/pxSize);

%These are cludges to account for those very rare occasions when a rod is
%over the edge of the final pixel - just subtract or add one to bring it within
%range.
xPxs(xPxs > size(inImg,2)) = xPxs(xPxs > size(inImg,2)) - 1;
yPxs(yPxs > size(inImg,1)) = yPxs(yPxs > size(inImg,1)) - 1;
xPxs(xPxs == 0) = 1;
yPxs(yPxs == 0) = 1;

majLenPxs = majors/pxSize;
minLenPxs = minors/pxSize;

halfWindow = round(max(majLenPxs)) + 3;
fullWindow = 2*halfWindow + 1; %Size of entire cutout

%Prepare image by padding it out.
outImg = inImg;
outImg = [outImg(end-halfWindow+1:end,:);outImg;outImg(1:halfWindow,:)];
outImg = [outImg(:,end-halfWindow+1:end),outImg,outImg(:,1:halfWindow)];

for i = 1:size(xs,1)
    xPx = xPxs(i);
    yPx = yPxs(i);
    majLenPx = majLenPxs(i);
    minLenPx = minLenPxs(i);
    phi = phis(i);
    
    %Generate a cut out bit of the coordinate grid for calculating the ellipse
    %over.
    [xGrid,yGrid] = meshgrid(xPx-halfWindow:xPx+halfWindow,yPx-halfWindow:yPx+halfWindow);
    
    %Transform x and y grids into canonical coordinates.
    xCan = (xGrid-xPx)*cosd(-phi) + (yGrid-yPx)*sind(-phi);
    yCan = -(xGrid-xPx)*sind(-phi) + (yGrid-yPx)*cosd(-phi);
    
    geometry = ((xCan.^2)/(majLenPx^2)) + ((yCan.^2)/(minLenPx^2));
    ellipseImg = geometry < 1.0; %Image of ellipse in 'small' coordinates

    %Do initial drawing
    minXbig = xPx + 2;
    maxXbig = xPx + fullWindow - 3;
    minYbig = yPx + 2;
    maxYbig = yPx + fullWindow - 3;
    
    %Don't need to worry about rods that are outside the actual field,
    %as the periodic boundary conditions will move them back in.
    outImg(minYbig:maxYbig,minXbig:maxXbig) = outImg(minYbig:maxYbig,minXbig:maxXbig) - ellipseImg(3:end-2,3:end-2)*intensities(i);
end

%Cut out padded edges and move them around to do wrap around
%x-dimension first...
leftChunk1 = outImg(:,halfWindow+1:fullWindow);
changeInds1 = leftChunk1 < 0.999;
jointChunk1 = outImg(:,end-halfWindow:end);
jointChunk1(changeInds1) = leftChunk1(changeInds1);
outImg(:,halfWindow+1:fullWindow) = jointChunk1;

rightChunk2 = outImg(:,end-fullWindow:end-halfWindow-2);
changeInds2 = rightChunk2 < 0.999;
jointChunk2 = outImg(:,1:halfWindow);
jointChunk2(changeInds2) = rightChunk2(changeInds2);
outImg(:,end-fullWindow:end-halfWindow-2) = jointChunk2;

%Then y
bottomChunk3 = outImg(halfWindow+1:fullWindow,:);
changeInds3 = bottomChunk3 < 0.999;
jointChunk3 = outImg(end-halfWindow:end,:);
jointChunk3(changeInds3) = bottomChunk3(changeInds3);
outImg(halfWindow+1:fullWindow,:) = jointChunk3;

topChunk4 = outImg(end-fullWindow:end-halfWindow-2,:);
changeInds4 = topChunk4 < 0.999;
jointChunk4 = outImg(1:halfWindow,:);
jointChunk4(changeInds4) = topChunk4(changeInds4);
    outImg(end-fullWindow:end-halfWindow-2,:) = jointChunk4;

outImg(:,end-halfWindow:end) = [];
outImg(:,1:halfWindow) = [];
outImg(end-halfWindow:end,:) = [];
outImg(1:halfWindow,:) = [];