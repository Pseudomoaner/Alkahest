function [] = simulateWensinkField(startField,fS,dS,ax)
%SIMULATEWENSINKFIELD performs the main simulation of Alakhest. The system
%is visualised and saved at each sampling timepoint.
%
%   INPUTS:
%       -startField: a WensinkField object
%       -fS: field settings associated with startField
%       -dS: display settings, indicating where output images should be
%       saved
%       -ax: handle to axis where output images should be displayed
%       -compiled: whether the energy gradient calculation function has
%       been compiled on this system
%
%   Author: Oliver J. Meacock

field = startField;

fC = 0;

%Actual simulation
for i = 1:fS.motileSteps
    progressbar(i/fS.motileSteps)
    field = field.stepModel(fS.motiledt);

    if rem(i,fS.FrameSkip) == 0
        outImg = field.drawField();
        pause(0.01)
        
        imPath1 = sprintf(dS.ImgPath,fC);
        fullImPath = [dS.imagedirectory, filesep, imPath1];
        
        imwrite(outImg,fullImPath)
        imshow(outImg,'Parent',ax)
        
        fC = fC + 1;
    end
end
progressbar(1)