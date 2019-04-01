function [fPar,fPerp,fRot] = calcGeomFactors(a)
%CALCGEOMFACTORS calculates the geometrical factors associated with the
%calculation of friction tensors. Based on the sedimentation dynamics of
%DNA fragments with variable lengths - for full details see following
%publications:
%
% Wensink, H. H., & Löwen, H. (2012). Emergent states in dense systems of 
% active rods: From swarming to turbulence. Journal of Physics Condensed 
% Matter, 24(46). https://doi.org/10.1088/0953-8984/24/46/464130
%
% Tirado, M. M., Martínez, C. L., & de la Torre, J. G. (2003). Comparison of
% theories for the translational and rotational diffusion coefficients of 
% rod?like macromolecules. Application to short DNA fragments. The Journal 
% of Chemical Physics, 81(4), 2047–2052. https://doi.org/10.1063/1.447827
%
%   INPUTS:
%       -a: Aspect ratios of rods. Nx1 vector
%
%   OUTPUTS:
%       -fPar: Frictional factor associated with rod moevement parallel to
%       its long axis
%       -fPerp: Frictional factor associated with rod movement
%       perpendicular to its long axis
%       -fRot: Frictional factor associated with rod rotation
%
%   Author: Oliver J. Meacock

%Calculates the geometric factors for this cell.
fPar = (2 * pi)./(log(a) - 0.207 + (0.980./a) - (0.133./(a.^2)));
fPerp = (4 * pi)./(log(a) + 0.839 + (0.185./a) + (0.233./(a.^2)));
fRot = ((a.^2) * pi)./(3*(log(a) - 0.662 + (0.917./a) - (0.050./(a.^2))));
