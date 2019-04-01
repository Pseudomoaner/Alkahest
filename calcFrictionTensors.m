function [fT,fR] = calcFrictionTensors(a,u,f0)
%CALCFRICTIONTENSORS calculates the friction tensors associated with the
%rods with input geometry given by a and u.
%
%   INPUTS:
%       -a: Aspect ratios of rods. Nx1 vector
%       -u: Unit orientation vector of each rod. Nx2 matrix
%       -f0: Stokesian friction coefficient
%
%   OUTPUTS:
%       -fT: Translational friction tensor
%       -fR: Rotational friction tensor
%
%   Author: Oliver J. Meacock

[fPar,fPerp,fRot] = calcGeomFactors(a);
fR = f0 * fRot;
fT = zeros(2,2,size(a,1));
for i = 1:size(a,1)
    uSq = u(i,:)'*u(i,:);
    fT(:,:,i) = f0 * (fPar(i)*uSq + fPerp(i)*(eye(2) - uSq)); %Sets the basis of the different parallel and perpendicular friction coefficients to be oriented along the cell axis.
end