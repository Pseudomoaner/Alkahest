function field = simulateWensinkFieldBurnIn(startField,fS,compiled)
%SIMULATEWENSINKFIELDBURNIN performs initial, short timestep simulations of
%the given input field to allow the elastic strain introduced by the
%initialisation steps to be released without causing the simulation to
%become numerically unstable. Results are not plotted.
%
%   INPUTS:
%       -startField: a WensinkField object
%       -fS: field settings associated with startField
%       -compiled: whether the energy gradient calculation function has
%       been compiled on this system
%
%   OUTPUTS:
%       -field: a WensinkField object
%
%   Author: Oliver J. Meacock

field = startField;

for i = 1:fS.burnInSteps
    progressbar(i/fS.burnInSteps)
    
    field = field.stepModel(fS.burnIndt,fS.f0,compiled);
end
progressbar(1)