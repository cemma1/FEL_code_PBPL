function [ sep ] = separatrix( psi,psir )
% Calculate height of separatrix at phase psi for resonant phase psir
% Detailed explanation goes here
sep = sqrt(abs(cos(psi)+cos(psir)+(psi+psir-pi*sign(psir))*sin(psir)));
end

