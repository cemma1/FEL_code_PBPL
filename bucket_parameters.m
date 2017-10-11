function [ psi1, psi2, bucket_height, capture_fraction, bucket_area, bunching, sinaverage] = bucket_parameters(psir)
% Calculate various bucket parameters as a function of psir
 
psi1 = pi - psir;
pond = @(psi,psir) cos(psi1)+psi1*sin(psir)-(cos(psi)+psi*sin(psir));
psi2 = fsolve(@(psi) pond(psi,psir),[-pi:psi1*0.99]);
psi2 = psi2(1);
bucket_height = sqrt(cos(psir)-(pi/2-psir)*sin(psir));
bucket_area = (1-sin(psir)) / (1+sin(psir));
capture_fraction = (psi1 - psi2)/2/pi;
potential = @(psi,psir) cos(psi)+psi*sin(psir);
bucket_sep = @(psi) sqrt( (-potential(psi1,psir)+potential(psi,psir))/2);
bucket_bun = @(psi) bucket_sep(psi).*exp(1i*psi);
bucket_sin = @(psi) bucket_sep(psi).*sin(psi);
bunching = abs(integral(bucket_bun,psi2,psi1) / integral(bucket_sep,psi2,psi1));
sinaverage = abs(integral(bucket_sin,psi2,psi1) / integral(bucket_sep,psi2,psi1));
end