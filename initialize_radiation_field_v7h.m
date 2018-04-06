% Initialize radiation field
radfield = zeros(param.Nsnap,param.nslices);
if param.fieldprofile
    dt = param.zsep*param.lambda0/c;
    tvector = [1:param.nslices].*dt-round((param.nslices+param.Nsnap)/2).*dt;
    radfield(1,:) = param.E0.*exp(-tvector.^2/2/param.sigmatfield^2);
else
    radfield(1,:)=param.E0;
end