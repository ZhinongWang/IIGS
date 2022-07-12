function dp = reduced_delta(fq,vs,bos,rho,thk,c)
% get the dispersion function value
%   bos--Poisson's ratio; vs--S-wave velocity;
%   rho--density; thk--thickness; fq--frequency; c--phase velocity

cn = length(vs);
vp = getvp(vs,bos);

x = [0 0 0 0 1];
[v,~] = get_dvdc(vs(cn),vp(cn),rho(cn),c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% begin the recursive process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ic = 1:cn-1
    T = get_T(fq,vs(ic),vp(ic),rho(ic),thk(ic),c);
    x = x*T;
end
dp = real(x*v);

