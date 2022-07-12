function vr = rayleighphase(FREQ,VS,BOS,RHO,THK,DV,varargin)
epsilon = 1e-6; 
VS = VS(:);
BOS = BOS(:);
RHO = RHO(:);
THK = THK(:);

VP = VS.*sqrt(2*(1-BOS)./(1-2*BOS));
cn = length(VP);

cvpmin = min(VP); cvpmax = max(VP);
cvsmin = min(VS); cvsmax = max(VS);
vrmin = homogeneous(cvpmin,cvsmin);
vrmin = vrmin/1.05;
vrmax = 1.05*cvsmax;

if nargin == 7, MODE_NUMBER = varargin{1};end  

FN = length(FREQ);
vr = zeros(FN,MODE_NUMBER);
for fi = 1:FN
    if cn == 1
        vr(fi,1) = homogeneous(VP,VS);
    else
        v0 = vrmin; v1 = v0+DV;
        for mi = 1:MODE_NUMBER       
            while v1 < vrmax
                dp1 = FastDelta(FREQ(fi),VS,BOS,RHO,THK,v0);
                if abs(dp1) <= eps
                    vr(fi,mi) = v0;
                    v0 = v0+0.5*DV; v1 = v0+DV;
                    break;
                end
                dp2 = FastDelta(FREQ(fi),VS,BOS,RHO,THK,v1);
                if abs(dp2) <= eps
                    vr(fi,mi) = v1;
                    v0 = v1+0.5*DV; v1 = v0+DV;
                    break;
                end
                if sign(dp1) ~= sign(dp2)
                    vr(fi,mi) = rayleigh_fzero_brent(@FastDelta,FREQ(fi),VS,BOS,RHO,THK,v0,v1); 
                    v0 = vr(fi,mi)+0.5*DV; v1 = v0+DV;
                    break;
                end
                v0=v1;v1=v0+DV;
            end
        end
    end
end
