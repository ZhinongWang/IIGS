function dpvalue = FastDelta( vp,vs,rou,d,fre,c )

w=2*pi*fre; 
k=w/c;
layernum=length(vs);
r=zeros(layernum,1);
rr=zeros(layernum,1);
s=zeros(layernum,1);
t=zeros(layernum,1);
miu=zeros(layernum,1);
Ca=zeros(layernum-1,1);
Sa=zeros(layernum-1,1);
Cb=zeros(layernum-1,1);
Sb=zeros(layernum-1,1);
%%
for ic=1:layernum
    r0=1-c*c/(vp(ic)*vp(ic));
    s0=1-c*c/(vs(ic)*vs(ic));
    if real(c)<=vp(ic)
        r(ic)=sqrt(r0);
    else
        r(ic)=1i*sqrt(-r0);
    end
    
    if real(c)<=vs(ic)
        s(ic)=sqrt(s0);
    else
        s(ic)=1i*sqrt(-s0);
    end
end
for ic=1:layernum-1
    r0=1-c*c/(vp(ic)*vp(ic));
    s0=1-c*c/(vs(ic)*vs(ic));
    if c<=vp(ic)
        Ca(ic)=cosh(k*r(ic)*d(ic));
        Sa(ic)=sinh(k*r(ic)*d(ic));
    else
        Ca(ic)=cos(k*sqrt(-r0)*d(ic));
        Sa(ic)=1i*sin(k*sqrt(-r0)*d(ic));
    end
    
    if c<=vs(ic)
        Cb(ic)=cosh(k*s(ic)*d(ic));
        Sb(ic)=sinh(k*s(ic)*d(ic));
    else
        Cb(ic)=cos(k*sqrt(-s0)*d(ic));
        Sb(ic)=1i*sin(k*sqrt(-s0)*d(ic));
    end
end
%%
t=2-c^2./(vs.^2);
miu=rou.*vs.^2;
rr=vs.^2/c^2;

X=zeros(1,6);
X(1,1)=2*t(1); X(1,2)= -t(1)^2; X(1,5)= -4; X(1,6)=2*t(1);
X=miu(1)^2*X;
for ic=1:layernum-1
    e=rou(ic+1)/rou(ic);
    h=2*(rr(ic)-e*rr(ic+1));
    a=e+h; a1=a-1;
    b=1-h; b1=b-1;
    
    p1=Cb(ic)*X(1,2)+s(ic)*Sb(ic)*X(1,3);
    p2=Cb(ic)*X(1,4)+s(ic)*Sb(ic)*X(1,5);
    p3=Sb(ic)*X(1,2)/s(ic)+Cb(ic)*X(1,3);
    p4=Sb(ic)*X(1,4)/s(ic)+Cb(ic)*X(1,5);
    
    q1=Ca(ic)*p1-r(ic)*Sa(ic)*p2;
    q2= -Sa(ic)*p3/r(ic)+Ca(ic)*p4;
    q3=Ca(ic)*p3-r(ic)*Sa(ic)*p4;
    q4= -Sa(ic)*p1/r(ic)+Ca(ic)*p2;
    
    y1=a1*X(1,1)+a*q1;
    y2=a*X(1,1)+a1*q2;
    z1=b*X(1,1)+b1*q1;
    z2=b1*X(1,1)+b*q2;
    
    X(1,1)=b1*y1+b*y2;
    X(1,2)=a*y1+a1*y2;
    X(1,3)=e*q3;
    X(1,4)=e*q4;
    X(1,5)=b1*z1+b*z2;
    X(1,6)=a*z1+a1*z2;
end
dpvalue=X(1,2)+s(layernum)*X(1,3)-r(layernum)*(X(1,4)+s(layernum)*X(1,5));
end

