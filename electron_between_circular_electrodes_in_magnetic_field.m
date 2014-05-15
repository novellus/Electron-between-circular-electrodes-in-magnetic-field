clear
clc

R=.025;
R0=.005;
V=30000;
I=1;
N=5000;
Rb=.06; %Radius of magnetic inducing coil.

dt=.000000000001;

u0=4*pi*10^-7;
q=-1.6021766*10^(-19);
m=9.10938*10^(-31);

y=0;
vy=0;
x=R0;
vx=0;

xg=[x];
yg=[y];

for i=1:10000,
    x=x+vx*dt;
    y=y+vy*dt;
    xg=[xg x];
    yg=[yg y];
    if(sqrt(x^2+y^2)>R || sqrt(x^2+y^2)<R0)
        break;
    end
    [K,E]=ellipke(sqrt(4*sqrt(x^2+y^2)/(Rb*(1+sqrt(x^2+y^2)/Rb)^2)));
    B=I*u0*N/(2*pi*Rb*(1+sqrt(x^2+y^2)/Rb))*(K+E*(1-(x^2+y^2)/Rb^2)/((1+sqrt(x^2+y^2)/Rb)^2-4*sqrt(x^2+y^2)/Rb));
    ax=q/m*(-V*x/(log(R/R0)*(x^2+y^2))+B*vy);
    ay=q/m*(-V*y/(log(R/R0)*(x^2+y^2))-B*vx);
    vx=vx+ax*dt;
    vy=vy+ay*dt;
end

c1x=[];
c1y=[];
c2x=[];
c2y=[];

for j=1:1000,
    c1x=[c1x R0*cos(2*pi*j/1000)];
    c1y=[c1y R0*sin(2*pi*j/1000)];
    c2x=[c2x R*cos(2*pi*j/1000)];
    c2y=[c2y R*sin(2*pi*j/1000)];
end
plot(xg,yg,c1x,c1y,c2x,c2y)
axis equal