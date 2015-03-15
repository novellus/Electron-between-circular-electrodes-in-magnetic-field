close all
clear
clc

R=.025; %Radius of outer electrode
R0=.005; %Radius of inner electrode
V=30000; %Differential voltage
I=1; %Magnetic coil current
N=5000; %Number of turns in magnetic coil
Rb=.06; %Radius of magnetic inducing coil.

dt=.000000000001; %Timestep in secs (electrons are nearly relativistic at these voltages)

%physical constants
u0=4*pi*10^(-7); %permeability of free space
e0=8.854187817*10^(-12); %permittivity of free space
q=-1.6021766*10^(-19); %electron charge
m=9.10938*10^(-31); %electron mass
c=299792458 ; %Speed of light

neps=(8*12/V)/-q; %Total electrons/sec of the arc, as approximated by power input and voltage output.

%Electron initial position and velocity
y=0;
vy=0;
x=R0;
vx=0;

%Position, time, power of emmitted radiation (single electron), and velocity as a fraction of the speed of light
xg=[x];
yg=[y];
tg=[0];
pg=[0];
vcg=[sqrt(vx^2+vy^2)/299792458];

for i=1:10000,
    [K,E]=ellipke(sqrt(4*sqrt(x^2+y^2)/(Rb*(1+sqrt(x^2+y^2)/Rb)^2)));
    B=I*u0*N/(2*pi*Rb*(1+sqrt(x^2+y^2)/Rb))*(K+E*(1-(x^2+y^2)/Rb^2)/((1+sqrt(x^2+y^2)/Rb)^2-4*sqrt(x^2+y^2)/Rb));
    ax=q/m*(-V*x/(log(R/R0)*(x^2+y^2))+B*vy);
    ay=q/m*(-V*y/(log(R/R0)*(x^2+y^2))-B*vx);
    vx=vx+ax*dt;
    vy=vy+ay*dt;
    x=x+vx*dt;
    y=y+vy*dt;
    tg=[tg dt+tg(end)];
    xg=[xg x];
    yg=[yg y];
    pg=[pg q^2*(ax^2+ay^2)/(6*pi*e0*c^3)];
    vcg=[vcg sqrt(vx^2+vy^2)/299792458];
    if(sqrt(x^2+y^2)>R || sqrt(x^2+y^2)<R0) %Stop if the electron collides with an electrode
        break;
    end
end

%Plot the electrodes for visual aids
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

figure('Position',[1200,620,700,500]);
plot(xg,yg,c1x,c1y,c2x,c2y);
title('position (m)');
axis equal;

figure('Position',[1200,50,700,500]);
s1=subplot(2,1,1);
plot(tg,pg);
xlabel('time (sec)');
ylabel('power radiated/e- (watts)');

s2=subplot(2,1,2);
plot(tg,vcg);
xlabel('time (sec)');
ylabel('v/c');
linkaxes([s1 s2],'x');

totalPowerRadiated_Watts=sum(dt*pg)*neps
exposureTimeToOneYearRadLimit_days=.05*70/totalPowerRadiated_Watts/60/60/24 %50mSv (Sv=J/kg) is US rad worker yearly limit
