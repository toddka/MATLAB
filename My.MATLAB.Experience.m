function My MATLAB Experience 
%% Calculating a particle's path through an electric field.

clear all
%general conditions
tstart=0;
tstop=1.2;
L=1;
d=5;
%Plotting and solving ODEs with different initial conditions
options=odeset('Reltol',.00001,'Event',@eventstop);

init_conditions=[0.5,0,0,1];
[times,sols]=ode45(@particle,[tstart,tstop],init_conditions,options);
figure()
plot(sols(:,2),sols(:,1),'Color',[0,0,1])

hold on
%second conditions
init_conditions2=[1,0,0,1];
[times,sols]=ode45(@particle,[tstart,tstop],init_conditions2,options);

plot(sols(:,2),sols(:,1),'Color',[0,1,0])
hold on
%third conditions
init_conditions3=[1.5,0,0,1];
[times,sols]=ode45(@particle,[tstart,tstop],init_conditions3,options);
plot(sols(:,2),sols(:,1),'Color',[1,0,0])
hold on
%graphing
legend('a','b','c')
xlabel({'z'});
ylabel({'r'});
title({'Particle Path Through Electric Field'});
%event 
    function [eval,stop_flag,edir]=eventstop(t,w)
            R=w(1);theta=w(2);vr=w(3);omega=w(4);
            eval= R;
            stop_flag=1;
            edir=-1;
        
    end
%Creating set-up for ode45 solver
function dwdt=particle(t,w)
       r=w(1); z=w(2); vr=w(3); vz=w(4);
       drdt=vr;
       dzdt=vz;
       if z<L
           c=5;
       else
           c=0;
       end
       dvrdt=-c*r/2;
       dvzdt=c*z;
       dwdt=[drdt;dzdt;dvrdt;dvzdt];
end
%% Calculating reentry conditions (back to Earth) at different angles
clear all
m=8600; %kg
A=2.5^2*pi ;%m
Cd=1 ;%drag coefficient
Re=6371*10^3; %m
d=5*10^3;%m
p0=1.2; %kg/m^
h0=122*10^3 ;%re-entry altitude (m)
v0=8888.89; %re-entry velocity m/s
g=9.8; %m/s^2
tstart=0;
tstop=240*60; %sec

%Phi= 3
phi=3*pi/180;
%Applying the event function
options=odeset('Event',@event);

%Solving with ode45
init_conditions=[Re+h0;0;-v0*sin(phi);v0*cos(phi)/(Re+h0)];
[times,sols]=ode45(@spaceballs,[tstart,tstop],init_conditions,options);

%Separating sols into separate 1 column matrices;
%Converting for x,y direction
for i=1:length(times)
    R=sols(i,1);
    theta=sols(i,2);
    vr=sols(i,3);
    omega=sols(i,4);
    x(i)=R*cos(theta);
    y(i)=R*sin(theta);
%Finding magnitude of acceleration for each time value
    ver(i)=vr;
    vetheta(i)=R*omega;
    magnitudev=sqrt((ver(i))^2+(vetheta(i))^2);
    
    timevalue=times(i);
    differentiatingagain=spaceballs(timevalue,[R,theta,vr,omega]);
    dvrdt=differentiatingagain(3);
    domegadt=differentiatingagain(4);
    
    %plugging into acceleration components in er and theta directions
    er_accel(i)=dvrdt-R*omega^2;
    etheta_accel(i)=R*domegadt+2*vr*omega;
    magnitude3_a(i)=sqrt(er_accel(i)^2+etheta_accel(i)^2);
    altitude3(i)=R-Re;
    times3=times;
end
%Use clearvars to prepare script to solve for next angle
clearvars -except times3 altitude3 magnitude3_a m A Cd Re s p0 h0 v0 g d tstart tstop


%Phi= 4.12
phi=4.12*pi/180;
%Applying the event function
options=odeset('Event',@event);

%Solving with ode45
init_conditions=[Re+h0;0;-v0*sin(phi);v0*cos(phi)/(Re+h0)];
[times,sols]=ode45(@spaceballs,[tstart,tstop],init_conditions,options);

%Separating sols into separate 1 column matrices;
%Converting for x,y direction
for i=1:length(sols)
    R=sols(i,1);
    theta=sols(i,2);
    vr=sols(i,3);
    omega=sols(i,4);
    x(i)=R*cos(theta);
    y(i)=R*sin(theta);
%Finding magnitude of acceleration for each time value
    ver(i)=vr;
    vetheta(i)=R*omega;
    magnitudev=sqrt((ver(i))^2+(vetheta(i))^2);
    
    timevalue=times(i);
    differentiatingagain=spaceballs(timevalue,[R,theta,vr,omega]);
    dvrdt=differentiatingagain(3);
    domegadt=differentiatingagain(4);
    %plugging into acceleration components in er and theta directions
    er_accel(i)=dvrdt-R*omega^2;
    etheta_accel(i)=R*domegadt+2*vr*omega;
    magnitude412_a(i)=sqrt(er_accel(i)^2+etheta_accel(i)^2);
    altitude412(i)=R-Re;
    times412=times;
end
clearvars -except times412 altitude412 magnitude412_a times3 altitude3 magnitude3_a m A Cd Re s p0 h0 v0 g d tstart tstop


%Phi= 5
phi=5*pi/180;
%Applying the event function
options=odeset('Event',@event);

%Solving with ode45
init_conditions=[Re+h0;0;-v0*sin(phi);v0*cos(phi)/(Re+h0)];
[times,sols]=ode45(@spaceballs,[tstart,tstop],init_conditions,options);

%Separating sols into separate 1 column matrices;
%Converting for x,y direction
for i=1:length(times)
    R=sols(i,1);
    theta=sols(i,2);
    vr=sols(i,3);
    omega=sols(i,4);
    x(i)=R*cos(theta);
    y(i)=R*sin(theta);
%Finding magnitude of acceleration for each time value
    ver(i)=vr;
    vetheta(i)=R*omega;
    magnitudev=sqrt((ver(i))^2+(vetheta(i))^2);
    
    timevalue=times(i);
    differentiatingagain=spaceballs(timevalue,[R,theta,vr,omega]);
    dvrdt=differentiatingagain(3);
    domegadt=differentiatingagain(4);
    %plugging into acceleration components in er and theta directions
    er_accel(i)=dvrdt-R*omega^2;
    etheta_accel(i)=R*domegadt+2*vr*omega;
    magnitude5_a(i)=sqrt(er_accel(i)^2+etheta_accel(i)^2);
    altitude5(i)=R-Re;
    times5=times;
end
clearvars -except times5 altitude5 magnitude5_a times412 altitude412 magnitude412_a times3 altitude3 magnitude3_a m A Cd Re s p0 h0 v0 g d tstart tstop

%Phi= 9
phi=9*pi/180;
%Applying the event function
options=odeset('RelTol',.000001,'Event',@event);

%Solving with ode45
init_conditions=[Re+h0;0;-v0*sin(phi);v0*cos(phi)/(Re+h0)];
[times,sols]=ode45(@spaceballs,[tstart,tstop],init_conditions,options);
sols
%Separating sols into separate 1 column matrices;
%Converting for x,y direction
for i=1:length(times)
    R=sols(i,1);
    theta=sols(i,2);
    vr=sols(i,3);
    omega=sols(i,4);
    x(i)=R*cos(theta);
    y(i)=R*sin(theta);
%Finding magnitude of acceleration for each time value
    ver(i)=vr;
    vetheta(i)=R*omega;
    magnitudev=sqrt((ver(i))^2+(vetheta(i))^2);
    
    timevalue=times(i);
    differentiatingagain=spaceballs(timevalue,[R,theta,vr,omega]);
    dvrdt=differentiatingagain(3);
    domegadt=differentiatingagain(4);
    %plugging into acceleration components in er and theta directions
    er_accel(i)=dvrdt-R*(omega)^2;
    etheta_accel(i)=R*domegadt+2*vr*omega;
    magnitude9_a(i)=sqrt((er_accel(i))^2+(etheta_accel(i))^2);
    altitude9(i)=R-Re;
    times9=times;
end
clearvars -except  times9 altitude9 magnitude9_a times5 altitude5 magnitude5_a times412 altitude412 magnitude412_a times3 altitude3 magnitude3_a m A Cd Re s p0 h0 v0 g d tstart tstop


%Plotting Accelerations
figure()
    plot(times3,magnitude3_a)
        title('Magnitude of Acceleration vs Time for 3 Degrees')
        xlabel('time')
        ylabel('magnitude of acceleration')
figure()
plot(times412,magnitude412_a)
        title('Magnitude of Acceleration vs Time for 4.12 Degrees')
        xlabel('time')
        ylabel('magnitude of acceleration')
figure()
plot(times5,magnitude5_a)
        title('Magnitude of Acceleration vs Time for 5 Degrees')
        xlabel('time')
        ylabel('magnitude of acceleration')
figure()
plot(times9,magnitude9_a)
        title('Magnitude of Acceleration vs Time for 9 Degrees')
        xlabel('time')
        ylabel('magnitude of acceleration')
%Plotting Altitudes
figure()
plot(times3,altitude3)
        title('Altitude vs Time for 3 Degrees')
        xlabel('time')
        ylabel('altitude')
figure()
plot(times412,altitude412)
        title('Altitude vs Time for 4.12 Degrees')
        xlabel('time')
        ylabel('altitude')
figure()
plot(times5,altitude5)
        title('Altitude vs Time for 5 Degrees')
        xlabel('time')
        ylabel('altitude')
figure()
plot(times9,altitude9)
        title('Altitude vs Time for 9 Degrees')
        xlabel('time')
        ylabel('altitude')
function dwdt=spaceballs(t,w)
        R=w(1);theta=w(2);vr=w(3);omega=w(4);
        B=g*Re^2/R^2;
        p=p0*exp(-(R-Re)/d);
        V=sqrt(vr^2+(R*omega)^2);
        C=1/2*p0*Cd*A*V/m*exp(-(R-Re)/d);
        
        dRdt=vr;
        dthetadt=omega;
        dvrdt=-B+R*omega^2-C*vr;
        domegadt=-2*vr*omega/R-C*omega;
        dwdt=[dRdt;dthetadt;dvrdt;domegadt];
    end
function [eval,stop_flag,edir]=event(t,w)
            R=w(1);theta=w(2);vr=w(3);omega=w(4);
            eval=R-Re;
            stop_flag=1;
            edir=-1;
end

%% Calculating a quadcopters path based on accelerometer data
clear all
m=.01; %kilograms
g=9.81 ;%m/s^2
L=.02 ;%meters
c=.0012; %Ns^2/m^2
rho=1.2 ;%kg/m^2
tstart=0; %secs
tstop=20;
init_conditions=[0,0,0];
h=1;
dhdt=0;
Kp=1;
Kd=0;
Ki=0;

%Solving ODE numerical and plotting A, B, and C
%A
[times,sols]=ode45(@quadcopter,[tstart,tstop],init_conditions);
figure()
plot(times,sols(:,1),'Color',[0,0,1])
hold on
Kp=1;
Kd=.2;
Ki=0;
%B
[times,sols]=ode45(@quadcopter,[tstart,tstop],init_conditions);
plot(times,sols(:,1),'Color',[1,0,0])
hold on
%C
Kp=1;
Kd=0.2;
Ki=0.2;
%C
[times,sols]=ode45(@quadcopter,[tstart,tstop],init_conditions);

plot(times,sols(:,1),'Color',[0,1,0])

legend('a','b','c')
xlabel({'Time (s)'});
% Create ylabel
ylabel({'Height (m)'});
% Create title
title({'Quadcopter Vertical Path'});
    function dwdt=quadcopter(t,w)
       y=w(1); vy=w(2); Q=w(3);
       P=Kp*(h-y)+Kd*(dhdt-vy)+Ki*Q;
       C=2*m*sqrt(m*g/(pi*rho*L^2));
       dydt=vy;
       dvydt=8*P/(vy*m+C)-(c*vy^2/m)-g;
       dQdt=h-y;
       dwdt=[dydt;dvydt;dQdt];
    end
end