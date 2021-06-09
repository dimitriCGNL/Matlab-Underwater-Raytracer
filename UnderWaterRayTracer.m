%% UNDERWATER RAT TRACER.

% Dimitri Croes
% Engineering physics student at Eindhoven, The Netherlands
% 2021

% paper used for making this tool
% TU Delft, “Underwater propagation Ray acoustics,” pp. 65–80.



%% Clear
clear
close all
clc
%% Setup SSP


%Depth, Positive is down
z=[0 102 204 306 408 537 667 796 926 1056 1185 1315 1444 1574 1704 1833 1963 2092 2222 2352 2481 2611 2740]/3.281;
%SSP
c=[5082.99 5084.57 5086.16 5087.75 5089.33 5089.05 5088.88 5088.84 5088.91 5089.09 5089.40 5089.82 5090.36 5091.02 5091.80 5092.69 5093.71 5094.84 5096.09 5097.45 5098.94 5100.54 5102.26]/3.281;

%interpolation of SSP
zz = 0:1:round(7300/3.281);
cc=interp1(z,c,zz,'pchip','extrap');

%% initial values
xmax=100000; %max value of x.
dx=50; %step size of x (lower = longer calculations and bigger files)
x0=0; %start x position
z0=300; %start depth positive is down
theta0=linspace(-4,4,20); %starting angle (can be array of angles)
Ylimit = [0 zz(end)];

%% run Simulation
arch = computer('arch');

if zz(1) ~= 0 %extend SSP to surface if that isn't already the case.
    zz = [0 zz(:)'];
    cc = [cc(1) cc(:)'];
end


g = diff(cc)./diff(zz); %Calculate gradient for every depth of zz


for m=1:length(theta0) % loops through starting angles
    %calculate first values by "Hand"
    x(1)=x0;
    Z(1)=z0;
    
    angle(1)=theta0(m);
    c0(1)=cc(Z(1));
    G=g(Z(1));
    R0(1)=-c0(1)/(cosd(angle(1))*G);
    
    i=2;
    while x<=xmax
        x(i)=x(i-1)+dx;
        angle(i)=asind(((x(i)-x(i-1))/R0(i-1))+sind(angle(i-1))); %new angle
        
        
        Z(i)=(R0(i-1)*(cosd(angle(i))-cosd(angle(i-1))))+Z(i-1); %new depth
        if Z(i)<=0 %if new depth is above surface (or out of bounds)
            Z(i)=0; %set depth to surface
            angle(i)=-angle(i-1); %reflect ray
            c0(i)=cc(1); %set speed to surface speed
            G=g(1); %set gradient to surface gradient
        elseif Z(i)>=zz(end) %if new depth is under seafloor (or out of bounds)
            Z(i)=zz(end);
            angle(i)=-angle(i-1);
            c0(i)=cc(end);
            G=g(end);
            
        else
            c0(i)=cc(floor(Z(i)+1))+(Z(i)+1-floor(Z(i)+1))*((cc(ceil(Z(i)+1))-cc(floor(Z(i)+1)))/(ceil(Z(i)+1)-floor(Z(i)+1))); %lineair interpolation of speed at current depth
            
            
            G=g(round(Z(i)+1)); %gradient at current depth (should also be interpolated like above but im lazy)
            
        end
        
        
        
        
        R0(i)=c0(i)/(cosd(angle(i))*G); %calculate curvature radius
        
        
        i=i+1;
        
    end
    X{m}=x; %store variables to cell
    ZZ{m}=Z;
    clear x Z R0 angle %clean variables for new run
end
%% plotting
f = figure(1);
if strcmp(arch,'maci')
    set(f,'Position',[10 1500 1300 400])
else
    set(f,'Position',[10 100 1300 400])
end
subplot(1,5,1)
plot(cc,zz);
axis ij
grid on
title('Sound Speed Profile')
xlabel('Sound Speed, m/s')
ylabel('Depth, m')
ylim(Ylimit)
subplot(1,5,2:5)


for m=1:length(theta0)
    plot(X{m}/1000,ZZ{m})
    hold on
end
hold off
xlim([0 xmax/1000])
ylim(Ylimit)
axis ij
grid on
title('Ray Trace, Start Depth: '+string(z0)+'m')
xlabel('Distance, km')
ylabel('Depth, m')
legend('Start Angle: '+string(theta0),'Location','southeast')