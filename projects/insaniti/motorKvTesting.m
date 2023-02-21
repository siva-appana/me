%% Simulation of DC Motor Systems
% By Varun Madabushi
% RoboJackets
% March 2021

close all

%% Weapon/Shaft Parameters
weaponLink = 'https://www.e-fliterc.com/product/power-10-brushless-outrunner-motor-1100kv-3.5mm-bullet/EFLM4010A.html'; % -Change
moment = 1.846/3417.114418; % Kg*m^2 (Input lb*in^2 to the left of divide)
motorDiameter = 1.4; %Input: Inches of Motor Diameter
pulleyMA = (2.25)/motorDiameter; %(Input: Inches of Pully Diameter)
J = moment * (1/pulleyMA); % Kg-m/s^2

%% DC Motor Parameters
% Based on: https://www.scorpionsystem.com/catalog/aeroplane/motors_1/s-30_v2/SII-3020-1110KV/
R = 0.043; % Ohms -Change
Kv = 1100 * (pi / 30); % Rad/s/Volt -Change
Ki = 1/Kv; % N-m/A
D = 0; % Viscous friction coefficient set to 0 bc hard to find

%% Initialize Arrays Holding State Variables
dt = 0.0001;
maxT = 5; % # of seconds
t = 0:dt:maxT;
theta = zeros(1,length(t));
omega = zeros(1,length(t));
current = zeros(1,length(t));
rpm = zeros(1,length(t));
KE = zeros(1,length(t));
V = 11.1; % Voltage of battery (Parameter)
r = Kv * V * (30 / pi);

%% Forward Euler Simulation
for k = 1:length(t)-1
    
    % Evaluate derivatives of all state variables
    % Speed is the integral of position
    thetadot = omega(k);
    
    % Current from Ohm's law
    current = (V - (omega(k)/Kv))/R;
    
    % Acceleration (omegadot)
    % Viscous friction model: friction torque proportional to speed
    omegadot = (1/J)*(Ki*current - D*omega(k)); % tau_L = D*omega
    
    % RPM
    rpm(k) = omega(k) * 30 / pi;
    
    % Kinetic Energy
    KE(k) = 0.5 * moment * omega(k)^2;
   
    % Forward Euler Method (First Order Taylor Series)
    theta(k+1) = theta(k) + dt*thetadot;
    omega(k+1) = omega(k) + dt*omegadot;
end

%% Plots
subplot(2,1,1)
plot(t,theta)
ylabel('Theta')
subplot(2,1,2)
hold on
plot(t,rpm)
plot(t, r*ones(1,length(t)), 'r--')
xlabel('Time')
ylabel('RPM')

%% Write to Excel
cA = readcell('Siva_WeaponMotors.xlsx');
for n = 1:length(t)-1
   if rpm(n) <= r - 1000
       spinUpTime = n * 0.0001;
   end
end
for n = 1:length(t)*0.67
   if rpm(n) <= r*0.67
       time67 = n * 0.0001;
   end
end
for n = 1:length(t)*0.33
   if rpm(n) <= r*0.33
       time33 = n * 0.0001;
   end
end
cA = [cA; {weaponLink},{pulleyMA},{moment},{r},{max(KE)},{r*0.67},{time67},{r*0.33},{time33},{spinUpTime}];
writecell(cA, 'Siva_WeaponMotors.xlsx');
