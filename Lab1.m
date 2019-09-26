%% Info:

% this code is post-data collection analysis for two tubes, one is open and
% one is closed, for more info, go to /Info

%% housekeeping

clear
clc
close all

%% define onstants:

G = 3.75 * 10 ^ 6 ; % psi (shear modulus);


% closed thin wall specimen:

De_closed = 3/4 ; % in, exterior daiamter
t_closed = 1/16 ; % in, thickness
L_closed = 13 ; % inches, measured in lab.
T_max_closed = 400; % (lbs-in), maximum torque
Tau_max_closed = 8620; % psi, max shear stress
FS_closed = mean([ 2.3 2.6]); % given in lab doc as range, so took the mean.



% open thin wall specimen
% everything same as closed, we assume cut is negligible :)

De_open = 3/4 ; % in, exterior daiamter
t_open = 1/16 ; % in, thickness
L_open = 13 ; % inches, measured in lab.
T_max_open = 20; % lbs -in
Tau_max_open = 7800; % psi
FS_open = 2.8; % factor of safety.

%% read data:

addpath('./Data'); % add path of data files

% to use importdata, fopen must be issued first
fopen('400inclosed.txt');
fopen('20inopen.txt');

% now import data safely
Data_closed = importdata('400inclosed.txt'); %import data
Data_open = importdata('20inopen.txt'); %import data


% close open handles
fclose('all');

% extract data

time_closed = Data_closed.data(:,1); % in sec
twist_angle_closed = Data_closed.data(:,2); % in deg
shear_strain_closed = Data_closed.data(:,3); % in deg
Torque_closed = Data_closed.data(:,4); % in-lbf
Axial_closed = Data_closed.data(:,5); % in
twist_angle_closed = twist_angle_closed - twist_angle_closed(1) ; % zero torsional angle

time_open = Data_open.data(:,1); % in sec
twist_angle_open = Data_open.data(:,2); % in deg
shear_strain_open = Data_open.data(:,3); % in deg
Torque_open = Data_open.data(:,4); % in-lbf
Axial_open = Data_open.data(:,5); % in
twist_angle_open = twist_angle_open - twist_angle_open(1) ; % zero total torsional angle


% zero time:
time_closed = time_closed - time_closed(1);
time_open = time_open - time_open(1);


Theortical_rigidity_closed_exact = G .* (1/2) * pi * ( (De_closed/2) - (De_closed - t_closed )).^2; % exact theortical rigidity for closed only

J_open = (1/3) * (2*pi*(De_open/2)) * (t_open)^3 ; % polar moment of inertia
J_Closed = ( 4 * ((pi)*(((De_closed - t_closed)/2)^2))^2 * t_closed ) ./ ( 2 * pi * ((De_closed - t_closed)/2))  ; % polar moment of inertia

Theortical_rigidity_CTW = G .* J_Closed; %theortical rigidity for closed thin wall theory
Theortical_rigidity_OTW = G .* J_open; % theortical rigidity for open thin wall theory


% convert units:

shear_strain_closed = deg2rad(shear_strain_closed);
twist_angle_closed = deg2rad(twist_angle_closed);

shear_strain_open = deg2rad(shear_strain_open);
twist_angle_open = deg2rad(twist_angle_open);



%% backgroun:

%{

No matter what's the theory we're using, it is always the case that we can
estimate the theortical rigitdy by knowing that for both OTW + CTW the rate
of total twist angle (derivative of Phi) is = T / (GJ); so dividing torque
by rate of total twist angle will estimate GJ (rigidity)..

it's also the case that for the exact method (I think it holds for both OTW
and CTW) that the shear strain = ( T*roh ) / (GJ) where roh is radial
distance, so if we have plots for T on Y and Shear strain on X if we take
dslope of that (i.e. T/Y)  * R it'll also be the rigidity :)
%}
%% Closed specimen:

% it's circular, so we can use exact, or thin wall theory.

% Plot the torque vs.  shear strain provided by the extensometer, as well as the torque vs.  shear
% strain calculated using the total rotation angle imposed by the testing machine.

% SS = Shear strain. C = closed
SS_C_Epsilon = shear_strain_closed;  % Shear strain, closed;
SS_C_Twist = (twist_angle_closed ./ L_closed) .* (De_closed/2 );

% Shear strain = SS, O = Open.
SS_O_Epsilon = shear_strain_open;  % Shear strain, closed;
SS_O_Twist = (twist_angle_open ./ L_open) .* (De_open/2 );


figure(1)

plot(Torque_closed,SS_C_Epsilon,'k');
hold on
plot(Torque_closed,SS_C_Twist,'r','LineWidth',2);
xlabel('Torque [lbs-in]');
ylabel('\gamma [unitless]');
title('Shear Strain of a closed circular cross section specimen');
legend('\gamma by extensometer','\gamma by torsional angle')


figure(2)

plot(Torque_open,SS_O_Epsilon,'k');
hold on
plot(Torque_open,SS_O_Twist,'r','LineWidth',2);
xlabel('Torque [lbs-in]');
ylabel('\gamma [unitless]');
title('Shear Strain of a closed circular cross section specimen');
legend('\gamma by extensometer','\gamma by torsional angle')


% do line fit:

[ p S ] = polyfit(SS_C_Epsilon,Torque_closed,1);
Rigidity_C_Epsilon = abs(p(1))*(De_closed/2);

[ p S ] = polyfit(SS_C_Twist,Torque_closed,1);
Rigidity_C_Twist = abs(p(1))*(De_closed/2);

[ p S ] = polyfit(SS_O_Epsilon,Torque_open,1);
Rigidity_O_Epsilon = abs(p(1))*(De_open/2);

[ p S ] = polyfit(SS_O_Twist,Torque_open,1);
Rigidity_O_Twist = abs(p(1))*(De_open/2);

%%