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

%%

