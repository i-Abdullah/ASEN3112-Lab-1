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

%addpath('./Data'); % add path of data files

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
Torsinal_angle_closed = Data_closed.data(:,2); % in deg
Epsilon_closed = Data_closed.data(:,3); % in deg
Torque_closed = Data_closed.data(:,4); % in-lbf
Axial_closed = Data_closed.data(:,5); % in


time_open = Data_open.data(:,1); % in sec
Torsinal_angle_open = Data_open.data(:,2); % in deg
Epsilon_open = Data_open.data(:,3); % in deg
Torque_open = Data_open.data(:,4); % in-lbf
Axial_open = Data_open.data(:,5); % in


%%