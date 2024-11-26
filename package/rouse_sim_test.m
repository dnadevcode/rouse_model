
addpath(genpath('C:\Users\Lenovo\git\dnadevcode\rouse_beads_project\package'));
%%
% Parameters
% Nb number of beads
Nb = 10;

zeta = 1.33*10^(-8); % Friction coefficient of a bead, Ns/m
ks = 2.88*10^(-6); % Spring constant N/m
Xs = 0.5400e-06; % equilibrium rest length of the spring m
kT = 4.11*10^(-21); % J
h=0.001; % timestep
% demo_flag = false;

demo_flag = true;


numFrames = 200;
frameLen = 0.04;       % time difference (s) between two frames
fTime = numFrames*frameLen; % simulation stop time

tic
[sol] = rouse_sim(Xs, Nb, zeta, ks, kT, h, fTime ,frameLen,numFrames);
toc

% mex rouseMex.cpp

% sol = rouseMex(Xs, Nb, zeta, ks, kT, h, fTime ,frameLen,numFrames)
%


figure,plot(sol,'black')

%%
% mex rouse.c
Nb=30;
tic
[sol2] = rouse(Xs, Nb, zeta, ks, kT, h, fTime ,frameLen,numFrames)';
toc
figure;
plot(sol2,'red')
