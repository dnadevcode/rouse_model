function tests = rouseTest
    tests = functiontests(localfunctions);
end


function testOne(testCase)
    hfig = rouse_gui('off');

        
    
    verifyEqual(testCase,hfig.Name,'Rouse model GUI v0.1');
%     close hfig;

end


function testMex(testCase)
    mex rouse.c;

end




function testTwo(testCase)

Nb = 10;

zeta = 1.33*10^(-8); % Friction coefficient of a bead, Ns/m
ks = 2.88*10^(-6); % Spring constant N/m
Xs = 0.5400e-06; % equilibrium rest length of the spring m
kT = 4.11*10^(-21); % J
h=0.001; % timestep

numFrames = 200;
frameLen = 0.04;       % time difference (s) between two frames
fTime = numFrames*frameLen; % simulation stop time


[sol2] = rouse(Xs, Nb, zeta, ks, kT, h, fTime ,frameLen,numFrames)';




        
    
    verifyEqual(testCase,size(sol2),[200 10]);
%     close hfig;

end