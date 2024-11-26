
function [sol] = rouse_sim(Xs, Nb, zeta, ks, kT, h, fTime ,frameLen,numFrames)
    % https://aip.scitation.org/doi/full/10.1063/1.4907552?casa_token=I-2ghjXaL90AAAAA:L9pax9TfDCxIOhst1-TjbVooI3Eh3Vq2bFYiaLOaEJPkE9-YdwND55H3zTWSzUEUXZdGn2_hbuY
    % from the SI
    % Modeling the relaxation of internal DNA segments
    % during genome mapping in nanochannel
    % Brownian dynamics simulation

    if nargin < 1
        % Nb number of beads
        Nb = 25;
        zeta = 1.33*10^(-8); % friction coefficient of a bead, Ns/m    
        ks = 2.88*10^(-6); % spring constant N/m
        Xs = 0.5400e-06;
        kT = 4.11*10^(-21); % J

        h=0.001; % timestep     
        demo_flag = false;
        numFrames=200;
        
        fTime = numFrames*frameLen;

    end
    %% input parameters
    steps = 0:h:fTime;
    
    % solution output: only the frames we want..
    sol = zeros(numFrames,Nb);
    
    rmatCur = linspace(0,(Nb-1)*Xs,Nb)'; % bead positions
    
    solFrames = floor(linspace(1,frameLen*numFrames/h,numFrames));
    valueFrames = zeros(length(steps),1);
    valueFrames(solFrames) = 1;
    
%     rmat = zeros(length(steps),Nb+2);
%     rmat(1,:) = linspace(-Xs,(Nb)*Xs,Nb+2); % first row..
%     
    
    %% random force acting on bead i, i.e. white noise
    M = zeros(1,Nb);                          % Mean vector 
%     var = h*ones(1,length(M));                 % Variance vector
%     Cov = var.*eye(Nb,Nb);                  % Covariance matrix / should make more mem. efficient
%     x = mvnrnd(M,Cov,length(steps));                  % MultiVariate random vector added to the input

    const = zeros(Nb,1);
    const(1)=-Xs;
    const(end)=Xs;
    sol(1,:) = rmatCur;
%     rmatNext = zeros(1,Nb);
%     rtemp = zeros(1,Nb);
    
    % create rmat with coefficients in the base R1,...,RBb
    rMat = zeros(Nb,Nb);
    rMat(1,1:2) = [-1 1];
    rMat(end,end-1:end) = [1 -1];
    for j=2:size(rMat,1)-1
        rMat(j,j-1:j+1) = [1 -2 1];
    end

    k=2;
    % Brownian Dynamics simulation corresponding to the Rouse-like model
    % S-23, gives stochastic Euler differential equation
    for t=2:length(steps)
%         t
%         rmatCur(1) =  rmatCur(2)-Xs;
%         rmatCur(end) =  rmatCur(end-1)+Xs;
        % go through all the steps
        curRnd = sqrt(2*kT/(zeta*h))*randn(Nb,1);
        % first step: Euler method
        %https://sites.me.ucsb.edu/~moehlis/APC591/tutorials/tutorial7/node2.html
        % Amin D. Computational and theoretical modelling of self-healable polymer materials (Doctoral dissertation, University of Reading).
        % because we simulate Wiener process, there is a sqrt(h) to relate
        % to normal dist
        f1 = (ks/zeta*(rMat*rmatCur+const)+curRnd);
        initialGuess = rmatCur+h*f1;
%         rmatNext = initialGuess;
%         sum((sol(1,:) -rmatNext').^2)
        % now use initial guess
        f2 = (ks/zeta*(rMat*initialGuess+const)+curRnd);
%         % predictor/corrector update:
        rmatNext = rmatCur+1/2*h*(f1+f2);
%           plot_match_iterate({20+rmatCur*10^9/117,20+initialGuess*10^9/117,20+rmatNext*10^9/117})
%                   plot_match(20+sol(1,:)*10^9/117,20+rmatNext*10^9/117)


       if valueFrames(t)==1
           sol(k,:) = rmatNext;
           k = k+1;
       end
       rmatCur = rmatNext;
    end
    

    
%     simSol = rmat(:,2:end-1);
    
    % take only the relevant timeframes, so h has to be small for this to
    % be long enough
%     sol = simSol(round(linspace(1,frameLen*numFrames/h,numFrames)),:);
  
    
end
