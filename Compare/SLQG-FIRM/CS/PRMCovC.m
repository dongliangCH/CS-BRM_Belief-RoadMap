function [CovCost, CollisonCost, K] = PRMCovC(PhatPrior0, PtildePrior0, PhatPriorf, PtildePriorf, N, Signal, param, Xbar, V, Cspace)

% Notation: Ak is single step of matrix, AA is 3d array with each Ak
% stacked in the 3rd dimension, A is block matrix

%% Problem setup

% Number of steps N
% Choose time step dt = 0.2

% Dynamics
dt = param.dt;
nx = 4;
nu = 2;
nw = 4;
Ak = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
Bk = [dt ^ 2 / 2 0; 0 dt ^ 2 / 2; dt 0; 0 dt];
% Gk = 0.1 * eye(nx);
Gk = diag([0.05 0.08 0.05 0.05]);

% Sequences of system matricies
AA = repmat(Ak, [1, 1, N]);
BB = repmat(Bk, [1, 1, N]);
GG = repmat(Gk, [1, 1, N]);

% Make block matricies for dynamics
A = makeStateMatrix(AA, N);
B = makeInputMatrix(AA, BB, N);

% Cost
Qk = 2 * blkdiag(2, 2, 2, 2);
Rk = 1 * blkdiag(2, 2);

% Make block matricies for cost
Q = [];
R = [];
for i = 1:N
    Q = blkdiag(Q, Qk);
    R = blkdiag(R, Rk);
end
Q = blkdiag(Q, zeros(nx));


%% Kalman filter
% Observation model

[SignalNum, ~] = size(Signal); 
ny = (nx-2) * SignalNum + 2;
Ck = repmat([1 0 0 0; 0 1 0 0], [SignalNum, 1]);
Ck = [Ck; 0 0 1 0; 0 0 0 1]; 
CC = repmat(Ck, [1, 1, N + 1]);
Dk = zeros(ny, ny);

vv = param.vv;
vp = param.vp;
vCons = param.vCons;

DD = repmat(Dk, [1, 1, N + 1]);
scale = param.scale;
Signal = Signal/scale;

for i = 1 : N + 1
    
    for j = 1:SignalNum
        dist = sqrt((Xbar(1+4*(i-1) : 2+4*(i-1)) - Signal(j,:)')'*(Xbar(1+4*(i-1) : 2+4*(i-1)) - Signal(j,:)'));
        Dk(1+(j-1)*2 : 2+(j-1)*2, 1+(j-1)*2 : 2+(j-1)*2) = vp * dist * diag([0.8 1]) + vCons * diag([0.9 0.8]); 
    end 
    Dk(1+j*2 : 2+j*2, 1+j*2 : 2+j*2) = vv * diag([0.8 1]);
    DD(:, :, i) = Dk;
    
end

% Initial state covariance
P0 = PhatPrior0 + PtildePrior0;

% Initial estimation error covariance  PtildePrior0
% Initial estimated state covariance  PhatPrior0 

% Target maximum final covariance
Pf = PhatPriorf + PtildePriorf;

% Solve for sequence of kalmain gains LL, error covariances PPtilde, and prior
% error covariances PPtm
[LL, PPtilde, PPtm] = solveKF(AA, GG, CC, DD, N, PtildePrior0);

%% Time-varying LQR tracking

S(:,:,N+1) = Qk;
for k = N:-1:1
    
    K(:,:,k) = (BB(:,:,k)' * S(:,:,k+1) * BB(:,:,k) + Rk) \ BB(:,:,k)' * S(:,:,k+1) * AA(:,:,k);
    S(:,:,k) = Qk + AA(:,:,k)' * S(:,:,k+1) * AA(:,:,k) - AA(:,:,k)' * S(:,:,k+1) * BB(:,:,k) * K(:,:,k);

end

%% State covariance and estimated state covariance

PP(:,:,1)= blkdiag(P0, P0 - PPtilde(:,:,1));
TrPP = zeros(1, N+1);
TrPP(1) = trace(PP(1:4, 1:4, 1));
% TVLQG covariance control cost
Jc = 0;

for k = 1:N
    
    Fk = [AA(:,:,k)                               -BB(:,:,k) * K(:,:,k);
          LL(:,:,k+1) * CC(:,:,k+1) * AA(:,:,k)   AA(:,:,k) - BB(:,:,k) * K(:,:,k) - LL(:,:,k+1) * CC(:,:,k+1) *AA(:,:,k)];
    [row, ~] = size(GG(:,:,k));
    [~, col] = size(LL(:,:,k+1) * DD(:,:,k+1));
    Ek = [GG(:,:,k)                               zeros(row, col);
          LL(:,:,k+1) * CC(:,:,k+1) * GG(:,:,k)   LL(:,:,k+1) * DD(:,:,k+1)];
    PP(:,:,k+1) = Fk * PP(:,:,k) * Fk' + Ek * Ek';
    TrPP(k+1) = trace(PP(1:4, 1:4, k+1));
    Jc = Jc + trace( PP(5:8, 5:8, k) * (Qk + K(:,:,k)' * Rk * K(:,:,k)) );
    
end

CovCost = dt * Jc;

CollisonCost = 10000;

    xbar = zeros(nx, N + 1);
    for k = 0:N
        Ek = [zeros(nx, k * nx) eye(nx) zeros(nx,(N - k) * nx)];
        xbar(:, k + 1) = Ek * Xbar;
    end

    %% Monte Carlo
    % Number of monte carlo trials to display (can be zero)
    MCnum = 100;

    figure(10); 
%     figure(20);
    hold on
    CollisonCost = montecarlo(MCnum, AA, BB, GG, CC, DD, LL, xbar, ...
        N, PhatPrior0, PtildePrior0, V, K, param.scale);

    %% Plot result

    % Step numbers to draw covariance ellipses at
    ellipse_steps = 0:1:N;

    figure(10);
%     figure(1);
    drawResult(xbar(:,end), xbar, Pf, PP, PPtilde, ellipse_steps, param.scale)

end

