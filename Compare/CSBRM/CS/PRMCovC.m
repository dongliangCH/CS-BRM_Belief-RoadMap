function [CovCost, CollisonCost, K, problem] = PRMCovC(PhatPrior0, PtildePrior0, PhatPriorf, PtildePriorf, N, Signal, param, Xbar, V, Cspace)

% Notation: Ak is single step of matrix, AA is 3d array with each Ak
% stacked in the 3rd dimension, A is block matrix

%% Problem setup

% Number of steps N

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

L = makeInputMatrixType2(AA, LL, N);

%% Covariance steering
% Innovation process covariance
PYtilde = zeros((N + 1) * ny, (N + 1) * ny);
for i = 0:N
    ndx = (i * ny) + (1:ny);
    PYtilde(ndx, ndx) = CC(:, :, i + 1) * PPtm(:, :, i + 1) ...
        * CC(:, :, i + 1)' + DD(:, :, i + 1) * DD(:, :, i + 1)';
end

% Feedback process covariance
PZ = A * PhatPrior0 * A' + L * PYtilde * L';

% Flag if should use history of disturbances in feedback
useHistory = false;

% Solve covariance steering with yalmip
% K = Kcov;
% Jc = 2.7734;
[K, CovCost, cpu_time2, problem] = solveCS(B, Q, R, N, nx, nu, PZ, PhatPrior0, PhatPriorf, useHistory);

CovCost = dt * CovCost;

CollisonCost = 10000;

if ~problem

    % Filtered state covariance (block matrix)
    IplusBK = eye((N + 1) * nx) + B * K;
    Phat = IplusBK * PZ * IplusBK';

    % Mean and covariance as normal sequences
    PPhat = zeros(nx, nx, N + 1);
    PP = zeros(nx, nx, N + 1);
    xbar = zeros(nx, N + 1);
    for k = 0:N
        Ek = [zeros(nx, k * nx) eye(nx) zeros(nx,(N - k) * nx)];
        PPhat(:, :, k + 1) = Ek * Phat * Ek';
        PP(:, :, k + 1) = PPhat(:, :, k + 1) + PPtilde(:, :, k + 1);
        xbar(:, k + 1) = Ek * Xbar;
    end

    %% Monte Carlo
    % Number of monte carlo trials to display (can be zero)
    MCnum = 100;

%     figure(20); 
    figure(10); 
    hold on
    CollisonCost = montecarlo(MCnum, AA, BB, GG, CC, DD, LL, xbar(:,1), ...
        N, PhatPrior0, PtildePrior0, V, K, useHistory);


    %% Plot result
    % Step numbers to draw covariance ellipses at
    ellipse_steps = 0:1:N;
    
%     figure(1);
    figure(10); 
    drawResult(MCnum, AA, BB, GG, CC, DD, LL, xbar(:,1), xbar(:,end), xbar, ...
        Pf, N, PhatPrior0, PtildePrior0, V, K, PP, PPhat, ...
        PPtilde, ellipse_steps, useHistory)

end

end

