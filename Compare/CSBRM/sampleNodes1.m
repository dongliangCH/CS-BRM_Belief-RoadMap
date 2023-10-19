function [samples, velRand, Covs] = sampleNodes1(Cspace, k, Sample, Signal, display, param)

%% Sample newSample
rng(0);
% i = 0;
% newSample = zeros(k, 2);
% while i < k  % iteratively add vertices
%     x = double(int32(rand(1, 2) .* size(Cspace)));
%     if checkpoint(x, Cspace)        
%         newSample(i+1, :) = x;
%         i = i + 1;
%         if display
%             rectangle('Position',[x(2)-5,x(1)-5,10,10],'Curvature',[1,1],'FaceColor','[0.5 0.5 0.5]'); 
%             hold on 
%         end
%     end
% end
newSample = [30, 320; 50, 50; 322, 31; 205, 170; 80, 220; 240, 270; 301, 369; 164, 458;  ... 
    41, 440; 153, 335; 347, 289; 335, 161; 459, 179; 460, 33; 454, 310; 458, 445];
for i = 1:size(newSample,1)
    rectangle('Position',[newSample(i,2)-5,newSample(i,1)-5,10,10],'Curvature',[1,1],'FaceColor','[0.5 0.5 0.5]'); hold on
end

newSample = [Sample; newSample];

% Sample velocity
n = 0;
for i = 1:size(newSample,1)    
    vec = newSample(i,1:2) - newSample(:,1:2);
    dist = sqrt(sum(vec.*vec, 2));
    dist(i) = dist(i) + size(Cspace,1) * 2;
    near_idx = find(dist < param.neb);
    for j = 1:size(near_idx, 1)
        if collisioncheck(newSample(i,1:2), newSample(near_idx(j),1:2),Cspace)
            
            n = n + 1;
            vel = newSample(near_idx(j),1:2)-newSample(i,1:2);
            vel = sqrt(4) * vel / norm(vel);    % rand * sqrt(2) * vel / norm(vel);
            Nodes(n, :) = [newSample(i,1:2), vel];                       
           
        end
    end
end


%% Kalman filter
% Dynamics
N = 20;
nx = 4;
nu = 2;
nw = 4;
dt = param.dt;
Ak = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
Gk = 0.1 * eye(nx);
% Sequences of system matricies
AA = repmat(Ak, [1, 1, N]);
GG = repmat(Gk, [1, 1, N]);

ny = 4;
Ck = eye(4);
CC = repmat(Ck, [1, 1, N + 1]);
power = 1/2;
scale = param.scale;
Signal = Signal/scale;
Nodes(:, 1:2) = Nodes(:, 1:2)/scale;

for i = 1:size(Nodes, 1)

    dd = 1 / (((Nodes(i,1:2) - Signal(1,:))*(Nodes(i,1:2) - Signal(1,:))'/2 + 0.000001)^(-power) + ...
              ((Nodes(i,1:2) - Signal(2,:))*(Nodes(i,1:2) - Signal(2,:))'/2 + 0.000001)^(-power) + ...
              ((Nodes(i,1:2) - Signal(3,:))*(Nodes(i,1:2) - Signal(3,:))'/2 + 0.000001)^(-power) + ...
              ((Nodes(i,1:2) - Signal(4,:))*(Nodes(i,1:2) - Signal(4,:))'/2 + 0.000001)^(-power) + ...
              ((Nodes(i,1:2) - Signal(5,:))*(Nodes(i,1:2) - Signal(5,:))'/2 + 0.000001)^(-power) + ...
              ((Nodes(i,1:2) - Signal(6,:))*(Nodes(i,1:2) - Signal(6,:))'/2 + 0.000001)^(-power)); 
    Dk = param.Dkparam * blkdiag(0.05, 0.05, 0.05, 0.05);
    Dk = dd^2 * Dk;
    DD = repmat(Dk, [1, 1, N + 1]);
    

% Initial estimation error covariance
PtildePrior0 = 2 * blkdiag(0.1, 0.1, 0.05, 0.05).* diag([0.8 0.8 0.7 0.7]);

% Solve for sequence of kalmain gains LL, error covariances PPtilde, and prior
% error covariances PPtm
[LL, PPtilde, PPtm] = solveKF(AA, GG, CC, DD, N, PtildePrior0);

PtildeN = PPtilde(:,:,end);
PPtmN = PPtm(:,:,end);

% Initial estimated state covariance
PhatPrior0 = (0 * dd + 3 * 0.5) * 2 * blkdiag(0.05, 0.05, 0.05, 0.05);
% PhatPrior0 = 6 * PPtmN + blkdiag(0.05, 0.05, 0.05, 0.05); 

% Initial state covariance
P0 = PhatPrior0 + 1.3 * PPtmN;

Covs{i,1} = PhatPrior0;
Covs{i,2} = 1.3 * PPtmN;

rl = 0.9973; % 3 Sigma   
hold on;
% Target covariance
plT = error_ellipse('C', P0(1:2, 1:2).*400, 'mu', [Nodes(i, 2); Nodes(i, 1)].*20, ...
    'conf', rl, 'style', 'k', 'linewidth', 2);

end

samples = Nodes(:, 1:2) * scale;
velRand = Nodes(:, 3:4);

end

