function [newSample, velRand, Covs] = sampleNodes(Cspace, k, Sample, Signal, display, param)

%% Sample points
rng(0);
% i = 0;
% newSample = zeros(k, 2);
% % Iteratively add vertices
% while i < k 
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

velRand = zeros(k+2, 2);
% for i = 1 : k+2
%     velRand(i,:) = [rand * 2 - 1, rand * 2 - 1];
% end

%% Kalman filter
% Dynamics
N = 30;
nx = 4;
nu = 2;
nw = 4;
dt = param.dt;
Ak = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
Bk = [dt ^ 2 / 2 0; 0 dt ^ 2 / 2; dt 0; 0 dt];
Gk = 0.10 * eye(nx);
% Sequences of system matricies
AA = repmat(Ak, [1, 1, N]);
GG = repmat(Gk, [1, 1, N]);

[SignalNum, ~] = size(Signal); 
ny = (nx-2) * SignalNum + 2;
Ck = repmat([1 0 0 0; 0 1 0 0], [SignalNum, 1]);
Ck = [Ck; 0 0 1 0; 0 0 0 1]; 
CC = repmat(Ck, [1, 1, N + 1]);
Dk = zeros(ny, ny);

% Cost
Qk = 2 * blkdiag(2, 2, 2, 2);
Rk = 1 * blkdiag(2, 2);

vv = param.vv;
vp = param.vp;
vCons = param.vCons;
scale = param.scale;
Signal = Signal/scale;
newSample = newSample/scale;

for i = 1:k+2

    for j = 1:SignalNum
        dist = sqrt((newSample(i,:) - Signal(j,:))*(newSample(i,:) - Signal(j,:))');
        Dk(1+(j-1)*2 : 2+(j-1)*2, 1+(j-1)*2 : 2+(j-1)*2) = vp * dist * eye(2) + vCons * eye(2); 
    end 
    Dk(1+j*2 : 2+j*2, 1+j*2 : 2+j*2) = vv * eye(2);
    DD = repmat(Dk, [1, 1, N + 1]);
    
    %% SLQG

    As = Ak; Bs = Bk; Gs = Gk; Cs = Ck; Ds = Dk; Qs = Qk; Rs = Rk;
    
    % Initial estimation error covariance
    PtildePrior0 = 2 * blkdiag(0.1, 0.1, 0.05, 0.05).* diag([0.8 0.8 0.7 0.7]);

    % Solve for sequence of kalmain gains LL, error covariances PPtilde, and prior
    % error covariances PPtm
    [LLs, PPtildes, PPtms] = solveKF(AA, GG, CC, DD, N, PtildePrior0);
    
    % SLQR
    [~, ~, Ks] = dare(As, Bs, Qs, Rs);

    %% State covariance and estimated state covariance

    PPs(:,:,1)= diag([0.2 0.2 0.2 0.2 0.1 0.1 0.1 0.1]);
    TrPPs = zeros(1, N + 1);
    TrPPs(1) = trace(PPs(1:4, 1:4, 1));
    for k = 1:N
        Fk = [As                         -Bs * Ks;
            LLs(:,:,k+1) * Cs * As     As - Bs * Ks - LLs(:,:,k+1) * Cs *As];      
        [row, ~] = size(Gs);
        [~, col] = size(LLs(:,:,k+1) * Ds);
        Ek = [Gs                         zeros(row, col);
            LLs(:,:,k+1) * Cs * Gs     LLs(:,:,k+1) * Ds];
        PPs(:,:,k+1) = Fk * PPs(:,:,k) * Fk' + Ek * Ek';
        TrPPs(k+1) = trace(PPs(1:4, 1:4, k+1));
    end

%     % Trace of the state covariance
%     figure(7)
%     plot(TrPPs, '.-');
%     x1 = xlabel('$Time \ steps$', 'interpreter', 'latex');
%     y1 = ylabel('$tr(P)$', 'interpreter', 'latex');
%     set(x1,'FontSize',18);
%     set(y1,'FontSize',18);
%     set(gca,'FontSize',16,'FontName','Times');

    PPtmNs = 1.2 * PPtms(:,:,end);
    
    % Initial estimated state covariance
    PhatPrior0 = 1.4 * PPs(5:8,5:8,end) + 1.0 * diag([0.05 0.05 0.04 0.04]);
    
    % Initial state covariance
    P0 = PhatPrior0 + PPtmNs;
    
    Covs{i,1} = PhatPrior0;
    Covs{i,2} = PPtmNs;

    rl = 0.9973; % 3 Sigma   

    figure(1); hold on
    % Target covariance
    plT = error_ellipse('C', P0(1:2, 1:2).*400, 'mu', [newSample(i, 2); newSample(i, 1)].*20, ...
        'conf', rl, 'style', 'k', 'linewidth', 2);

end

newSample = scale * newSample;
    
end

