

function drawResult(MCnum, AA, BB, GG, CC, DD, LL, xbar0, xbarf, xbar, ...
    Pf, N, PhatPrior0, PtildePrior0, V, K, PP, PPhat, PPtilde, ...
    ellipse_steps, useHistory)

rl = 0.9545; % 2 Sigma 
hold on;

% Target covariance
plT = error_ellipse('C',Pf(1:3, 1:3),'mu', [xbarf(1); xbarf(2); xbarf(3)], ...
    'conf', rl, 'style', 'k', 'linewidth', 2);

for k = ellipse_steps

    % Total covariance
    pl1 = error_ellipse('C', PP(1:3, 1:3, k + 1), 'mu', xbar(1:3, k + 1), ...
        'conf', rl, 'style', 'k', 'linewidth', 0.567);
    
%     % Filtered state covariance
%     pl2 = error_ellipse('C', PP(5:6, 5:6, k + 1), 'mu', xbar(1:2, k + 1), ...
%         'conf', rl, 'style', 'k', 'linewidth', 0.567);
       
    % State error covariance
    pl3 = error_ellipse('C', PPtilde(1:3, 1:3, k + 1), 'mu', ...
        xbar(1:3, k + 1), ...
        'conf', rl, 'style', 'r', 'linewidth', 0.567);
    
end

end