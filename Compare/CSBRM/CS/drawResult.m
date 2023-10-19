

function drawResult(MCnum, AA, BB, GG, CC, DD, LL, xbar0, xbarf, xbar, ...
    Pf, N, PhatPrior0, PtildePrior0, V, K, PP, PPhat, PPtilde, ...
    ellipse_steps, useHistory)

rl = 0.9973; % 3 Sigma   
hold on;

% Target covariance
plT = error_ellipse('C',Pf(1:2, 1:2).*400,'mu', [xbarf(2); xbarf(1)].*20, ...
    'conf', rl, 'style', 'k', 'linewidth', 2);

for k = ellipse_steps

    % Total covariance
    pl1 = error_ellipse('C', PP(1:2, 1:2, k + 1).*400, 'mu', [xbar(2, k + 1); xbar(1, k + 1)].*20, ...
        'conf', rl, 'style', 'k', 'linewidth', 0.567);
    
%     % Filtered state covariance
%     pl2 = error_ellipse('C', PPhat(1:2, 1:2, k + 1), 'mu', xbar(1:2, k + 1), ...
%         'conf', rl, 'style', 'k', 'linewidth', 0.567);
%     set(pl2, 'Color', 0.6 *  ones(3, 1))
%     
%     % State error covariance
%     pl3 = error_ellipse('C', PPtilde(1:2, 1:2, k + 1), 'mu', ...
%         xbar(1:2, k + 1), ...
%         'conf', rl, 'style', 'k:', 'linewidth', 0.567);
%     set(pl3, 'Color', 0.2 *  ones(3, 1))
    
end

% lh = legend([pl1, p12, p13, plT], '$P_k$', '$\hat{P}_k$', ...
%     '$\tilde{P}_k$', '$P_f$', 'location', 'southwest');
% lh = legend([pl1, plT], '$P_k$', '$P_f$', 'location', 'southwest');
% set(lh, 'Interpreter', 'latex');
% hold off;
% grid on;
% axis equal;
x1 = xlabel('$x_1$', 'interpreter', 'latex');
y1 = ylabel('$x_2$', 'interpreter', 'latex');
set(x1,'FontSize',18);
set(y1,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
xlim([-10 510])
ylim([-10 510])

end

function b = vec(a)

b = a(:);

end