

function drawResult(xbar, PP, PPtilde, ellipse_steps)

rl = 0.9545; % 2 Sigma   
hold on;

for k = ellipse_steps

%     % Total covariance
%     pl1 = error_ellipse('C', PP(1:2, 1:2, k + 1), 'mu', xbar(1:2, k + 1), ...
%         'conf', rl, 'style', 'k', 'linewidth', 0.567);
    
%     % Estimated state covariance
%     pl2 = error_ellipse('C', PP(5:6, 5:6, k + 1), 'mu', xbar(1:2, k + 1), ...
%         'conf', rl, 'style', 'k', 'linewidth', 0.567);
       
    % State error covariance
    pl3 = error_ellipse('C', PPtilde(1:3, 1:3, k + 1), 'mu', ...
        xbar(1:3, k + 1), ...
        'conf', rl, 'style', 'r', 'linewidth', 0.567);
%     set(pl3, 'Color', [1 0 0])
    
end

% lh = legend([pl1, pl3], '$P_k$', '$\tilde{P}_k$', 'location', 'southwest');
% set(lh, 'Interpreter', 'latex');
% hold off;
% grid on;
% axis equal;

% x1 = xlabel('$x_1$', 'interpreter', 'latex');
% y1 = ylabel('$x_2$', 'interpreter', 'latex');
% set(x1,'FontSize',18);
% set(y1,'FontSize',18);
% set(gca,'FontSize',16,'FontName','Times');

end