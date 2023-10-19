function plotCovariance(xbarf, Pf)
    
    rl = 0.9545; % 2 Sigma   
    error_ellipse('C', Pf(1:3, 1:3), 'mu', xbarf(1:3), ...
        'conf', rl, 'style', 'k', 'linewidth', 0.567);
    plot3(xbarf(1), xbarf(2), xbarf(3),'Marker','.','MarkerSize',18,'MarkerEdgeColor','[0.8500 0.3250 0.0980]')
    
end