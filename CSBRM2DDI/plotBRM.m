plotWorld(world, dim); 
for i = 1:size(EdgeTraj,1)
    for j = 1:size(EdgeTraj,2)
        meanTraj = EdgeTraj{i,j};
        if ~isempty(meanTraj)
            plot3(meanTraj(1,:), meanTraj(2,:), meanTraj(3,:), 'color', [0 0.4470 0.7410], 'LineWidth', 1);
            plot3(meanTraj(1,end), meanTraj(2,end), meanTraj(3,end), 'Marker','.','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]');
            plot3(meanTraj(1,1), meanTraj(2,1), meanTraj(3,1), 'Marker','.','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]');
        end
    end
end