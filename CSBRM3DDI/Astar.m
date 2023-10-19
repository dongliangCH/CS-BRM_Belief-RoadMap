function [count_vertex, Path, cost] = Astar(Vertices, Edges, EdgesCost, start_idx, goal_idx, weight)
N = size(Vertices,1);
costToCome = ones(N,1) * 1e6;
costToCome(start_idx) = 0; 

% h_cost = Vertices(:, 1:3) - Vertices(goal_idx, 1:3);
% h_cost = weight * sqrt(sum(h_cost.*h_cost,2));
h_cost = zeros(N,1);  % Dijkstra

QueueVertex =  [start_idx 0+h_cost(start_idx)];
back_pointer = zeros(N,1);

count_vertex = 0;
while ~isempty(QueueVertex)
    [~, pop_idx] = min(QueueVertex(:, 2));
    vertex_idx = QueueVertex(pop_idx, 1);
    if vertex_idx == goal_idx
        cost = costToCome(goal_idx);
        Path = PlotResult(back_pointer, start_idx, goal_idx);
        return
    end    
    count_vertex = count_vertex + 1;
    for i = 1:Vertices(vertex_idx, 7)
        child_idx = Edges(vertex_idx, i);
        new_cost = costToCome(vertex_idx) + EdgesCost(vertex_idx, i);
        if new_cost < costToCome(child_idx)
            costToCome(child_idx) = new_cost;
            back_pointer(child_idx) = vertex_idx;
            idx = QueueVertex(:, 1) == child_idx;
            if sum(idx)
                QueueVertex(idx, 2) = new_cost+h_cost(child_idx);
            else
                QueueVertex(end+1, :) = [child_idx new_cost+h_cost(child_idx)];
            end
        end
    end
    QueueVertex(pop_idx, :) = [];
end
end

function Path = PlotResult(back_pointer, start_idx, goal_idx)
Path = [goal_idx];
current_idx = goal_idx;
while current_idx ~= start_idx
    current_idx = back_pointer(current_idx);
    Path = [current_idx Path];
end
end