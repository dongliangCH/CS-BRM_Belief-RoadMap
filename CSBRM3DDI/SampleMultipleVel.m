function [points_new, Covs_new, DistM] = SampleMultipleVel(points, Covs, dim, l)

n = 1;
points_new(1, :) = points(1, :); % for i=1, n=1, the initial node is not changed
Covs_new(1, 1) = Covs(1, 1);
Covs_new(1, 2) = Covs(1, 2);
DistM = 1;

for i = 2:size(points,1)-l    
    vec = points(i,1:dim) - points(:,1:dim);
    dist = sqrt(sum(vec.*vec, 2));
    dist(i) = dist(i) + 10;
    neb = 1;
    near_idx = find(dist < neb);
    for j = 1:size(near_idx, 1)
%     if collisioncheck(points(i,1:2), points(near_idx(j),1:2),Cspace)            
        n = n + 1;
        vel = points(near_idx(j),1:dim)-points(i,1:dim);
        vel = vel / norm(vel);
        points_new(n, :) = [points(i,1:dim), vel]; 
        Covs_new(n, 1) = Covs(i, 1);
        Covs_new(n, 2) = Covs(i, 2);  
        DistM = [DistM neb];
%     end
    end
end

for i = (size(points,1)-l+1):size(points,1)
    neb = 1;
    n = n + 1;
    points_new(n, :) = points(i,:); 
    Covs_new(n, 1) = Covs(i, 1);
    Covs_new(n, 2) = Covs(i, 2);  
    DistM = [DistM neb];
end