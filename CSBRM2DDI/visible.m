function visible_flag = visible(x, y, z, world)

camera = [2 5 3; 4.5 0 3; 8 5 3];
visible_flag = 0;
dim = 3;

for i = 1:size(camera,1)
    
    vec_norm = norm(camera(i,:) - [x, y, z]);
    pathx = linspace(x, camera(i,1), vec_norm/0.2);
    pathy = linspace(y, camera(i,2), vec_norm/0.2);
    pathz = linspace(z, camera(i,3), vec_norm/0.2);
    if ~MeanCollisionCheck([pathx; pathy; pathz], world, dim)
        visible_flag = 1;
        return;
    end
        
end

end