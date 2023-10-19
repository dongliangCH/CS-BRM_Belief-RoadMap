function collision_flag = MeanCollisionCheck(MeanTraj, world, dim)

collision_flag = 0;

if collision_flag == 0 && dim == 2
    
    for j = 2 : size(MeanTraj, 2)-1
    p = MeanTraj(:, j);
    
    for i=1:dim
        if (p(i)>world.endcorner(i))||(p(i)<world.origincorner(i))
            collision_flag = 1;
            return;    
        end
    end
    
      % check each obstacle
      for i=1:world.NumObstacles
        if p(1) >= world.ox(i) && p(1) <= world.ox(i)+world.oa(i) && p(2) >= world.oy(i) && p(2) <= world.oy(i)+world.ob(i)
            collision_flag = 1;
            return;    
        end
      end
    end

%%%%% dim=3 case %%%%%
elseif collision_flag == 0 && dim ==3
    
    for j = 2 : size(MeanTraj, 2)-1
    p = MeanTraj(:, j);
    
    for i=1:dim
        if (p(i)>world.endcorner(i))||(p(i)<world.origincorner(i))
            collision_flag = 1;
            return;   
        end
    end
    
      % check each obstacle
      for i=1:world.NumObstacles
        if p(1) >= world.ox(i) && p(1) <= world.ox(i)+world.oa(i) && p(2) >= world.oy(i) && p(2) <= world.oy(i)+world.ob(i) && p(3) >= world.oz(i) && p(3) <= world.oz(i)+world.oc(i)
            collision_flag = 1;
            return;    
        end
      end
    end
    
end
end