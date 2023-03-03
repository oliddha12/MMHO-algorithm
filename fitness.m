function fx = fitness(x,G)
    S = [1 1];       % start point
    E = size(G);     % Endpoint
    route = [S(1) x E(1)];
    dim = length(route);
    nB = 0;        
    route=round(route);
    for j = 2 : dim-1
       if G(route(j),j) == 1
           nB = nB + 1;
       end
    end
    if nB == 0      
        path=GenerateRoute(route,G); 
        path=GenerateSmoothPath(path,G);
        path=GenerateSmoothPath(path,G);
        fx = 0;
        for i = 1:size(path,1)-1
            fx = fx + sqrt((path(i+1,1)-path(i,1))^2 + (path(i+1,2)-path(i,2))^2);
        end
    else
    fx = E(1)*E(2) * nB;
    end