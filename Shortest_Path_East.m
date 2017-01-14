% ------------
% Shortest path --- Australian East Coast
clear all;
% Load data (Aus_Dist_Mat, Des, N_towns,...):
load Aus_Coords_East

% --- Set parameters
Tol = 1; % Tolerance of 1 km for shortest path bounds
Dmax = 110; % Maximum running distance on graph
Tg = 44; % Maximum total number of towns allowed to visit
% Tg = 47; For 2D supergradient
Th = 1; % Maximum number of desert towns
% Singly-constrained RCSP feasible with Tg >= 40
% Multiply-constrained RCSP feasble with Tg >= 43, Th = 1

% --- Setup graph
% Adjacency matrix D (distance matrix):
D = Aus_Dist_Mat_East;
% Number of nodes
N = size(D,1);
% Ensure distance matrix is symmetric
D = triu(D,1);
D = D + D';

% --- Additional constraints
% i) Upper limit on number of towns
Dg = ones(N);
% ii) Penalty for stopping in desert towns
Dh = D.*Des; % Des = indicator matrix for towns in desert

% --- Construct graph connectivity:
% Remove edges that are not walkable, and sparsify distance matrices
I = (D>Dmax);
D(I) = 0;
D = sparse(D);
Dg = triu(Dg,1);
Dg = Dg + Dg';
Dg(I) = 0;
Dg = sparse(Dg);
Dh = tril(Dh,1);
Dh = Dh + Dh';
Dh(I) = 0;
Dh = sparse(Dh);

% --- Create graphs
G = graph(D);
Gg = graph(Dg);
Gh = graph(Dh);

% --- Solve unconstrained SPs:
[path,d] = shortestpath(G, N_Start, N_End,'Method','positive');
%[path_g] = shortestpath(Gg, N_Start, N_End,'Method','positive');
%[path_h] = shortestpath(Gh, N_Start, N_End,'Method','positive');
% Print output information
fprintf('\n\nUnconstrained SP:')
fprintf('\nPath length: %i', length(path))
fprintf('\nDist: %.0f', d)
fprintf('\n')





% --- Lagrangian dual --- Cutting-plane algorithm:
[path_LD, dist_LD, u, LB, UB, its_cp] = ...
    Lagrange_Dual_Cutting_Plane(D, Dg, Tg, N_Start, N_End);
% Print output information
fprintf('\n\nCutting-Plane:')
fprintf('\nPath length: %i', length(path_LD))
if UB - LB <= 1
    fprintf('\nPath is RCSP optimal')
else
    fprintf('\nPath is not RCSP optimal')
    fprintf('\nBounds: [%.0f,%.0f]', LB, UB)
    fprintf('\nDuality gap = %.0f km', UB - LB)
end
fprintf('\nDist: %.0f', dist_LD)
fprintf('\nIterations: %i', its_cp)
fprintf('\n')




% --- Lagrangian dual --- Supergradient algorithm:
[path_LD_sup, dist_LD_sup, u_sup, LB_sup, UB_sup, its_sup] = ...
    Lagrange_Dual_Supergradient(D, Dg, Tg, N_Start, N_End);
% Print output information
fprintf('\n\nSupergradient:')
fprintf('\nPath length: %i', length(path_LD_sup))
if UB_sup - LB_sup <= 1
    fprintf('\nPath is RCSP optimal')
else
    fprintf('\nPath is not RCSP optimal')
    fprintf('\nBounds: [%.0f,%.0f]', LB_sup, UB_sup)
    fprintf('\nDuality gap = %.0f km', UB_sup - LB_sup)
end
fprintf('\nDist: %.0f', dist_LD_sup)
fprintf('\nIterations: %i', its_sup)
fprintf('\n')





% --- Multiply constrained Lagrangian dual --- Supergradient algorithm:
[path_LD_2D_sup, dist_LD_2D_sup, u_2D_sup, v_2D_sup, LB_2D_sup,...
    UB_2D_sup, its_2D_sup] = ...
    Lagrange_Dual_Supergradient_2D(D, Dg, Dh, Tg, Th, N_Start, N_End);
% Print output information
fprintf('\n\nSupergradient 2D:')
fprintf('\nPath length: %i', length(path_LD_2D_sup))
if UB_2D_sup - LB_2D_sup <= 1
    fprintf('\nPath is RCSP optimal')
else
    fprintf('\nPath is not RCSP optimal')
    fprintf('\nBounds: [%.0f,%.0f]', LB_2D_sup, UB_2D_sup)
    fprintf('\nDuality gap = %.0f km', UB_2D_sup - LB_2D_sup)
end
fprintf('\nDist: %.0f', dist_LD_2D_sup)
fprintf('\nIterations: %i', its_2D_sup)
fprintf('\n')


% --- Solution visialisation
figure()
landareas = shaperead('landareas.shp','UseGeoCoords',true);
%axesm ('pcarree', 'Frame', 'on', 'Grid', 'on');
geoshow(landareas,'FaceColor',[0.5 1.0 0.5],'EdgeColor',[.6 .6 .6]);
geoshow(Lats,Longs, 'DisplayType', 'point')
geoshow(Lats(path),Longs(path),'LineWidth',3)
%geoshow(Lats(path_g),Longs(path_g),'LineWidth',1,'Color','g')
%geoshow(Lats(path_h),Longs(path_h),'LineWidth',1,'Color','b')
%geoshow(Lats(path_LD),Longs(path_LD),'LineWidth',3,'Color','k')
%geoshow(Lats(path_LD_sup),Longs(path_LD_sup),'LineWidth',2,'Color','y')
geoshow(Lats(path_LD_2D_sup),Longs(path_LD_2D_sup),'LineWidth',3,'Color','k')
%legend('Towns', 'Unconstrained SP', '1st Constraint', '2nd Constraint',...
%    'LD - Cutting plane', 'LD - Subgradient','LD - 2D Subgradient',...
%    'Location','west')
title('Australian East-Coast towns')
xlabel('Longitude')
ylabel('Latitude')