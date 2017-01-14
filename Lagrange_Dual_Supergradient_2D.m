function [path, dist, ustar, vstar, LB, UB, its] =...
    Lagrange_Dual_Supergradient_2D(Df,Dg,Dh,Tg,Th, s, t)
% Lagrangian dual solver via supergradient methods
%   --- one Lagrange multiplier.
% INPUTs:
%   Df   ::IP distance matrix --- symmetric with non-negative entries
%   Dg  ::Additional constraint matrix --- symmetric with non-negative
%           entries
%   T   ::Maximum number of towns
%   s   ::Start node number
%   t   ::End node number
% OUTPUTs:
%   path::Dual-optimal path vector
%   dist::Distance along optimal path, accordind to distance matrix D
%   ustar::Optimal Lagrange multiplier
%   LB,UB:: Lower and upper bounds on optimal distance for IP

path = [];
Xbest = [];
dist = 0;
ustar = 0;
LB = 0;
UB = 10000;
Tol = 1; % Tolerance of 1 km for dual gap
Maxits = 200;

% --- Implicit functions:
    function Dval = Dist(X,Dmat)
    Dval = 0;
        for i = 1:length(X)-1
            j1 = X(i); j2 = X(i+1);
            Dval = Dval + Dmat(j1,j2);
        end
    end

% --- Initialise Lagrange multiplier
u=0; v = 0; its = 0;
% Commence iterations
while(its < Maxits)
    % Solve SP with constraint matrix Df + u*Dg + v*Dh
    Guv = graph(Df + u*Dg + v*Dh);
    Xuv = shortestpath(Guv, s, t,'Method','positive');
    Luv = Dist(Xuv, Df + u*Dg + v*Dh) - u*Tg - v*Th;
    % Check if improvement
    if Luv > LB % Better lower bound
        LB = Luv;
    end
    % Check RCSP-feasibility
    if Dist(Xuv, Dg) <= Tg && Dist(Xuv, Dh) <= Th
        % Xuv is RCSP feasible (but maybe not optimal). Tighten
        % the upper bound, and store best feasible solution
        cuv = Dist(Xuv, Df);
        if UB > cuv
            UB = cuv; Xbest = Xuv; ustar = u; vstar = v;
            %fprintf('\nSupergradient 2D: Feasible solution exists.\nNew UB = %.0f at iteration %i\n', UB, its)
        end
        % Now that bound is tightened, check for optimality
        if LB >= UB
            % Feasible solution is optimal
            %fprintf('\nSupergradient 2D: Optimal solution after %i its', its)
            path = Xuv; dist = cuv; ustar = u; vstar = v; LB = UB;
            return;
        end
    end
    % If this point is reached, we do not have an optimal feasible solution.
    % Supergradient procedure:
    Gamma_u = - Dist(Xuv, Dg);
    gamma_u = Dist(Xuv, Dg) - Tg;
    Gamma_v = - Dist(Xuv, Dh);
    gamma_v = Dist(Xuv, Dh) - Th;
    %7) Determine stepsize
    Tu = 2*(UB-LB)/(Gamma_u^2 + gamma_u^2 + Gamma_v^2 + gamma_v^2);
    Tv = 2*(UB-LB)/(Gamma_u^2 + gamma_u^2 + Gamma_v^2 + gamma_v^2);
    % Update Lagrange multipliers
    u = max(0, u + Tu*gamma_u);
    v = max(0, v + Tv*gamma_v);
    % Continue to next iteration
    its = its + 1;
end
% Maximum iterations reached an no optimal solution found. Return best
% iterate with bounds, or flag that no feasible solution was found.
%fprintf('\nSupergradient 2D: its = %i', its)
if isempty(Xbest)
    % No feasible point found --- return bounds
    fprintf('\nSupergradient 2D output not feasible. Dist bounds: [%.0f,%.0f]',...
        LB, UB)
    path = Xuv;
else
    %fprintf('\nSupergradient 2D output feasible. Dist bounds: [%.0f,%.0f]',...
    %    LB, UB)
    path = Xbest;
end
dist = Dist(path, Df);
end