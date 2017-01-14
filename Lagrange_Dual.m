function [path, dist, ustar, LB, UB, its] =...
    Lagrange_Dual_Cutting_Plane(D,Dg,T, s, t)
% Lagrangian dual solver via cutting-plane method
%   --- one Lagrange multiplier.
% INPUTs:
%   D   ::IP distance matrix --- symmetric with non-negative entries
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

% --- Initialise:
    path = [];
    dist = 0;
    ustar = 0;
    its = 0;
    maxits = 100;
    LB = 0;
    UB = Inf;
    Tol = 1; % Tolerance of 1 km for LD
    
    % --- Implicit functions:
    function fval = f(X)
    fval = 0;
        for i = 1:length(X)-1
            j1 = X(i); j2 = X(i+1);
            fval = fval + D(j1,j2);
        end
    end

    function gval = g(X)
    gval = - T;
        for i = 1:length(X)-1
            j1 = X(i); j2 = X(i+1);
            gval = gval + Dg(j1,j2);
        end
    end

% --- Algorithm initialisation stage
% Create graphs
Gf = graph(D);
Gg = graph(Dg);
% Get shortest paths using Dijkstra's algorithm
[Xp,cp] = shortestpath(Gf, s, t,'Method','positive');
gp = g(Xp);
    % Check if SP solves RCSP
    if gp <= 0
        path = Xp;
        dist = cp;
        LB = cp;
        UB = cp;
        return;
    end
    
% SP does not solve RCSP. Check feasiblity of RCSP
Xm = shortestpath(Gg, s, t,'Method','positive');
cm = f(Xm); gm = g(Xm);
    if gm > 0 % RCSP is infeasible
       fprintf('\n\nProblem is intractible with T = %i', T);
       return;
    end
    % RCSP is feasible. Move onto Lagrangian dual stage (Algorithm 1)

% Algorithm 1 --- Cutting-plane for Lagrangian dual
exit_flag = 0;
UB = cm; % Initialise upper bound
its = 0;
    while(~exit_flag && its <= maxits)
        its = its + 1;
        u = (cm-cp)/(gp-gm);
        ru = cp + u*gp;
        Gu = graph(D + u*Dg); % Lagrangian dual-weighted graph
        % Solve SP with graph Gu
        Xu = shortestpath(Gu, s, t,'Method','positive');
        gu = g(Xu); cu = f(Xu); Lu = cu + u*gu;       
            if abs(gu) < Tol % i.e. g - lambda == 0
                % Xu solves contsrained IP optimally. LB = UB.
                path = Xu;
                dist = cu;
                UB = cu;
                LB = cu;
                ustar = u;
                return;
            elseif abs(Lu - ru) < Tol && gu < 0
                % Xu is feasible for RCSP, and algorithm gives no more
                % improvement in LB. Update LB and UB, then stop.
                LB = Lu;
                UB = min(UB, cu);
                path = Xu;
                exit_flag = 1;
            elseif abs(ru - Lu) < Tol && gu > 0
                % Xu infeasible for RCSP, and algorithm gives no more
                % improvement in LB.
                LB = ru;
                UB = cm;
                exit_flag = 1;
            elseif ru < Lu && gu > 0
                % X infeasible for RCSP.
                cp = cu;
                gp = gu;
            else
                % Xu is feasible for RCSP, but algorithm still giving
                % improvements on LB. Update bounds and continue.
                %Xm = X;
                cm = cu;
                gm = gu;
                UB = min(UB,cu);
                LB = max(LB, Lu);
            end
    end
   
end
