classdef PotentialWaveSolver < BaseWaveSolver
    % PotentialWaveSolver 解决:
    % d²u/dt² = c² d²u/dx² - α u
    %
    % 与普通波方程相比，多了一个 - α u 项。
    
    properties
        alpha = 1.0 % 势项系数 α > 0
    end
    
    methods
        function obj = PotentialWaveSolver(length, c, alpha, nx, dt, t_steps)
            if nargin >= 1; obj.length = length; end
            if nargin >= 2; obj.c = c; end
            if nargin >= 3; obj.alpha = alpha; end
            if nargin >= 4; obj.nx = nx; end
            if nargin >= 5; obj.dt = dt; end
            if nargin >= 6; obj.t_steps = t_steps; end
            
            obj.dx = obj.length / obj.nx;
            obj.x = linspace(0, obj.length, obj.nx);
            obj.u = zeros(1, obj.nx);
            obj.u_prev = zeros(1, obj.nx);
            obj.u_next = zeros(1, obj.nx);
            
            obj.checkStability();
        end
        
        function checkStability(obj)
            % 基本波方程的稳定性条件仍是参考 c*dt/dx <= 1
            % 有势项时，这个条件依然是一个基本要求，但可能更严格。
            if obj.c * obj.dt / obj.dx > 1
                error('Stability condition violated: c * dt / dx <= 1 for potential wave equation.');
            end
        end
        
        function step(obj)
            cdt_dx = (obj.c * obj.dt / obj.dx)^2;
            alpha_dt2 = obj.alpha * (obj.dt)^2;
            for i = 2:obj.nx-1
                obj.u_next(i) = 2*obj.u(i) - obj.u_prev(i) ...
                    + cdt_dx*(obj.u(i+1)-2*obj.u(i)+obj.u(i-1)) ...
                    - alpha_dt2*obj.u(i);
            end
            obj.applyBoundaryConditions();
            [obj.u_prev, obj.u] = deal(obj.u, obj.u_next);
        end
        
        function applyBoundaryConditions(obj)
            % 固定边界条件
            obj.u_next(1) = 0;
            obj.u_next(end) = 0;
        end
    end
end
