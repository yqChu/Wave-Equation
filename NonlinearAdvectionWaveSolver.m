classdef NonlinearAdvectionWaveSolver < BaseWaveSolver
    % NonlinearAdvectionWaveSolver 实现具有非线性项 u * ∂u/∂x 的波动方程:
    %
    % d²u/dt² = c² d²u/dx² + α u (∂u/∂x)
    %
    % 与线性方程相比，这里有非线性对流项，使数值稳定性和精度更具挑战性。
    
    properties
        alpha = 1.0 % 非线性项系数 α
    end
    
    methods
        function obj = NonlinearAdvectionWaveSolver(length, c, alpha, nx, dt, t_steps)
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
            % 基本的线性CFL条件为 c*dt/dx <= 1。
            % 非线性项存在时需更加严格选择dt和dx。
            if obj.c * obj.dt / obj.dx > 1
                error('Stability condition violated: c * dt / dx <= 1 (nonlinear advection may require even stricter conditions).');
            end
        end
        
        function step(obj)
            cdt_dx = (obj.c * obj.dt / obj.dx)^2;
            half = 1/(2*obj.dx);
            alpha_dt2 = obj.alpha * (obj.dt)^2;
            
            for i = 2:obj.nx-1
                u_mid_x = (obj.u(i+1) - obj.u(i-1)) * half; % ∂u/∂x 的近似
                obj.u_next(i) = 2*obj.u(i) - obj.u_prev(i) ...
                    + cdt_dx*(obj.u(i+1)-2*obj.u(i)+obj.u(i-1)) ...
                    + alpha_dt2 * obj.u(i) * u_mid_x;
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
