classdef NonlinearPotentialWaveSolver < BaseWaveSolver
    % NonlinearWaveSolver 实现非线性波动方程:
    % d²u/dt² = c² d²u/dx² - α (u^2)
    %
    % 该方程是非线性的，不再满足叠加原理。数值求解时需谨慎选择时间步长和空间步长。
    
    properties
        alpha = 1.0 % 非线性项系数 α > 0
    end
    
    methods
        function obj = NonlinearPotentialWaveSolver(length, c, alpha, nx, dt, t_steps)
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
            % 对于非线性方程，简单的线性CFL条件不再严格适用。
            % 这里仍要求 c*dt/dx <= 1 作为基本稳定性参考。
            % 实际中可能需要更小的 dt 来确保非线性稳定。
            if obj.c * obj.dt / obj.dx > 1
                error('Stability condition violated: c * dt / dx <= 1 (nonlinear case may require stricter conditions).');
            end
        end
        
        function step(obj)
            cdt2_dx2 = (obj.c * obj.dt / obj.dx)^2;
            alpha_dt2 = obj.alpha * (obj.dt)^2;
            for i = 2:obj.nx-1
                obj.u_next(i) = 2*obj.u(i) - obj.u_prev(i) ...
                    + cdt2_dx2*(obj.u(i+1)-2*obj.u(i)+obj.u(i-1)) ...
                    - alpha_dt2*(obj.u(i))^2;
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
