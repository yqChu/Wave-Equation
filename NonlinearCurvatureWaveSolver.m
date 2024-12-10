classdef NonlinearCurvatureWaveSolver < BaseWaveSolver
    % NonlinearCurvatureWaveSolver 解决:
    % d²u/dt² = c² d²u/dx² + α u d²u/dx²
    %
    % 该方程是非线性的，其更新规则依赖于当前的 u 值。
    
    properties
        alpha = 1.0 % 非线性项系数
    end
    
    methods
        function obj = NonlinearCurvatureWaveSolver(length, c, alpha, nx, dt, t_steps)
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
            % 基本波方程CFL条件: c*dt/dx <= 1
            % 非线性情况通常需要更严格条件，这里仍以此为参考。
            if obj.c * obj.dt / obj.dx > 1
                error('Stability condition violated: c*dt/dx <= 1 for nonlinear curvature wave.');
            end
        end
        
        function step(obj)
            cdt_dx2 = (obj.c * obj.dt / obj.dx)^2;
            alpha_dt2 = obj.alpha * (obj.dt)^2;
            for i = 2:obj.nx-1
                laplacian = (obj.u(i+1) - 2*obj.u(i) + obj.u(i-1));
                obj.u_next(i) = 2*obj.u(i) - obj.u_prev(i) ...
                    + (cdt_dx2 + alpha_dt2 * obj.u(i)/((obj.dx)^2)) * laplacian;
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
