classdef DampedWaveSolver < BaseWaveSolver
    % DampedWaveSolver 实现带耗散项的波动方程:
    % d²u/dt² + γ du/dt = c² d²u/dx²
    %
    % 在此示例中，我们使用类似于无耗散方程的稳定性条件 c*dt/dx <= 1。
    % 实际稳定性条件可能更复杂，但本示例为演示用。
    
    properties
        gamma = 0.1 % 耗散强度，γ > 0
    end
    
    methods
        function obj = DampedWaveSolver(length, c, gamma, nx, dt, t_steps)
            if nargin >= 1; obj.length = length; end
            if nargin >= 2; obj.c = c; end
            if nargin >= 3; obj.gamma = gamma; end
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
            % 与无耗散情况类似，我们仍要求 c*dt/dx <= 1 作为基本要求
            if obj.c * obj.dt / obj.dx > 1
                error('Stability condition violated: c * dt / dx <= 1 for damped wave');
            end
        end
        
        function step(obj)
            % 带耗散项的时间步进更新公式:
            % u^(n+1)_i = 2u^(n)_i - u^(n-1)_i
            %           + (c² (dt/dx)²)(u^(n)_{i+1}-2u^(n)_i+u^(n)_{i-1})
            %           - γ dt (u^(n)_i - u^(n-1)_i)
            for i = 2:obj.nx-1
                obj.u_next(i) = (2*obj.u(i) - obj.u_prev(i)) ...
                    + (obj.c * obj.dt / obj.dx)^2 * (obj.u(i+1) - 2*obj.u(i) + obj.u(i-1)) ...
                    - obj.gamma * obj.dt * (obj.u(i) - obj.u_prev(i));
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
