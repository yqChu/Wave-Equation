classdef SimpleWaveSolver < BaseWaveSolver
    % SimpleWaveSolver 实现无耗散一维波动方程:
    % d²u/dt² = c² d²u/dx²
    %
    % 稳定性条件: c*dt/dx <= 1
    
    methods
        function checkStability(obj)
            if obj.c * obj.dt / obj.dx > 1
                error('Stability condition violated: c * dt / dx <= 1');
            end
        end
        
        function step(obj)
            for i = 2:obj.nx-1
                obj.u_next(i) = 2*obj.u(i) - obj.u_prev(i) + ...
                    (obj.c * obj.dt / obj.dx)^2 * (obj.u(i+1) - 2*obj.u(i) + obj.u(i-1));
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
