classdef BaseWaveSolver < handle
    % BaseWaveSolver 抽象类：定义波动方程数值求解的基本框架和接口
    properties
        length = 1.0      % 绳子的长度
        c = 1.0           % 波速
        nx = 100          % 空间网格点数
        dt = 0.001        % 时间步长
        t_steps = 500     % 时间步数
        dx                % 空间步长
        x                 % 空间坐标
        u                 % 当前时刻波形
        u_prev            % 前一时刻波形
        u_next            % 下一时刻波形
    end
    
    methods
        function obj = BaseWaveSolver(length, c, nx, dt, t_steps)
            if nargin > 0
                obj.length = length;
            end
            if nargin > 1
                obj.c = c;
            end
            if nargin > 2
                obj.nx = nx;
            end
            if nargin > 3
                obj.dt = dt;
            end
            if nargin > 4
                obj.t_steps = t_steps;
            end
            
            obj.dx = obj.length / obj.nx;
            obj.x = linspace(0, obj.length, obj.nx);
            obj.u = zeros(1, obj.nx);
            obj.u_prev = zeros(1, obj.nx);
            obj.u_next = zeros(1, obj.nx);
            
            obj.checkStability();
        end
        
        function setInitialCondition(obj, func, du_dt)
            if nargin < 3
                du_dt = @(x) zeros(size(x));
            end
            obj.u = func(obj.x);
            % 利用 u_prev 来包含初始速度信息:
            % u(t-dt) = u(t) - dt * du/dt
            obj.u_prev = obj.u - obj.dt * du_dt(obj.x);
        end
        
        function runSimulation(obj, output_file)
            if nargin < 2
                output_file = 'wave_animation.gif';
            end
            
            fig = figure;
            ax = axes(fig);
            lineH = plot(ax, obj.x, obj.u, 'LineWidth', 2);
            ax.XLim = [0 obj.length];
            ax.YLim = [-1.1 1.1];
            title(ax, '1D Wave Equation Simulation');
            xlabel(ax, 'x');
            ylabel(ax, 'u(x,t)');
            timeText = text(ax, 0.05, 0.9, '', 'Units', 'normalized');
            
            % 创建GIF动画的第一帧
            frame = getframe(fig);
            [A,map] = rgb2ind(frame2im(frame),256,'nodither');
            imwrite(A,map,output_file,'GIF','LoopCount',Inf,'DelayTime',0.02);
            
            % 循环迭代模拟
            for n = 1:obj.t_steps
                obj.step();
                set(lineH, 'YData', obj.u);
                current_time = n * obj.dt;
                set(timeText, 'String', sprintf('t = %.4f s', current_time));
                
                drawnow;
                frame = getframe(fig);
                [A,map] = rgb2ind(frame2im(frame),256,'nodither');
                imwrite(A,map,output_file,'GIF','WriteMode','append','DelayTime',0.02);
            end
            
            % 将最后一帧保存为PNG图片
            % 再次获取当前帧（此时是最后时刻）
            frame = getframe(fig);
            finalImage = frame2im(frame);
            imwrite(finalImage, 'final_frame.png');
            
            close(fig);
        end

    end
    
    methods (Abstract)
        checkStability(obj)
        step(obj)
        applyBoundaryConditions(obj)
    end
end
