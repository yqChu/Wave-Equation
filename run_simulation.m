% run_simulation.m
% 使用 SimpleWaveSolver 来进行模拟，并生成GIF动画
clear all; close all; clc;
% 定义初始条件
gaussian_func = @(x) exp(-100*(x - 0.5).^2);
zero_velocity = @(x) zeros(size(x));

% 创建求解器实例
solver = SimpleWaveSolver(1.0, 1.0, 200, 0.001, 200);
solver.setInitialCondition(gaussian_func, zero_velocity);

% 运行模拟并生成GIF
solver.runSimulation('simple_wave');

% run_simulation.m (更新内容)
% 使用 DampedWaveSolver 来进行模拟
damped_solver = DampedWaveSolver(1.0, 1.0, 2.0, 200, 0.001, 200);
damped_solver.setInitialCondition(gaussian_func, zero_velocity);
damped_solver.runSimulation('damped_wave');