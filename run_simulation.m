% run_simulation.m
% 使用 SimpleWaveSolver 来进行模拟，并生成GIF动画
clear all; close all; clc;
% 定义初始条件
gaussian_func = @(x) exp(-100*(x - 0.5).^2);
zero_velocity = @(x) zeros(size(x));

% 创建求解器实例
solver = SimpleWaveSolver(1.0, 1.0, 200, 0.001, 500);
solver.setInitialCondition(gaussian_func, zero_velocity);
solver.runSimulation('simple_wave');

% run_simulation.m (更新内容)
% 使用 DampedWaveSolver 来进行模拟
damped_solver = DampedWaveSolver(1.0, 1.0, 2.0, 200, 0.001, 500);
damped_solver.setInitialCondition(gaussian_func, zero_velocity);
damped_solver.runSimulation('damped_wave');
% 
solver = PotentialWaveSolver(1.0, 1.0, 50.0, 200, 0.001, 500);
solver.setInitialCondition(gaussian_func, zero_velocity);
solver.runSimulation('potential_wave');

solver = NonlinearPotentialWaveSolver(1.0, 1.0, 50.0, 200, 0.001, 500);
solver.setInitialCondition(gaussian_func, zero_velocity);
solver.runSimulation('nonlinear_potential_wave');

solver = NonlinearAdvectionWaveSolver(1.0, 1.0, 10.0, 200, 0.001, 500);
solver.setInitialCondition(gaussian_func, zero_velocity);
solver.runSimulation('nonlinear_advection_wave');

solver = NonlinearCurvatureWaveSolver(1.0, 1.0, 1.0, 200, 0.001, 500);
solver.setInitialCondition(gaussian_func, zero_velocity);
solver.runSimulation('nonlinear_curvature_wave');
% This seems to be not physical, this term break down the wave.