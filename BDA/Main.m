%___________________________________________________________________%
%  Binary Dragonfly Algorithm (BDA) source codes demo version 1.0   %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper:                                                     %
%                                                                   %
%   S. Mirjalili, Dragonfly algorithm: a new meta-heuristic         %
%   optimization technique for solving single-objective, discrete,  %
%   and multi-objective problems, Neural Computing and Applications % 
%   DOI: http://dx.doi.org/10.1007/s00521-015-1920-1                %
%___________________________________________________________________%

clear all 
close all
clc

CostFunction=@(x) MyCost(x); % Modify or replace Mycost.m according to your cost funciton

Max_iteration=500; % Maximum number of iterations
N=60; % Number of particles
nVar=40;

% BDA with a v-shaped transfer function
[Best_pos, Best_score ,Convergence_curve]=BDA(N, Max_iteration, nVar, CostFunction);

plot(Convergence_curve,'DisplayName','BDA','Color', 'r');
hold on

title(['Convergence curve']);
xlabel('Iteration');ylabel('Average Best-so-far');
legend('BDA',1);
box on
axis tight

display(['The best solution obtained by BDA is : ', num2str(Best_pos')]);
display(['The best optimal value of the objective funciton found by BDA is : ', num2str(Best_score)]);

