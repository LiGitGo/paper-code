%% 
clear all;clc;close all;
addpath('algorithm_general') 
addpath('algorithm_preCondition') 
addpath('algorithm_linesearch') 
addpath('utils') 
addpath('test') 
rng('default')
rng(100);

%% Problem parameters
p=200;
q=2000;
s=10;
u=0.1;

% Generate K
K = randn(p, q);
[U, S, V] = svd(K);
S(1, 1) = 1500;   %100,300,500,1000,1500,2000
K = U * S * V';

% Generate b
xstar=20*rand(size(K,2),1)-10;  
zeroIndex=randperm(size(K,2),size(K,2)-s); 
xstar(zeroIndex)=0; 
nu=sqrt(0.1) * randn(size(K,1),1);
b=K*xstar+nu;

%% get the algorithm step size

L=normest(K,1e-6);

ga_pre=5;
a=(2-sqrt(2))/2; 
ga=ga_pre;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
step=h_a;   
gama1=sqrt(1/step)/ga;
gama2=sqrt(1/step)*ga;
M1=gama1*sum(abs(K))';
M2=gama2*sum(abs(K),2);
tauM=1./M1;
sigmaM=1./M2;
ga_gen=sqrt(mean(tauM)/mean(sigmaM))




%% get avg(tau)*avg(sigma)
format long;
mean(tauM)*mean(sigmaM)

a=(2-sqrt(2))/2;
ga=ga_gen;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
tau = sqrt(h_a)/L*ga;
sigma = sqrt(h_a)/L/ga;
tau*sigma


