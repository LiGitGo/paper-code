%% Problem parameters
clc;clear all;close all
M=512;
N=512;
MN = M*N;

% Addpath
addpath('functions')
addpath('algorithm_general') 
addpath('algorithm_preCondition') 
addpath('algorithm_linesearch') 
addpath('test') 

% load a test problem 
[A,b,x_true,theta,~,R,d] = fancurvedtomo(N, 0:10:359);
% Remove the all-zero rows from A along with the corresponding rows in b.
nonZeroRowsA = any(A, 2);
A = A(nonZeroRowsA, :);
b=b(nonZeroRowsA, :);
[rA,cA]=size(A);

% matrix B is the discrete gradient operator
[B1,B2] = generate_B_Neumann(M,N);
B=[B1;B2];
mu=1;
% Remove the all-zero rows from B
nonZeroRowsB = any(B, 2);
B = B(nonZeroRowsB, :);


%% get the algorithm step size
K=[A;B];
L=normest(K); 

ga_pre=1/10;
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


%% get f_optval
load 'optval512.mat'; 

