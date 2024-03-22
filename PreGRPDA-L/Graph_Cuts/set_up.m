%%
clc;clear all;close all;
addpath('functions')
addpath('toolbox_general')
addpath('toolbox_signal')
addpath('utils')
addpath('algorithm_gen')
addpath('algorithm_linesearch')
%%
name = 'blue_flower';
% name = 'yellow_flower';
I = rescale(load_image(name));

m = size(I,1);
n =size(I,2);
figure(1)
imshow(I);

Tol=1e-7; %1e-8

%% Compute w0 and w1. Compute and display the segmentation
if strcmp(name,'blue_flower')  
    c0 = [0;0;1];
    c1 = [0;1;0];
elseif strcmp(name,'yellow_flower')
    c0 = [1;1;0];
    c1 = [0;0;0];     
end

compute_w = @(I,c)sum( (I - repmat(reshape(c,[1 1 3]), [m n 1])).^2, 3);
w0 = compute_w(I,c0);
w1 = compute_w(I,c1);
w = w0-w1;

beta=10;
Diag=zeros(m,n,2);

Diag(:,:,1)=exp(-beta*(abs(I([2:m,1],:,1)-I(:,:,1))+...
    abs(I([2:m,1],:,2)-I(:,:,2))+abs(I([2:m,1],:,3)-I(:,:,3))));
Diag(:,:,2)=exp(-beta*(abs(I(:,[2:n,1],1)-I(:,:,1))+...
    abs(I(:,[2:n,1],2)-I(:,:,2))+abs(I(:,[2:n,1],3)-I(:,:,3))));

lambda = 1;
alpha = 1/lambda;

%%
[B1,B2]=generate_B_Neumann(m,n);
B=[B2;B1];
d1=reshape(Diag(:,:,1),[m*n,1]); %w^{b,1}
d2=reshape(Diag(:,:,2),[m*n,1]); %w^{b,2}
w_r=reshape(w,[m*n,1]);

Diag1=diag(sparse(d1));
Diag2=diag(sparse(d2));

K=[Diag1*B2;Diag2*B1];
% Remove the all-zero rows from K
nonZeroRowsK = any(K, 2);
K = K(nonZeroRowsK, :);
L=normest(K);

%% obtain the optval by cvx
% cvx_begin 
%  variable x(m*n)
%  minimize norm(K*x,1)+alpha*sum(x.*w_r)
%  subject to
%           0<=x<=1
% cvx_end 

if strcmp(name,'blue_flower')  
    optval=-1.651641644308036e+05;
elseif strcmp(name,'yellow_flower')
    optval=-1.541507760851508e+05;
end

%% get the algorithm step size
ga_pre=10;
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





