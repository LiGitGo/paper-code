clc; 
clear; 
close all;
addpath Images;
addpath algorithm;
addpath utils;
rng(100)

ImageName = 'housergb.png';
% ImageName = 'fruits.png';

I = double(imread(ImageName))/255;   n = length(I); 
figure(1);imshow(I);
z= imnoise(I,'gaussian',0,0.0004); 
if strcmp(ImageName,'housergb.png')  
    S =double(imread('mask.bmp'))/255; S =floor(S); 
elseif strcmp(ImageName,'fruits.png')
    S = ones(n);
    rowsToZero = randperm(n, round(0.8*n));
    S(rowsToZero, :) = 0; 
end
A(:,:,1) = S; A(:,:,2) = S; A(:,:,3) = S; S =A;
z =S.*z;
figure(2);imshow(z);

z=reshape(z, [], 3);
I=reshape(I, [], 3);

%%get A: discrete gradient operator
[A1,A2] = generate_B_Neumann(n,n);
A=[A1;A2];
L=normest(A);

%get B: mask operator
B(:,1)=reshape(S(:,:,1), n*n,1);
B(:,2)=reshape(S(:,:,2), n*n,1);
B(:,3)=reshape(S(:,:,3), n*n,1);

%Shared parameters
Tol = 1e-4;   
lamda=100;

