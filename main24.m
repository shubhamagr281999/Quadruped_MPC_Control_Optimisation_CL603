clc
clear
Ai=zeros(12,12,10);
Bi=zeros(12,13,10);
for i =1:10
Ai(:,:,i)=eye(12,12);
Bi(:,:,i)=eye(12,13);
end
x0 = zeros(12,1); 
Xref = zeros(12*3,1);
% Xref(1:12)=20;
% Xref(13:24)=32;
% Xref(25:36)=50;
n_steps = 3;
mu = 0.6;
fmin = 10;
fmax = 660;
Q = eye(12,12)*100;
R = eye(12,12)*10^-2;
onGround = ones(4,3);
onGround(2,1)=0;
onGround(3,2)=0;
onGround(4,3)=0;
load("A_B.mat")
[a,b,c]  = Active_set(x0,Xref,n_steps,mu,fmin,fmax,Q,R,onGround,Bi,Ai);