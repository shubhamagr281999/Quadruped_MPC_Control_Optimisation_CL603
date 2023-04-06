clc
clear all
k = 1;
xinitial = zeros(24*k,1);
xinitial([1:36],1)=0.1;
[xfinal,x] = Active_set(xinitial);