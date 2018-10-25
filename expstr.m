function [v n]=expstr(x)
v=x*10.^floor(1-log10(abs(x)))  ;
n=floor(log10(abs(x)));