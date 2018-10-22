clear;
clc;
format compact;
a=12.45566
x=strfind(num2str(a), '.')
anew=a*10^-(x-2)
a=9.45566
x=strfind(num2str(a), '.')
anew=a*10^-(x-2)
a=91239.45566
x=strfind(num2str(a), '.')
anew=a*10^-(x-2)

b=0.13456
x=strfind(num2str(b), '.')
bnew=b*10^(x-1)