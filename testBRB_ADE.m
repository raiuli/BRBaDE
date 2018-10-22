clear;
clc;
addpath(pwd+"/DE");
addpath(pwd+"/BRBADE");
d11=0:0.1:1;
d12=0:0.1:1;
d21=2*d11;
d22=2*d12;
F=[]
CR=[]
D11=[]
D12=[]
for i=1 :length(d11)
    for j=1 :length(d12)
        [newF,newCR]=brb_optimization_customeAWRW(d11(i),d12(j), d21(i),d22(j));
        fprintf(1,'%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f\n',d11(i),d12(j), d21(i),d22(j),newF,newCR);
        F=horzcat(F,newF);
        CR=horzcat(CR,newCR);
        D11=horzcat(D11,d11(i));
        D12=horzcat(D12,d12(j));
    end
end
x=1:1:121;
m=ones(121,1)*0.5;
plot(x,D11,'r--*',x,D12,'g--*',x,CR,'b--*',x,m,'b','LineWidth',1,...
    'MarkerSize',5);
legend('D11','D12','CR');
ruleweight=[1	1	0	0 0	0	0 0 0];
 chr=mat2str(ruleweight);
title(chr)
grid on;
 axis on;
    %set(gca, 'YScale', 'log')
    s=strcat(chr,'[ 1 1]ruleweight_linechart.png');
    saveas(gcf,s);