% Stat clear previous data
clear;
clc;
format compact;
fclose('all');

% filename = 'SecDataset.xls';
% [num,txt,raw]= xlsread(filename);
% for i=1: size(num,2)
%     [h,stats]=cdfplot(num(:,i))
%     filename = strcat('Barchart',num2str(i),'.png');
%     hold on
%     title(strcat('Empirical CDF',txt(i)));
%     hold off
%     saveas(gcf,filename);
% end


filename = 'JISC_Dataset_Paper_refined-2a.xls';
[num,txt,raw]= xlsread(filename);
corrcoef(num)
for i=1: size(num,2)
    [h,stats]=cdfplot(num(:,i))
    ecdf(num(:,i),'alpha',0.05,'Bounds','on')
    filename = strcat('Barchart',num2str(i));
    hold on
    title(strcat('Empirical CDF',txt(i)));
    hold off
    saveas(gcf,filename,'epsc');
end
y = datasample(num,100)
csvwrite('test.csv',y);
for i=1: size(y,2)
    [h,stats]=cdfplot(num(:,i))
    ecdf(num(:,i),'alpha',0.05,'Bounds','on')
    filename = strcat('Barchart_sample',num2str(i));
    hold on
    title(strcat('Empirical CDF',txt(i)));
    hold off
    saveas(gcf,filename,'epsc');
end
for i=1: size(num,2)
    find(y==num(1,:))
end