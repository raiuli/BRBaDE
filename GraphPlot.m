clear;
clc;
format compact;
fclose('all');
fname_v=["test01FunOne","test02FunTwo","test03FunThree","test04FunFour","test05FunFive",...
    "test06FunSix","test07FunSeven","test08FunEight","test09FunNine","test10FunTen"];
title_v=["Test Function 01","Test Function 02","Test Function 03","Test Function 04","Test Function 05",...
    "Test Function 06","Test Function 07","Test Function 08","Test Function 09","Test Function 10"];
% data(iter,1)=PC;
% data(iter,2)=FC;
% data(iter,3)=d11;
% data(iter,4)=d12;
% data(iter,5)=d21;
% data(iter,6)=d22;
% data(iter,7)=F;
% data(iter,8)=CR;
% data(iter,9)=bestval;
% data(iter,10)=iter;
for i=1:10
    
%     if(i==6)
%         continue
%     end    
    s=strcat(fname_v(i),'_brbinsde.mat');
    load(s,'data');
    brbinsde_data=data;
    s=strcat(fname_v(i),'_fade.mat');
    load(s,'data');
    fade_data=data;
    s=strcat(fname_v(i),'_de.mat');
    load(s,'data');
    de_data=data;
    s=strcat(fname_v(i),'_gade.mat');
    load(s,'data');
    gade_data=data;
    s=strcat(fname_v(i),'_pso.mat');
    load(s,'data');
    pso_data=data;
    %plot(brbinsde_data(:,10),brbinsde_data(:,9),'-ok',fade_data(:,10),fade_data(:,9),'--+r',...
    %de_data(:,10),de_data(:,9),'-.*b','markersize',2),xlabel('iterantions'), ylabel('f(x)'),...
    %plot(brbinsde_data(:,10),brbinsde_data(:,9),'-ok',fade_data(:,10),fade_data(:,9),'--+r',de_data(:,10),de_data(:,9),'-.*b',pso_data(:,10),pso_data(:,9),'-.*g'),xlabel('iterantions'), ylabel('f(x)'),...
    %plot(brbinsde_data(:,10),brbinsde_data(:,9),'-k',fade_data(:,10),fade_data(:,9),'--r',de_data(:,10),de_data(:,9),':b',pso_data(:,10),pso_data(:,9),'-.g'),xlabel('iterantions'), ylabel('f(x)'),...
     plot(brbinsde_data(:,10),brbinsde_data(:,9),'-k',...
          fade_data(:,10),fade_data(:,9),'--r', ...
          de_data(:,10),de_data(:,9),':b',...
          gade_data(:,10),gade_data(:,9),'-.g',...
          'LineWidth',2) 
          xlabel('iterantions'), ylabel('f(x)')
         legend('BRBADE', 'FADE','DE','GADE')
         title(title_v(i))
    %legend('BRBADE', 'FADE','DE','PSO'), title(title_v(i))
    
    grid on;
    %axis on;
    %set(gca, 'YScale', 'log')
    s=strcat(fname_v(i),'_linechart.png');
    saveas(gcf,s);
    s=sprintf('%s,%5.5f,%5.5f,%5.5f,%5.5f,%5.5f--%5.5f,%5.5f,%5.5f,%5.5f,%5.5f',...
        title_v(i),brbinsde_data(size(brbinsde_data,1),9),fade_data(size(fade_data,1),9),de_data(size(de_data,1),9),pso_data(size(pso_data,1),9),gade_data(size(gade_data,1),9)...
        ,brbinsde_data(size(brbinsde_data,1),10),fade_data(size(fade_data,1),10),de_data(size(de_data,1),10),pso_data(size(pso_data,1),10),gade_data(size(gade_data,1),10));
    disp(s);
end
% plot(brbinsde_data(:,10),brbinsde_data(:,9),'og','markersize',10);
% hold on
% plot(fade_data(:,10),fade_data(:,9),'ob','markersize',10);
% hold on
% plot(de_data(:,10),de_data(:,9),'*b','markersize',10);
% hold on
% plot(brbinsde_data(:,10),brbinsde_data(:,9),'g',fade_data(:,10),fade_data(:,9),'r',...
%     de_data(:,10),de_data(:,9),'b','markersize',10), legend('brbDe', 'fade','DE');
% grid on;
% axis on;
% set(gca, 'YScale', 'log')
% saveas(gcf,'Barchart01.png')
