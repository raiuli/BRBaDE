function [F,CR] = FLCS(d11,d12,d21,d22)

fis = newfis('F_FLC');
%add inputs
fis = addvar(fis,'input','d11',[0 1]);
fis = addmf(fis,'input',1,'small','gaussmf',[0.25 0.05]);
fis = addmf(fis,'input',1,'medium','gaussmf',[0.25 0.5]);
fis = addmf(fis,'input',1,'big','gaussmf',[0.25 0.9]);

fis = addvar(fis,'input','d12',[0 1]);
fis = addmf(fis,'input',2,'small','gaussmf',[0.35 0.01]);
fis = addmf(fis,'input',2,'medium','gaussmf',[0.35 0.5]);
fis = addmf(fis,'input',2,'big','gaussmf',[0.35 0.9]);

fis = addvar(fis,'output','F',[0 2]);
fis = addmf(fis,'output',1,'small','gaussmf',[0.5 0.3]);
fis = addmf(fis,'output',1,'medium','gaussmf',[0.5 0.6]);
fis = addmf(fis,'output',1,'big','gaussmf',[0.5 0.9]);

ruleList = [1 1 1 1 1;
            1 2 2 1 1;
            1 3 3 1 1;
            
            2 1 2 1 1;
            2 2 2 1 1;
            2 3 3 1 1;
            
            3 1 3 1 1;
            3 2 3 1 1;
            3 3 3 1 1];
       
fis = addrule(fis,ruleList);
F=evalfis([d11 d12],fis);

fis = newfis('CR_FLC');
%add inputs
fis = addvar(fis,'input','d21',[0 2]);
fis = addmf(fis,'input',1,'small','gaussmf',[0.5 0.1]);
fis = addmf(fis,'input',1,'medium','gaussmf',[0.5 0.8]);
fis = addmf(fis,'input',1,'big','gaussmf',[0.5 1.5]);

fis = addvar(fis,'input','d22',[0 2]);
fis = addmf(fis,'input',2,'small','gaussmf',[0.5 0.1]);
fis = addmf(fis,'input',2,'medium','gaussmf',[0.5 0.8]);
fis = addmf(fis,'input',2,'big','gaussmf',[0.5 1.5]);

fis = addvar(fis,'output','CR',[0 1]);
fis = addmf(fis,'output',1,'small','gaussmf',[0.35 0.4]);
fis = addmf(fis,'output',1,'medium','gaussmf',[0.35 0.7]);
fis = addmf(fis,'output',1,'big','gaussmf',[0.35 1.0]);

ruleList = [1 1 1 1 1;
            1 2 2 1 1;
            1 3 3 1 1;
            
            2 1 2 1 1;
            2 2 2 1 1;
            2 3 3 1 1;
            
            3 1 3 1 1;
            3 2 3 1 1;
            3 3 3 1 1];
        
fis = addrule(fis,ruleList);
CR=evalfis([d21 d22],fis);
% d11_s = gaussmf(d11,[0.25 0.05])
% d11_m = gaussmf(d11,[0.25 0.5])
% d11_b = gaussmf(d11,[0.25 0.9])
% 
% d12_s = gaussmf(d12,[0.35 0.01])
% d12_m = gaussmf(d12,[0.35 0.5])
% d12_b = gaussmf(d12,[0.35 0.9])
% 
% d21_s = gaussmf(d21,[0.5 0.1])
% d21_m = gaussmf(d21,[0.5 0.8])
% d21_b = gaussmf(d21,[0.5 1.5])
% 
% 
% d22_s = gaussmf(d22,[0.5 0.1])
% d22_m = gaussmf(d22,[0.5 0.8])
% d22_b = gaussmf(d22,[0.5 1.5])
% 
% F_s = gaussmf(F,[0.5 0.3])
% F_m = gaussmf(F,[0.5 0.6])
% F_b = gaussmf(F,[0.5 0.9])
% 
% CR_s = gaussmf(CR,[0.35 0.4])
% CR_m = gaussmf(CR,[0.35 0.7])
% CR_b = gaussmf(CR,[0.35 1.0])
%F=0.8;
%C=0.8;
