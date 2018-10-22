% fis = readfis('tipper.fis')
% getfis(fis)
% getfis(fis,'Inlabels')
% getfis(fis,'input',1)
%getfis(fis,'output',1)

fis = newfis('tipper');
%add inputs
fis = addvar(fis,'input','service',[0 10]);
fis = addmf(fis,'input',1,'poor','gaussmf',[1.5 0]);
fis = addmf(fis,'input',1,'good','gaussmf',[1.5 5]);
fis = addmf(fis,'input',1,'excellent','gaussmf',[1.5 10]);

fis = addvar(fis,'input','food',[0 10]);
fis = addmf(fis,'input',2,'rancid','trapmf',[-2 0 1 3]);
fis = addmf(fis,'input',2,'delicious','trapmf',[7 9 10 12]);

fis = addvar(fis,'output','tip',[0 30]);
fis = addmf(fis,'output',1,'cheap','trimf',[0 5 10]);
fis = addmf(fis,'output',1,'average','trimf',[10 15 20]);
fis = addmf(fis,'output',1,'generous','trimf',[20 25 30]);

ruleList = [1 1 1 1 2;
            2 0 2 1 1;
            3 2 3 1 2];
fis = addrule(fis,ruleList);
evalfis([1 2],fis)