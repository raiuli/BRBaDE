clear;
clc;
format compact;
fclose('all');

addpath(pwd+"/DE");
addpath(pwd+"/BRBADE");
addpath(pwd+"/Fuzzy");
formatOut = 'yyyy-mmm-dd_HH_MM_SS';
dateString = datestr(datetime('now'),formatOut);
filename = strcat('randomConf',dateString,'.mat');
load('randomConf2018-Oct-19_18_32_46.mat', 's');
rng(s);
%s=rng();
%save(filename,'s')

%addpath(pwd+"/logging4matlab")
% Initialization and run of differential evolution optimizer.
% A simpler version with fewer explicit parameters is in run0.m
%
% Here for Rosenbrock's function
% Change relevant entries to adapt to your personal applications
%
% The file ofunc.m must also be changed 
% to return the objective function
%

% VTR		"Value To Reach" (stop when ofunc < VTR)
		VTR = 1.e-6; 
        %VTR = 1.e-60; 
% D		number of parameters of the objective function 
		D = 50; 
        %D=2;
% XVmin,XVmax   vector of lower and bounds of initial population
%    		the algorithm seems to work well only if [XVmin,XVmax] 
%    		covers the region where the global minimum is expected
%               *** note: these are no bound constraints!! ***
		XVmin = [-2 -1]; 
		XVmax = [2 3];
        XVmin= -5.12.*ones(D,1)';
        XVmax= 5.12.*ones(D,1)';
        
        XVmin= -100*ones(D,1)';
        XVmax= 100*ones(D,1)';
        XVmin= -200*ones(D,1)';
        XVmax= 200*ones(D,1)';
        
        XVmin= -600*ones(D,1)';
        XVmax= 600*ones(D,1)';
%         XVmin = [-10 -1]; 
% 		XVmax = [10 1];
% y		problem data vector (remains fixed during optimization)
		y=[]; 

% NP            number of population members
		NP = 100; 

% itermax       maximum number of iterations (generations)
		itermax = 5000; 

% F             DE-stepsize F ex [0, 2]
		F = 0.8; 

% CR            crossover probabililty constant ex [0, 1]
		CR = 0.8; 

% strategy       1 --> DE/best/1/exp           6 --> DE/best/1/bin
%                2 --> DE/rand/1/exp           7 --> DE/rand/1/bin
%                3 --> DE/rand-to-best/1/exp   8 --> DE/rand-to-best/1/bin
%                4 --> DE/best/2/exp           9 --> DE/best/2/bin
%                5 --> DE/rand/2/exp           else  DE/rand/2/bin

		strategy = 7;

% refresh       intermediate output will be produced after "refresh"
%               iterations. No intermediate output will be produced
%               if refresh is < 1
refresh = 5; 
%[X,Y]=meshgrid(-2:0.1:3);
%Z=100*(Y-X.^2).^2+(ones(size(X))-X).^2;
%Z=100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%figure(1)
%surf(X,Y,Z)
%hold on
%contour(X,Y,Z)
%hold on
%[x,f,nf] = BRBinsDE('rosen',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)

%[x,f,nf] = BRBinsDE('testFunOne',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[x,f,nf] = BRBinsDE('testFunThree',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[x,f,nf] = BRBinsDE('testFunSix',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[x,f,nf] = BRBinsDE('testFunSeven',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[x,f,nf] = BRBinsDE('testFunEight',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[x,f,nf] = BRBinsDE('testFunNine',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)

[x,f,nf] = BRBADEv3('testFunTen',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%load('testFunTen_brbinsde.mat','data');
%brbinsde_data=data;
%plot3(x(1,1),x(1,2),f,'or','markersize',10);
%plot(x(1,1),x(1,2),'or','markersize',10);
%grid on;
%hold on
%[x,f,nf] = FADE('rosen',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[x,f,nf] = FADE('testFunOne',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[x,f,nf] = FADE('testFunThree',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[x,f,nf] = FADE('testFunSix',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[x,f,nf] = FADE('testFunSeven',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[x,f,nf] = FADE('testFunEight',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[x,f,nf] = FADE('testFunNine',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[x,f,nf] = FADE('testFunTen',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%load('testFunTen_fade.mat','data');
%fade_data=data;
%[y,f,nf] = devec3('rosen',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[y,f,nf] = devec3('testFunOne',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[y,f,nf] = devec3('testFunThree',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[y,f,nf] = devec3('testFunSix',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[y,f,nf] = devec3('testFunSeven',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[y,f,nf] = devec3('testFunEight',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[y,f,nf] = devec3('testFunNine',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%[y,f,nf] = devec3('testFunTen',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
%load('testFunTen_de.mat','data');
%de_data=data;
%plot3(x(1,1),x(1,2),f,'og','markersize',10);
display('finished')
figure 1
plot(brbinsde_data(:,10),brbinsde_data(:,9),'og','markersize',10);
hold on
plot(fade_data(:,10),fade_data(:,9),'ob','markersize',10);
hold on
saveas(gcf,'Barchart.png')
x = [0 : 0.01: 10];
y = sin(x);
g = cos(x);
plot(x, y, x, g, '.-'), legend('Sin(x)', 'Cos(x)')
