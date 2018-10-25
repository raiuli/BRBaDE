%parpool(10)
%parfor i=1:3, c(:,i) = eig(rand(1000)); end
clear;
clc;
format compact;
fclose('all');

addpath(pwd+"/DE");
addpath(pwd+"/BRBADE");
addpath(pwd+"/Fuzzy");
addpath(pwd+"/LSHADE-EpSin");
addpath(pwd+"/GADE");
%set variables
VTR = 1.e-6; 
D_v(1) = 50; 
D_v(2) = 50;
D_v(3) = 50;
D_v(4) = 30;
D_v(5) = 2;
D_v(6) = 50; 
D_v(7) = 50;
D_v(8) = 2;
D_v(9) = 2;
D_v(10) = 50;

NP_v(1) = 500; 
NP_v(2) = 500;
NP_v(3) = 500;
NP_v(4) = 300;
NP_v(5) = 20;
NP_v(6) = 500; 
NP_v(7) = 500;
NP_v(8) = 20;
NP_v(9) = 20;
NP_v(10) = 500;
 
itermax_v(1) = 5000; 
itermax_v(2) = 5000;
itermax_v(3) = 5000;
itermax_v(4) = 5000;
itermax_v(5) = 100;
itermax_v(6) = 5000; 
itermax_v(7) = 5000;
itermax_v(8) = 200;
itermax_v(9) = 50;
itermax_v(10) = 5000;

fname_v=["test01FunOne","test02FunTwo","test03FunThree","test04FunFour","test05FunFive",...
    "test06FunSix","test07FunSeven","test08FunEight","test09FunNine","test10FunTen"];
XVmin = [-2 -1]; XVmax = [2 3];

XVmin_v(1,:)= -5.12;XVmax_v(1,:)= 5.12;
XVmin_v(2,:)= -2.048;XVmax_v(2,:)= 2.047;        
XVmin_v(3,:)= -5.12;XVmax_v(3,:)= 5.12;        
XVmin_v(4,:)= -1.28;XVmax_v(4,:)= 1.28;        
XVmin_v(5,:)= -65.536;XVmax_v(5,:)= 65.536;        
XVmin_v(6,:)= -5.12;XVmax_v(6,:)= 5.12;        
XVmin_v(7,:)= -32.768;XVmax_v(7,:)= 32.768;        
XVmin_v(8,:)= -100;XVmax_v(8,:)= 100;        
XVmin_v(9,:)= -2;XVmax_v(9,:)= 2;        
XVmin_v(10,:)= -600;XVmax_v(10,:)= 600;        
        
F = 0.9;
CR = 0.8; 
strategy = 7;
refresh = 100;
tic
y=[];
for i = 1:1
    D=D_v(i);
    NP=NP_v(i);
    itermax=itermax_v(i);
    XVmin=XVmin_v(i).*ones(D,1)';
    XVmax=XVmax_v(i).*ones(D,1)';
    disp(fname_v(i));
    disp('BRBinsDE');
    [x,f,nf] = BRBADEcustomeRWAW(fname_v(i),VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh);
    
    %[x,f,nf] = BRBinsDE_orig(fname_v(i),VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh);
  % strategy = 10;
   
   %[x,f,nf] = BRBinsDE_CPbest(fname_v(i),VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh);
   %disp('GADE');
   [x,f,nf] = GADE(fname_v(i),VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh);
    
%     
    % disp('FADE');
     %[x,f,nf] = FADE(fname_v(i),VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh);
     %disp('DE');
%     
     %[x,f,nf] = DiffEv(fname_v(i),VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh);
     %disp('Particle Swarm');
     lb = XVmin;
     ub = XVmax;
     nvars = D;
     fun=str2func(fname_v(i));
    options = optimoptions('particleswarm','SwarmSize',25,'HybridFcn',@fmincon,'Display','iter','OutputFcn',@pswplotranges);

    %options = optimoptions('particleswarm','SwarmSize',25,'Display','iter','OutputFcn',@pswplotranges, ...
    %    'FunctionTolerance',VTR,'MaxIterations',itermax);
   % [x,fval,exitflag,output] = particleswarm(fun,nvars,lb,ub,options)
    
    %% Problem Definiton

%  problem.CostFunction = str2func(fname_v(i));  % Cost Function
%  problem.nVar = D;       % Number of Unknown (Decision) Variables
%  problem.VarMin =  XVmin;  % Lower Bound of Decision Variables
%  problem.VarMax =  XVmax;   % Upper Bound of Decision Variables

%% Parameters of PSO

%  params.MaxIt = itermax;        % Maximum Number of Iterations
%  params.nPop = 25;           % Population Size (Swarm Size)
%  params.w = 0.6;               % Intertia Coefficient
%  params.wdamp = 0.99;        % Damping Ratio of Inertia Coefficient
%  params.c1 = 1.8;              % Personal Acceleration Coefficient
%  params.c2 = 1.8;              % Social Acceleration Coefficient
%  params.ShowIterInfo = true; % Flag for Showing Iteration Informatin
%   out = PSO(problem, params);
    ss=sprintf('pso_data.mat');
    load(ss,'data');
    ss=strcat(fname_v(i),'_pso.mat');
    %data=out.result;
    save(ss,'data');
    %disp('LSHADE-EpSin');
     %[x,f,nf] = LSHADE_EpSin(fname_v(i),VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh);

end
toc