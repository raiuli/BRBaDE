% Stat clear previous data
clear;
clc;
format compact;
fclose('all');
% % //addpath('/Users/raiuli/Google Drive/MathLabBRBOptimization/MyMatlab/BRBHmodifRefValueParallelAWRWBD_DE/DE');
addpath(pwd+"/DE");
addpath(pwd+"/BRBADE");
%delete(gcp('nocreate'))
%parpool('local',4)
global input outputOpti observedOutput...
    transformedRefVal conseQuentRef...
    rulebase sizeOfData...
    numOfVariables numOfconRefval numOfAttrWeight numOfRuleWeight numOfbeliefDegrees ...
    fid_x1 fid_f1 data_rule Aeq beq d;

formatOut = 'yyyy-mmm-dd_HH_MM_SS';
dateString = datestr(datetime('now'),formatOut);
s = strcat('Log/x1_AWRWBD_',dateString,'.txt');
fid_x1 = fopen (s, 'w');
%fid_x1 = fopen ('x1.txt', 'w');
fid_f1 = fopen ('f1.txt', 'w');
fid_tp = fopen ('trainedParam.txt', 'w');

dateString = datestr(datetime('now'));
fprintf(fid_x1,'Starting program %s \n',dateString);
%read input file
%fid = fopen ('SecDatasetsmf.txt', 'r');
fid = fopen ('inputJBY.txt', 'r');
numberOfInputData=0;
%input=zeros(12,5);
while ~feof(fid)
    numberOfInputData=numberOfInputData+1;
    line=fgetl(fid);
    if numberOfInputData==1
        keySet=split(line,',');
        %valueSet = {ones};
        % mapObj = containers.Map(keySet,valueSet);
        
    else
        valueSet(numberOfInputData-1,:)=str2num(line);
        %svalueSet=str2num(split(line,','));
    end
    
    
end
fclose(fid);
numberOfInputData=numberOfInputData-1;
sizeOfData=numberOfInputData;
%t_data_rule=zeros(1,sizeOfData);
fivid = fopen ('trainedParamAWRWBD.txt', 'r');
numberOfInitialVal=0;
%input=zeros(12,5);
while ~feof(fivid)
    numberOfInitialVal=numberOfInitialVal+1;
    line=fgetl(fivid);
    %initialVal(numberOfInitialVal,:)=str2num(line);
    initialVal{numberOfInitialVal}=str2num(line);
end
fclose(fivid);



keySet=cellstr(keySet);
valueSet=num2cell(valueSet,1);
mapObj = containers.Map(keySet,valueSet);

%For shahadat pue data
brbTree(1).antecedent=cellstr(['x1';'x2']);
brbTree(1).antRefval=[ 17.0 10.0 3.0 ;
                       26.0 23.0 20.0];
brbTree(1).consequent=cellstr('x3');
brbTree(1).conRefval=[2.0 1.5 1];
brbTree(1).rulebaseFile=['rulebaseX3.txt'];

brbTree(1).antecedent=cellstr(['x1';'x2']);
brbTree(1).antRefval=[ 17.0  13.5  10.0  6.5  3.0 ;
                      26.0 24.5 23.0 21.5 20.0];
brbTree(1).consequent=cellstr('x3');
brbTree(1).conRefval=[2.0 1.75 1.5 1.25 1];
brbTree(1).rulebaseFile=['rulebaseX3.txt'];

% brbTree(1).antecedent=cellstr(['x22';'x23']);
% brbTree(1).antRefval=[110 56 2;
%                       11 5.75 0.5];
% brbTree(1).consequent=cellstr('x08');
% brbTree(1).conRefval=[90 70 50];
% brbTree(1).rulebaseFile=['rulebaseX08.txt'];
% brbTree(1).antecedent=cellstr(['x27';'x28';'x29';'x30';'x31']);
% brbTree(1).antRefval=[1 0.5 0;
%                       1 0.5 0;
%                       1 0.5 0;
%                       1 0.5 0;
%                       1 0.5 0];
% brbTree(1).consequent=cellstr('x12');
% brbTree(1).conRefval=[1 0.5 0];
% brbTree(1).rulebaseFile=['rulebaseX12.txt'];
tic;
for brdTreeID=1:size(brbTree,2)
    conseQuentRef=brbTree(brdTreeID).conRefval;
    numOfAttrWeight=size(brbTree(brdTreeID).antRefval,1);
    
    numOfconRefval=size(brbTree(brdTreeID).conRefval,2);
    %read initial rule base for subrule base 1
    %rulebase=readRuleBase(brbTree(brdTreeID).rulebaseFile);
    rule=calculateInitialRulebase(brbTree(brdTreeID).antRefval,brbTree(brdTreeID).conRefval);
    rulebase=struct;
    d=prod(size(brbTree(brdTreeID).antRefval))
    for i=1:size(rule,1)
        rulebase(i).conse=rule(i,size(brbTree(brdTreeID).antRefval,1)+1:end);
        rulebase(i).ruleweight=1;
    end
    %rulebase.
    size(brbTree(brdTreeID).antecedent,1);
    observedOutput=mapObj(brbTree(brdTreeID).consequent{1});
    transformedRefVal={};
    fprintf(fid_x1,'\nAntecedents:');
    fprintf('\nAntecedents:');
    for antecedentID=1:size(brbTree(brdTreeID).antecedent,1)
        fprintf(fid_x1,' %s(',brbTree(brdTreeID).antecedent{antecedentID});
        fprintf(' %s(',brbTree(brdTreeID).antecedent{antecedentID});
        in=mapObj(brbTree(brdTreeID).antecedent{antecedentID,1});
        antcedentRefVal=brbTree(brdTreeID).antRefval(antecedentID,:);
        fprintf(fid_x1,'%2.2f ',antcedentRefVal);
        fprintf(fid_x1,')');
        fprintf('%2.2f ',antcedentRefVal);
        fprintf(')');
        
        tmp=inputTransform(in,antcedentRefVal,numberOfInputData);
        transformedRefVal(antecedentID,:)={tmp};
        %      transformedRefVal(:,antecedentID)=tmp
        attrWeight(antecedentID)=1;
    end
    fprintf(fid_x1,'=>%s (',brbTree(brdTreeID).consequent{1});
    fprintf(fid_x1,'%2.2f ',brbTree(brdTreeID).conRefval);
    fprintf(fid_x1,')\n');
    
    fprintf('=>%s (',brbTree(brdTreeID).consequent{1});
    fprintf('%2.2f ',brbTree(brdTreeID).conRefval);
    fprintf(')\n');
    
    numOfRuleWeight=size(rulebase,2);
    numOfbeliefDegrees=numOfRuleWeight*numOfconRefval;
    %numOfVariables=numOfAttrWeight+numOfRuleWeight+numOfbeliefDegrees;
    numOfVariables=numOfconRefval+numOfAttrWeight+numOfRuleWeight+numOfbeliefDegrees;
    fprintf(fid_x1,'Number of Varaibles: %d=%d(CR)+%d(AW)+%d(RW)+%d(BD)\n',numOfVariables,numOfconRefval,numOfAttrWeight,numOfRuleWeight,numOfbeliefDegrees);
    fprintf('Number of Varaibles: %d=%d(CR)+%d(AW)+%d(RW)+%d(BD)\n',numOfVariables,numOfconRefval,numOfAttrWeight,numOfRuleWeight,numOfbeliefDegrees);
    %initialiaze the x0
    initialValAttrWeight=ones([1,numOfAttrWeight]);
    initialValRuleWeight=ones([1,numOfRuleWeight]);
    initialValConsequent=ones([1,numOfconRefval]);
    betam=[];
    for i=1:numOfRuleWeight
        betam(i,:)=rulebase(i).conse;
    end
    z1=betam';
    betam=z1(:)';
    x0=horzcat(initialValAttrWeight,initialValRuleWeight,betam,initialValConsequent);
    %x0=initialVal{brdTreeID};
    %x0=horzcat(x0,initialValConsequent);
    %initialiaze the constraints
    lb = zeros(1,numOfVariables);
    ub =ones(1,numOfVariables);
    %initialVal{brdTreeID};
    fprintf(fid_x1,'\nIntial value\n');
    fprintf (fid_x1,'Attribute Weights\n');
    fprintf (fid_x1,'%d ', x0(1:numOfAttrWeight) );
    fprintf (fid_x1,'\nRuleWeights\n');
    fprintf (fid_x1,'%d ', x0(numOfAttrWeight+1:numOfAttrWeight+numOfRuleWeight) );
    fprintf (fid_x1,'\nBelief Degrees\n');
    z=x0(numOfAttrWeight+numOfRuleWeight+1:numOfVariables-numOfconRefval);
    j=1;
    for i=1:length(z)/length(conseQuentRef)
        fprintf (fid_x1,'%2.2f ', z(j:j+length(conseQuentRef)-1));
        fprintf (fid_x1,'\n');
        j=j+length(conseQuentRef);
    end
    fprintf (fid_x1,'\nConsequent\n');
    fprintf (fid_x1,'%d ', x0(numOfVariables-numOfconRefval+1:numOfVariables) );
    Aeq=zeros(numOfRuleWeight,numOfbeliefDegrees);
    Aeq1=zeros(numOfRuleWeight,numOfAttrWeight+numOfRuleWeight);
    Aeq2=zeros(numOfRuleWeight,numOfconRefval);
    j=1;
    for i=1:numOfRuleWeight
        Aeq([i],j:j+numOfconRefval-1)=ones(1,numOfconRefval);
        j=j+numOfconRefval;
    end
    Aeq=horzcat(Aeq1,Aeq,Aeq2);
    fprintf(fid_x1,'\nAeq\n');
    %fprintf('%d ',Aeq)
    fprintf(fid_x1,[repmat('%2.2f\t', 1, size(Aeq, 2)) '\n'], Aeq');
    %   Aeq= [0     0     0     0     1     1     1     0     0     0     0     0     0;
    %       0     0     0     0     0     0     0     1     1     1     0     0     0;
    %       0     0     0     0     0     0     0     0     0     0     1     1     1]
    beq = ones(numOfRuleWeight,1);
    
    A=zeros(numOfconRefval-1,numOfVariables);
    %A(numOfVariables-numOfconRefval:numOfVariables)
%     A(1,numOfVariables-numOfconRefval+1:numOfVariables)=[1 -1 0 ];
%     A(2,numOfVariables-numOfconRefval+1:numOfVariables)=[0 1 -1 ];
    
    
    A(1,numOfVariables-numOfconRefval+1:numOfVariables)=[1 -1 0 0 0];
    A(2,numOfVariables-numOfconRefval+1:numOfVariables)=[0 1 -1 0 0];
    A(3,numOfVariables-numOfconRefval+1:numOfVariables)=[0 0  1 -1 0];
    A(4,numOfVariables-numOfconRefval+1:numOfVariables)=[0 0  0  1 -1];
    A(5,numOfVariables-numOfconRefval+1:numOfVariables)=[0 0  1 -1 0];
    B=zeros(numOfconRefval-1,1);
    save('dumpGlobalVariable.mat','input', 'outputOpti', 'observedOutput',...
        'transformedRefVal', 'conseQuentRef', 'rulebase', 'sizeOfData',...
        'numOfVariables', 'numOfconRefval', 'numOfAttrWeight', 'numOfRuleWeight', 'numOfbeliefDegrees','Aeq','beq','d');
    
    BRBES.input=input;
    BRBES.outputOpti=outputOpti
    display('starting Optimization------------------------------------------')
    tic
    % implementing DE
    
    % Initialization and run of differential evolution optimizer.
    
    
    % VTR		"Value To Reach" (stop when ofunc < VTR)
    VTR = 1.e-6;
    
    % D		number of parameters of the objective function
    D = numOfVariables;
    
    % XVmin,XVmax   vector of lower and bounds of initial population
    %    		the algorithm seems to work well only if [XVmin,XVmax]
    %    		covers the region where the global minimum is expected
    %               *** note: these are no bound constraints!! ***
    
    XVmin = lb;
    XVmax = ub;
    % y		problem data vector (remains fixed during optimization)
    y=[];
    
    % NP            number of population members
    %NP = 10*numOfVariables;
    NP = 500;
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
    
    strategy = 1;
    
    % refresh       intermediate output will be produced after "refresh"
    %               iterations. No intermediate output will be produced
    %               if refresh is < 1
    refresh = 1;
    tic
    [x,f,nf] = BRBADEcustomeRWAWv3('objFunAllParallel',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
    %[x,f,nf] = deEva('objFunAllParallel',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
    toc
    %defEva(fname              ,VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
    %[x,f,nf] = deAckely('Ackley',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)
    
    %fsurf(objFunAll(x0),[0,0],'ShowContours','on')
    %objFunAll(x0)
    %sqp
    %options = optimoptions('fmincon','Display','iter','Algorithm','sqp','PlotFcn',{@optimplotx,...
    %    @optimplotfval,@optimplotfirstorderopt});
    % options = optimoptions('fmincon','Display','iter','Algorithm','sqp','UseParallel',true);
    %options = optimoptions('fmincon','Display','iter','PlotFcn',{@optimplotx,...
    %    @optimplotfval,@optimplotfirstorderopt});
    
    %options.MaxFunctionEvaluations=6000
    %options.MaxIterations=60;
    % [ x, fval, exitflag, output ]=fmincon ( @objFunAllParallel, x0, ...
    %      A, B, Aeq, beq,lb, ub,emptyNolinearConstraints,options);
    
    %   [ x, fval, exitflag, output ]=fmincon ( @objFunAll, x0, ...
    %       [], [], Aeq, beq,lb, ub,option
    %    if size(x,2)==numOfVariables
    %objFunAll(x)
    t=toc;
    
    objFunAlleva(x);
    % t_data_rule=vertcat(t_data_rule,data_rule);
    x0=x;
    fprintf(fid_tp,'%2.2f ',x0);
    fprintf(fid_tp,'\n');
    
    mapObj(brbTree(brdTreeID).consequent{1})=outputOpti;
    fprintf(fid_x1,'\nOptimizied value\n');
    fprintf (fid_x1,'Attribute Weights\n');
    fprintf (fid_x1,'%2.2f ', x0(1:numOfAttrWeight) );
    fprintf (fid_x1,'\nRuleWeights\n');
    fprintf (fid_x1,'%2.2f ', x0(numOfAttrWeight+1:numOfAttrWeight+numOfRuleWeight) );
    fprintf (fid_x1,'\nBelief Degrees\n');
    z=x0(numOfAttrWeight+numOfRuleWeight+1:numOfVariables);
    j=1;
    for i=1:length(z)/length(conseQuentRef)
        fprintf (fid_x1,'%2.2f ', z(j:j+length(conseQuentRef)-1));
        fprintf (fid_x1,'\n');
        j=j+length(conseQuentRef);
    end
    %         s=sprintf('  aw%d',[1:numOfAttrWeight]);
    %         %s=strcat(s,{' ');
    %         s=strcat(s, sprintf('  rw%d',[1:numOfRuleWeight]));
    %         %s=strcat(s,{' '});
    %         s=strcat(s, sprintf('  bd%d',[1:numOfbeliefDegrees]));
    fprintf (fid_x1,'Output vals:');
    fprintf (fid_x1,'%2.5f ', outputOpti );
    %         fprintf (fid_x1,'\nfval: %2.5f \n', fval );
    %          fprintf (fid_x1,'iterations: %2.5f \n', output.iterations );
    %          fprintf (fid_x1,'funcCount: %2.5f \n', output.funcCount );
    fprintf(fid_x1,'Elapsed time is %2.4f sec\n',t);
    fprintf('Elapsed time is %2.4f sec\n',t);
    %         fprintf ( fid_x1,'%s\n', s );
    %         fprintf ( fid_x1,'%2.2f ', x );
    %         fprintf ( fid_x1,'\n');
    
    %      else
    %         display('Variable  do not match------------------------------------------')
    %     end
end
dateString = datestr(datetime('now'));
fprintf(fid_x1,'Ending program %s',dateString);
fprintf('Ending program %s',dateString);
% fprintf ( fid_x1,'____________________________\n');
% fprintf ( fid_x1,[repmat('%2.2f\t', 1, size(t_data_rule, 2)) '\n'], t_data_rule );
fclose(fid_x1);
fclose(fid_f1);
fclose(fid_tp);
