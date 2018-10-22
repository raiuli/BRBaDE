% Stat clear previous data
clear;
clc;
format compact;
fclose('all');

addpath(pwd+"/DE");
addpath(pwd+"/BRBADE");
%load('randomConf.mat','s')
formatOut = 'yyyy-mmm-dd_HH_MM_SS';
dateString = datestr(datetime('now'),formatOut);
%filename = strcat('randomConf',dateString,'.mat');
%load('randomConf2018-Oct-04_12_05_55.mat', 's');
%rng(s);
%s=rng();
%save(filename,'s')
delete(gcp('nocreate'))
parpool('local',6)
global input outputOpti observedOutput...
    transformedRefVal conseQuentRef...
    rulebase sizeOfData...
    numOfVariables numOfconRefval numOfAttrWeight numOfRuleWeight numOfbeliefDegrees ...
    fid_x1 fid_f1 data_rule Aeq beq d;

s = strcat('Log/simulation_result',dateString,'.txt');
fid_x1 = fopen (s, 'w');

fid_f1 = fopen ('Log/f1.txt', 'w');

dateString = datestr(datetime('now'));
fprintf(fid_x1,'Starting program %s \n',dateString);
%read input file
%fid = fopen ('JISC_Dataset_Paper_refined-2.csv', 'r');
fid = fopen ('SecDatasetstiny.txt', 'r');

%fid = fopen ('inputJBY.txt', 'r');
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
        
        allvalueSet(numberOfInputData-1,:)=str2num(line);
        %svalueSet=str2num(split(line,','));
    end
    
    
end
fclose(fid);

indices = crossvalind('Kfold',allvalueSet(:,4),5);

for counter =1:5
    test = (indices == counter);
    train = ~test;
    
    valueSet=allvalueSet(train,:);
    sizeOfData=size(valueSet,1);
    numberOfInputData=sizeOfData;
    keySet=cellstr(keySet);
    valueSet=num2cell(valueSet,1);
    mapObj = containers.Map(keySet,valueSet);
    %numberOfInputData=numberOfInputData-1;
    %sizeOfData=numberOfInputData;
    brbTree(1).antecedent=cellstr(['x2';'x3';'x4']);
    brbTree(1).antRefval=[17 10 3;
        27 23 21;
        1188897 608186 27475];
    
    brbTree(1).consequent=cellstr('x1');
    brbTree(1).conRefval=[5 2 0];
    
    fid_tp = fopen ('Log/trainedParam.txt', 'a');
    
    
    for brdTreeID=1:size(brbTree,2)
        conseQuentRef=brbTree(brdTreeID).conRefval;
        numOfAttrWeight=size(brbTree(brdTreeID).antRefval,1);
        
        numOfconRefval=size(brbTree(brdTreeID).conRefval,2);
        %read initial rule base for subrule base 1
        %rulebase=readRuleBase(brbTree(brdTreeID).rulebaseFile);
        rule=calculateInitialRulebase(brbTree(brdTreeID).antRefval,brbTree(brdTreeID).conRefval);
        rulebase=struct;
        d=prod(size(brbTree(brdTreeID).antRefval));
        for i=1:size(rule,1)
            rulebase(i).conse=rule(i,size(brbTree(brdTreeID).antRefval,1)+1:end);
            rulebase(i).ruleweight=1;
        end
        %rulebase.
        size(brbTree(brdTreeID).antecedent,1);
        %observedOutput_old=mapObj(brbTree(brdTreeID).consequent{1});
        observedOutput=cell2mat(valueSet(find(strcmp(keySet,brbTree(brdTreeID).consequent{1}))));
        transformedRefVal={};
        fprintf(fid_x1,'\nAntecedents:');
        fprintf('\nAntecedents:');
        for antecedentID=1:size(brbTree(brdTreeID).antecedent,1)
            fprintf(fid_x1,' %s(',brbTree(brdTreeID).antecedent{antecedentID});
            fprintf(' %s(',brbTree(brdTreeID).antecedent{antecedentID});
            % in_old=mapObj(brbTree(brdTreeID).antecedent{antecedentID,1});
            in=cell2mat(valueSet(find(strcmp(keySet,brbTree(brdTreeID).antecedent{antecedentID,1}))));
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
        fprintf (fid_x1,'%2.2f ',z);
        fprintf (fid_x1,'\n');
        %     j=1;
        %     for i=1:length(z)/length(conseQuentRef)
        %         fprintf (fid_x1,'%2.2f ', z(j:j+length(conseQuentRef)-1));
        %         fprintf (fid_x1,'\n');
        %         j=j+length(conseQuentRef);
        %     end
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
        %fprintf(fid_x1,'\nAeq\n');
        %fprintf('%d ',Aeq)
        %fprintf(fid_x1,[repmat('%2.2f\t', 1, size(Aeq, 2)) '\n'], Aeq');
        %   Aeq= [0     0     0     0     1     1     1     0     0     0     0     0     0;
        %         0     0     0     0     0     0     0     1     1     1     0     0     0;
        %         0     0     0     0     0     0     0     0     0     0     1     1     1]
        beq = ones(numOfRuleWeight,1);
        
        A=zeros(numOfconRefval-1,numOfVariables);
        A(numOfVariables-numOfconRefval:numOfVariables)
        A(1,numOfVariables-numOfconRefval+1:numOfVariables)=[1 -1 0 ];
        A(2,numOfVariables-numOfconRefval+1:numOfVariables)=[0 1 -1 ];
        
        
        %     A(1,numOfVariables-numOfconRefval+1:numOfVariables)=[1 -1 0 0 0];
        %     A(2,numOfVariables-numOfconRefval+1:numOfVariables)=[0 1 -1 0 0];
        %     A(3,numOfVariables-numOfconRefval+1:numOfVariables)=[0 0  1 -1 0];
        %     A(4,numOfVariables-numOfconRefval+1:numOfVariables)=[0 0  0  1 -1];
        %     A(5,numOfVariables-numOfconRefval+1:numOfVariables)=[0 0  1 -1 0];
        B=zeros(numOfconRefval-1,1);
        %     save('dumpGlobalVariable.mat','input', 'outputOpti', 'observedOutput',...
        %         'transformedRefVal', 'conseQuentRef', 'rulebase', 'sizeOfData',...
        %         'numOfVariables', 'numOfconRefval', 'numOfAttrWeight', 'numOfRuleWeight', 'numOfbeliefDegrees','Aeq','beq','d');
        %    brbConfigdata.currentbrbTree=currentbrbTree;
        numOfAntecedentsRefVals=0;
        brbConfigdata.conseQuentRef=conseQuentRef;
        brbConfigdata.numOfAttrWeight=numOfAttrWeight;
        brbConfigdata.numOfconRefval=numOfconRefval;
        brbConfigdata.input=input;
        brbConfigdata.numOfAntecedentsRefVals=numOfAntecedentsRefVals;
        brbConfigdata.outputOpti=outputOpti;
        brbConfigdata.observedOutput=observedOutput;
        brbConfigdata.transformedRefVal=transformedRefVal;
        brbConfigdata.rulebase=rulebase;
        brbConfigdata.sizeOfData=sizeOfData;
        brbConfigdata.numOfVariables=numOfVariables;
        brbConfigdata.numOfRuleWeight=numOfRuleWeight;
        brbConfigdata.numOfbeliefDegrees=numOfbeliefDegrees;
        
        display('starting Optimization------------------------------------------')
        
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
        NP = 10*numOfVariables;
        %NP=5;
        if (NP>500)
            NP=500;
        end
        %NP = 100;
        % itermax       maximum number of iterations (generations)
        itermax = 5000;
        
        % F             DE-stepsize F ex [0, 2]
        F = 0.8;
        
        % CR            crossover probabililty constant ex [0, 1]
        CR = 0.5;
        
        % strategy       1 --> DE/best/1/exp           6 --> DE/best/1/bin
        %                2 --> DE/rand/1/exp           7 --> DE/rand/1/bin
        %                3 --> DE/rand-to-best/1/exp   8 --> DE/rand-to-best/1/bin
        %                4 --> DE/best/2/exp           9 --> DE/best/2/bin
        %                5 --> DE/rand/2/exp           else  DE/rand/2/bin
        
        strategy = 1;
        
        % refresh       intermediate output will be produced after "refresh"
        %               iterations. No intermediate output will be produced
        %               if refresh is < 1
        refresh =1000;
        tic
        trainparameter=zeros(5,numOfVariables);
        fvalue=zeros(1,1);
        nofof=zeros(1,1);
        % s = struct([])
       for i=1:1
            v(i)= rng(i);
            [x,f,nf] = BRBaDE('objFunAllParallelv3',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh,brbConfigdata,i)
            trainparameter(i,:)=x';
            fvalue(i)=f;
            nofof(i)=nf;
        end
        for i=1:1
            filename = strcat('randomConf_',num2str(i),'.mat');
            ss=v(i);
            save(filename,'ss');
            fprintf(fid_x1,'\nF=%2.5f\n',fvalue(i));
            fprintf (fid_x1,'Number of Function call= %2.2f\n',nofof(i));
        end
        bestFval=find(fvalue==min(fvalue));
        x=trainparameter(bestFval,:);
        fprintf(fid_x1,'\nbest seed=%2.0f\n',bestFval*100);
        fprintf(fid_x1,'best F=%2.5f\n',min(fvalue));
        fprintf(fid_x1,'X=>');
        fprintf(fid_x1,'%2.4f ',x);
        t=toc;
        
        % t_data_rule=vertcat(t_data_rule,data_rule);
        x0=x;
        fprintf(fid_tp,'%s=',brbTree(brdTreeID).consequent{1});
        fprintf(fid_tp,'%2.5f,',x0);
        fprintf(fid_tp,'\n');
        
        mapObj(brbTree(brdTreeID).consequent{1})=outputOpti;
        fprintf(fid_x1,'\nOptimizied value\n');
        fprintf (fid_x1,'Attribute Weights\n');
        fprintf (fid_x1,'%2.2f ', x0(1:numOfAttrWeight) );
        fprintf (fid_x1,'\nRuleWeights\n');
        fprintf (fid_x1,'%2.2f ', x0(numOfAttrWeight+1:numOfAttrWeight+numOfRuleWeight) );
        fprintf (fid_x1,'\nBelief Degrees\n');
        z=x0(numOfAttrWeight+numOfRuleWeight+1:numOfVariables);
        fprintf (fid_x1,'%2.2f ',z);
        fprintf (fid_x1,'\n');
        
       
        
        valueSet=allvalueSet(test,:);
        sizeOfData=size(valueSet,1);
        numberOfInputData=sizeOfData;
        keySet=cellstr(keySet);
        valueSet=num2cell(valueSet,1);
        mapObj = containers.Map(keySet,valueSet);
        transformedRefVal={};
        for antecedentID=1:size(brbTree(brdTreeID).antecedent,1)
            in=mapObj(brbTree(brdTreeID).antecedent{antecedentID,1});
            antcedentRefVal=brbTree(brdTreeID).antRefval(antecedentID,:);
            tmp=inputTransform(in,antcedentRefVal,numberOfInputData);
            transformedRefVal(antecedentID,:)={tmp};
            %      transformedRefVal(:,antecedentID)=tmp
            attrWeight(antecedentID)=1;
        end
         observedOutput=cell2mat(valueSet(find(strcmp(keySet,brbTree(brdTreeID).consequent{1}))));
         brbConfigdata.observedOutput=observedOutput;
         brbConfigdata.transformedRefVal=transformedRefVal;
         brbConfigdata.sizeOfData=sizeOfData;
        [f,outputOpti]=objFunAllevav2(x,brbConfigdata);
        fprintf (fid_x1,'Output vals:');
        fprintf (fid_x1,'%2.5f ', outputOpti );
        fprintf (fid_x1,'\nF vals:');
        fprintf (fid_x1,'%2.5f ', f );
        fprintf(fid_x1,'\nElapsed time is %2.4f sec\n',t);
    end
end
dateString = datestr(datetime('now'));
fprintf(fid_x1,'Ending program %s\n',dateString);
fprintf('Ending program %s\n',dateString);
% fprintf ( fid_x1,'____________________________\n');
% fprintf ( fid_x1,[repmat('%2.2f\t', 1, size(t_data_rule, 2)) '\n'], t_data_rule );
fclose(fid_x1);
fclose(fid_f1);
fclose(fid_tp);
