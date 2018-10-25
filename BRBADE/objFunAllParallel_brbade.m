function f=objFunAllParallel_brbade(x1,brbdeConfigdata)
%fprintf('o');
% global  input outputOpti observedOutput conseQuentRef ...
%     transformedRefVal ...
%     noOfRules rulebase sizeOfData...
%     numOfVariables numOfconRefval numOfAttrWeight numOfRuleWeight numOfbeliefDegrees ...
%     fid_x1 fid_f1;

% load('dumpGlobalVariable.mat','input', 'outputOpti', 'observedOutput',...
%     'transformedRefVal', 'conseQuentRef', 'rulebase', 'sizeOfData',...
%     'numOfVariables', 'numOfconRefval', 'numOfAttrWeight', 'numOfRuleWeight', 'numOfbeliefDegrees');

% load('BRBADE_dumpGlobalVariable.mat',...
%     'fid_x3','transformedRefVal', 'conseQuentRef', 'rulebase', 'numberOfInputData',...
%     'numOfVariables', 'numOfconRefval', 'numOfAttrWeight', 'numOfRuleWeight', 'numOfbeliefDegrees');
transformedRefVal=brbdeConfigdata.transformedRefVal;
conseQuentRef=brbdeConfigdata.conseQuentRef;
rulebase=brbdeConfigdata.rulebase;
numberOfInputData=brbdeConfigdata.numberOfInputData;
numOfVariables=brbdeConfigdata.numOfVariables;
numOfconRefval=brbdeConfigdata.numOfconRefval;
numOfAttrWeight=brbdeConfigdata.numOfAttrWeight;
numOfRuleWeight=brbdeConfigdata.numOfRuleWeight;
numOfbeliefDegrees=brbdeConfigdata.numOfbeliefDegrees;
sizeOfData=numberOfInputData;

% formatOut = 'yyyy-mmm-dd_HH_MM_SS';
% dateString = datestr(datetime('now'),formatOut);
% s = strcat('Log/crisp_',dateString,'.txt');
% fid_nonC1=fopen(s,'w');
%fid_nonC1=fid_x3;
%fprintf ( fid_x1,'%f ', x1 );
%fprintf ( fid_x1,'\n');
attrWeight=x1(1:numOfAttrWeight);


for trsId=1:size(transformedRefVal,1)
    transformedRefValM(:,:,trsId)=cell2mat(transformedRefVal(trsId));
end
crispValue=zeros(sizeOfData,1);
for data_id=1:sizeOfData
    %    size(transformedRefVal,1);
    %    size(transformedRefVal,2);
    %fprintf ( fid_nonC1,'Transformed input\n');
    for trsId=1:size(transformedRefVal,1)
        transformedRefVal(trsId);
     %   fprintf ( fid_nonC1,'%f,',transformedRefValM(data_id  ,:,trsId));
     %   fprintf(fid_nonC1,'\n');
    end
    
    matchingDegree=calMatchingDegree(transformedRefValM(data_id  ,:,1:size(transformedRefValM,3)),attrWeight);
    
    % Assigning RuleWeights from main data
    for i=1:numOfRuleWeight
        ruleweight(i,:)=rulebase(i).ruleweight;
    end
    % Assigning RuleWeights from fmincon x1
    j=numOfAttrWeight ;
    for i=1:numOfRuleWeight
        j=j+1;
        ruleweight(i)=x1(j);
        
    end
    z=sum(ruleweight.*matchingDegree);
    activationWeight=(ruleweight.*matchingDegree)./z;
    for i=1:numOfRuleWeight
        rulebase(i).activationWeight=activationWeight(i,1);
    end
    %getRulebase -- belief degree of consequents
    beta=zeros(numOfRuleWeight:size(rulebase(1).conse,2));
    % Assigning Belief Degrees from main data
    for i=1:numOfRuleWeight
        beta(i,:)=rulebase(i).conse;
    end
    
    j=numOfAttrWeight +numOfRuleWeight+1;
    % Assigning Belief Degrees from fmincon x1
    for i=1:numOfRuleWeight
        %j=j+1
        beta(i,:)=x1(j:j+numOfconRefval-1);
        j=j+numOfconRefval;
    end
    
    beta;
    %findMN
    MN=beta.*activationWeight;
    %MN = transpose(MN)
    %findMD
    total=sum(beta,2);
    MD=1-activationWeight.*total;
    rowsum=prod(MN+MD);
    rowsum_total=sum(rowsum);
    mh=prod(MD);
    kn=rowsum_total-(2*mh);
    kn1=1/kn;
    m=kn1*(rowsum-mh);
    mhn=kn1*mh;
    aggregatedValues=m/(1-mhn);
%     fprintf ( fid_nonC1,'matchingDegree[');
%     fprintf ( fid_nonC1,'%2.2f ', matchingDegree );
%     fprintf ( fid_nonC1,']\n');
%     fprintf ( fid_nonC1,'activationWeight[');
%     fprintf ( fid_nonC1,'%2.2f ', activationWeight );
%     fprintf ( fid_nonC1,']\n');
%     fprintf ( fid_nonC1,'=========================');
%     fprintf ( fid_nonC1,'\n');
%     
%     fprintf ( fid_nonC1,'beta matchingDegree activationWeight\n');
%     fprintf ( fid_nonC1,[repmat('%2.2f\t', 1, size(beta, 2)) '%2.2f\t%2.2f\n'],[ beta, matchingDegree, activationWeight ]');
%     fprintf ( fid_nonC1,'\n');
%     fprintf ( fid_nonC1,'aggregatedValues[');
%     fprintf ( fid_nonC1,'%2.2f ', aggregatedValues );
%     fprintf ( fid_nonC1,']\n');
    
    
    crispValue(data_id)=sum(aggregatedValues.*conseQuentRef,2);
    if isnan(crispValue(data_id))
        crispValue(data_id)=0;
    end
    f=crispValue(data_id);
end
% fprintf ( fid_nonC1,'____________________________\n');
% fprintf ( fid_nonC1,'Crisp value=>');
% fprintf ( fid_nonC1,'%f ', crispValue );
% fprintf ( fid_nonC1,'\n');
%fclose(fid_nonC1);
return
end
