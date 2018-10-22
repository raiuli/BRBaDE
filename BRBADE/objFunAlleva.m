function f=objFunAlleva(x1)

global  input outputOpti observedOutput conseQuentRef ...
    transformedRefVal ...
    noOfRules rulebase sizeOfData...
    numOfVariables numOfconRefval numOfAttrWeight numOfRuleWeight numOfbeliefDegrees ...
    fid_x1 fid_f1 data_rule usedRules;

formatOut = 'yyyy-mmm-dd_HH_MM_SS';
dateString = datestr(datetime('now'),formatOut);
s = strcat('Log/data_rule_',dateString,'.txt');
fid_dr = fopen (s, 'w');


attrWeight=x1(1:numOfAttrWeight);

data_rule_o=zeros(numOfRuleWeight,sizeOfData);
countOfActivateRule_o=zeros(numOfRuleWeight,1);
data_rule=zeros(numOfRuleWeight,sizeOfData);
countOfActivateRule=zeros(numOfRuleWeight,1);
for trsId=1:size(transformedRefVal,1)
    transformedRefValM(:,:,trsId)=cell2mat(transformedRefVal(trsId));
end 
crispValue=zeros(sizeOfData,1);

for data_id=1:sizeOfData
    size(transformedRefVal,1);
    size(transformedRefVal,2);
    matchingDegree=calMatchingDegree(transformedRefValM(data_id  ,:,1:size(transformedRefValM,3)),attrWeight);
    %matchingDegree=matchingDegree(find(used))
    rulebase.ruleweight;
    for i=1:numOfRuleWeight
        ruleweight(i,:)=rulebase(i).ruleweight;
    end
    j=numOfAttrWeight ;
    for i=1:numOfRuleWeight
        j=j+1;
        ruleweight(i)=x1(j);
        
    end
    ruleweight.*matchingDegree;
    z=sum(ruleweight.*matchingDegree);
    activationWeight=(ruleweight.*matchingDegree)./z;
    for i=1:numOfRuleWeight
        rulebase(i).activationWeight=activationWeight(i,1);
    end
%getRulebase -- belief degree of consequents

     for i=1:numOfRuleWeight
         beta(i,:)=rulebase(i).conse;
     end
%     beta;
    j=numOfAttrWeight +numOfRuleWeight+1;
    for i=1:numOfRuleWeight
        %j=j+1
        beta(i,:)=x1(j:j+numOfconRefval-1);
        j=j+numOfconRefval;
    end
    beta;
    %findMN
    MN=beta.*activationWeight;
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
%     fprintf ( fid_x1,'matchingDegree[');
%     fprintf ( fid_x1,'%2.2f ', matchingDegree );
%     fprintf ( fid_x1,']\n');
%     fprintf ( fid_x1,'activationWeight[');
%    fprintf ( fid_x1,'%2.2f ', activationWeight );
%    fprintf ( fid_x1,']\n');
    
%     countOfActivateRule_o=countOfActivateRule_o+activationWeight;
%     tmp=find(activationWeight);
%     countOfActivateRule(tmp)=1;
%     data_rule_o(:,data_id)=activationWeight;
%     data_rule(tmp,data_id)=1;
%     fprintf ( fid_dr,'=========================\n');
%     fprintf ( fid_dr,'\n');
%     
%     fprintf ( fid_dr,'beta matchingDegree activationWeight\n');
%     fprintf ( fid_dr,[repmat('%2.2f\t', 1, size(beta, 2)) '%2.2f\t%2.2f\n'],[ beta, matchingDegree, activationWeight ]');
%      fprintf ( fid_dr,'\n');
%     fprintf ( fid_dr,'aggregatedValues[');
%     fprintf ( fid_dr,'%2.2f ', aggregatedValues );
%     fprintf ( fid_dr,']\n');

%    data_id;
    crispValue(data_id)=sum(aggregatedValues.*conseQuentRef,2);
    if isnan(crispValue(data_id))
        crispValue(data_id)=0;
    end
   
end
    fprintf ( fid_dr,'____________________________\n');
    fprintf ( fid_dr,'%f ', x1 );
    fprintf ( fid_dr,'\n');
    fprintf ( fid_dr,'____________________________\n');
    for i=1:size(crispValue,2)
     fprintf ( fid_dr,'%f ', crispValue(i) );
     fprintf ( fid_dr,'\n');
    end
%    fclose(fid_crisp1);
%outputOpti=crispValue;
%  fprintf('%2.2f ',outputOpti);
%  fprintf('\n');
% f_v=zeros(data_id,1);
% for data_id=1:sizeOfData
%     f_v(data_id)=sum((crispValue(data_id)-observedOutput(data_id))^2);
% end
% % sum_f=sum(f_v,2);
% % f=sum_f/sizeOfData;
% size(crispValue);
% size(observedOutput);
% f=sqrt(sum( (crispValue(:)-observedOutput(:)).^2))/numel(crispValue)
% if isnan(f)==1
%     fprintf('non nan');
%     size(x1)
%     size(crispValue)
%     size(observedOutput)
% end
%fprintf ( fid_f1,'%f ', f );
%fprintf ( fid_f1,'\n');
% countOfActivateRule;
% countOfActivateRule_o;
% data_rule;
% 
% data_rule_o;
return 
fclose(fid_dr);
end
