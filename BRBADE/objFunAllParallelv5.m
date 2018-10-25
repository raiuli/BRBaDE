function f=objFunAllParallelv5(x1,brbConfigdata)

conseQuentRef=brbConfigdata.conseQuentRef;
numOfAttrWeight=    brbConfigdata.numOfAttrWeight;
numOfconRefval=    brbConfigdata.numOfconRefval;
%input= brbConfigdata.input;
%numOfAntecedentsRefVals=    brbConfigdata.numOfAntecedentsRefVals;
%outputOpti=  brbConfigdata.outputOpti;
observedOutput=    brbConfigdata.observedOutput;
transformedRefVal=   brbConfigdata.transformedRefVal;
rulebase=  brbConfigdata.rulebase;
sizeOfData=   brbConfigdata.sizeOfData;
numOfVariables=   brbConfigdata.numOfVariables;
numOfRuleWeight=   brbConfigdata.numOfRuleWeight;
numOfbeliefDegrees=    brbConfigdata.numOfbeliefDegrees;
formatOut = 'yyyy-mmm-dd_HH_MM_SS';
%dateString = datestr(datetime('now'),formatOut);
%s = strcat('Log/crisp.txt_',dateString,'.txt');
%fid_nonC1=fopen('Log/crisp.txt','w');
%fprintf ( fid_x1,'%f ', x1 );
%fprintf ( fid_x1,'\n');
attrWeight=x1(1:numOfAttrWeight);
for trsId=1:size(transformedRefVal,1)
    transformedRefValM(:,:,trsId)=cell2mat(transformedRefVal(trsId));
end 
crispValue=zeros(sizeOfData,1);
for data_id=1:sizeOfData
    matchingDegree=calMatchingDegree(transformedRefValM(data_id  ,:,1:size(transformedRefValM,3)),attrWeight);
     
    % Assigning RuleWeights from fmincon x1
    ruleweight=ones(numOfRuleWeight,1);
    ruleweight=x1(numOfAttrWeight+1:numOfAttrWeight+numOfRuleWeight)';
    z=sum(ruleweight.*matchingDegree);
    activationWeight=(ruleweight.*matchingDegree)./z;
    for i=1:numOfRuleWeight
        rulebase(i).activationWeight=activationWeight(i,1);
    end
%getRulebase -- belief degree of consequents
    % Assigning Belief Degrees from fmincon x1
     beta=zeros(numOfRuleWeight:size(rulebase(1).conse,2));
     j=numOfAttrWeight +numOfRuleWeight+1;
     a=x1(numOfAttrWeight +numOfRuleWeight+1:numOfVariables-numOfconRefval);
     beta=reshape(a,[numOfconRefval,numOfRuleWeight])';
     
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
   conseQuentRef=x1(numOfVariables-numOfconRefval+1:numOfVariables);
   % data_id;
    crispValue(data_id)=sum(aggregatedValues.*conseQuentRef,2);
   if isnan(crispValue(data_id))
        crispValue(data_id)=0;
    end
end
 %   fprintf ( fid_nonC1,'____________________________\n');
 %   fprintf ( fid_nonC1,'x=>');
 %   fprintf ( fid_nonC1,'%f ', x1 );
 %   fprintf ( fid_nonC1,'\n');
 %   fprintf ( fid_nonC1,'Crisp value=>');
 %   fprintf ( fid_nonC1,'%f ', crispValue );
 %   fprintf ( fid_nonC1,'\n');
 %   fprintf ( fid_nonC1,'observedOutput=>');
 %   fprintf ( fid_nonC1,'%f ', observedOutput );
 %   fprintf ( fid_nonC1,'\n');
    %fclose(fid_nonC1);
%      fid_crisp1 = fopen ('crisp1.txt', 'a');
%      fprintf ( fid_crisp1,'____________________________\n');
%     for i=1:size(crispValue,2)
%      fprintf ( fid_crisp1,'%f ', crispValue(i) );
%      fprintf ( fid_crisp1,'\n');
%     end
%     fclose(fid_crisp1);
% outputOpti=crispValue;
%  fprintf('%2.2f ',outputOpti);
%  fprintf('\n');
% f_v=zeros(data_id,1);
% for data_id=1:sizeOfData
%     f_v(data_id)=sum((crispValue(data_id)-observedOutput(data_id))^2);
% end
f_v=sum((crispValue(:)-observedOutput(:)).^2);
% sum_f=sum(f_v,2);
f=f_v/sizeOfData;
%size(crispValue);
%size(observedOutput);
%f_sqrt=sqrt(sum( (crispValue(:)-observedOutput(:)).^2))/numel(crispValue);
if isnan(f)==1
    fprintf('non nan');
%    size(x1);
%    size(crispValue);
%    size(observedOutput);
end
%fprintf ( fid_f1,'%f ', f );
%fprintf ( fid_f1,'\n');
%figure(2)
%plot3(bestmem(1,1),bestmem(1,2),0,'dk','markersize',10);
%grid on
%    fprintf ( fid_nonC1,'f_sqrt= %f ', f_sqrt );
%    fprintf ( fid_nonC1,'\n');
%    fprintf ( fid_nonC1,'f= %f ', f );
%    fprintf ( fid_nonC1,'\n');
   
    outputOpti=crispValue;
%save('dumpGlobalVariable.mat','outputOpti');  
%fprintf(1,'i:%d output f: %f:\n',yyy,f);
%fprintf(fid_nonC1,'i:%d output f: %f:\n',yyy,f);
%fclose(fid_nonC1);
% toc
return
end

% T=length(crispValue);
% MAE=sum(abs(crispValue(:)-observedOutput(:)))/T;
% r=0.1;
% part1=(d*log(2)-log(r));
% part2=sqrt((part1*(conseQuentRef(3)-conseQuentRef(1))^2)/(2*T));
% ge=MAE+part2;