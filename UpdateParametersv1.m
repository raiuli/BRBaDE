function [lb,ub,brbconfigdata]=UpdateParametersv1(brbconfigdata,fid_x1,lbCU,ubCU)

brbTree=brbconfigdata.brbTree;
conseQuentRef=brbTree.conRefval;
in=brbconfigdata.input;
numOfAttrWeight=size(brbTree.antRefval,1);
numOfconRefval=size(brbTree.conRefval,2);
%read initial rule base for subrule base 1

for i=1:size(brbTree.antecedent,1)
    if i==1
        xb3=combvec(cell2mat(brbTree.antRefval(i)));
    else
        xb3=combvec(cell2mat(brbTree.antRefval(i)),xb3);
    end
end

%         rule=calculateInitialRulebasev2(cell2mat(brbTree.antRefval),brbTree.conRefval);
%         rulebase=struct;
%         for i=1:size(rule,1)
%            rulebase(i).conse=rule(i,size(brbTree(brdTreeID).antRefval,1)+1:end);
%            rulebase(i).ruleweight=1;
%         end
%rulebase.

%observedOutput_old=mapObj(brbTree(brdTreeID).consequent{1});
%observedOutput=cell2mat(valueSet(find(strcmp(keySet,brbTree(brdTreeID).consequent{1}))));
%lbCU=horzcat(lbCU,ones(1,length(conseQuentRef))*min(observedOutput));
%ubCU=horzcat(ubCU,ones(1,length(conseQuentRef))*max(observedOutput));
transformedRefVal={};
lbAU=[];
ubAU=[];
numOfAntecedentsRefVals=0;
fprintf(fid_x1,'\nAntecedents:');
fprintf('\nAntecedents:');
for antecedentID=1:size(brbTree.antecedent,1)
    fprintf(fid_x1,' %s(',brbTree.antecedent{antecedentID});
    fprintf(' %s(',brbTree.antecedent{antecedentID});
    antcedentRefVal=cell2mat(brbTree.antRefval(antecedentID,:));
    lbAU=horzcat(lbAU,ones(1,length(antcedentRefVal))*min(in(antecedentID,:)));
    ubAU=horzcat(ubAU,ones(1,length(antcedentRefVal))*max(in(antecedentID,:)));
    numOfAntecedentsRefVals=numOfAntecedentsRefVals+length(antcedentRefVal);
    fprintf(fid_x1,'%2.2f ',antcedentRefVal);
    fprintf(fid_x1,')');
    fprintf('%2.2f ',antcedentRefVal);
    fprintf(')');
    attrWeight(antecedentID)=1;
end
fprintf(fid_x1,'=>%s (',brbTree.consequent{1});
fprintf(fid_x1,'%2.2f ',brbTree.conRefval);
fprintf(fid_x1,')\n');

fprintf('=>%s (',brbTree.consequent{1});
fprintf('%2.2f ',brbTree.conRefval);
fprintf(')\n');

numOfRuleWeight=size(xb3',1);
numOfbeliefDegrees=numOfRuleWeight*numOfconRefval;
%numOfVariables=numOfAttrWeight+numOfRuleWeight+numOfbeliefDegrees;
numOfVariables=numOfconRefval+numOfAttrWeight+numOfRuleWeight+numOfbeliefDegrees+numOfAntecedentsRefVals;
fprintf(fid_x1,'Number of Varaibles: %d=%d(CR)+%d(AW)+%d(RW)+%d(BD)+%d(Arefv)\n',numOfVariables,numOfconRefval,numOfAttrWeight,numOfRuleWeight,numOfbeliefDegrees,numOfAntecedentsRefVals);
fprintf('Number of Varaibles: %d=%d(CR)+%d(AW)+%d(RW)+%d(BD)+%d(Arefv)\n',numOfVariables,numOfconRefval,numOfAttrWeight,numOfRuleWeight,numOfbeliefDegrees,numOfAntecedentsRefVals);

%initialiaze the constraints
lb = zeros(1,numOfVariables);
ub =ones(1,numOfVariables);
lb(numOfAttrWeight+numOfRuleWeight+numOfbeliefDegrees+1:numOfVariables-numOfAntecedentsRefVals)=lbCU;
ub(numOfAttrWeight+numOfRuleWeight+numOfbeliefDegrees+1:numOfVariables-numOfAntecedentsRefVals)=ubCU;
lb(numOfAttrWeight+numOfRuleWeight+numOfbeliefDegrees+numOfconRefval+1:numOfVariables)=lbAU;
ub(numOfAttrWeight+numOfRuleWeight+numOfbeliefDegrees+numOfconRefval+1:numOfVariables)=ubAU;

%brbConfigdata.conseQuentRef=conseQuentRef;
brbconfigdata.numOfAttrWeight=numOfAttrWeight;
%brbConfigdata.numOfconRefval=numOfconRefval;
%brbConfigdata.input=in;
brbconfigdata.numOfAntecedentsRefVals=numOfAntecedentsRefVals;
%brbConfigdata.outputOpti=outputOpti;
%brbConfigdata.observedOutput=observedOutput;
%brbConfigdata.transformedRefVal=transformedRefVal;
%brbConfigdata.rulebase=rulebase;
%brbConfigdata.sizeOfData=sizeOfData;
brbconfigdata.numOfVariables=numOfVariables;
brbconfigdata.numOfRuleWeight=numOfRuleWeight;
brbconfigdata.numOfbeliefDegrees=numOfbeliefDegrees;
%brbConfigdata.brbTree=brbTree;

end