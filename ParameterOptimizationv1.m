function [rulebase,f,nf]=UpdateParametersv1(brbconfigdata,fid_x1)
        
        brbTree=brbconfigdata.brbTree;
        conseQuentRef=brbTree.conRefval;
        
        
        numOfAttrWeight=size(brbTree.antRefval,1);
        numOfconRefval=size(brbTree.conRefval,2);
        %read initial rule base for subrule base 1
        
        rule=calculateInitialRulebase(cell2mat(brbTree(brdTreeID).antRefval),brbTree(brdTreeID).conRefval);
        rulebase=struct;
        for i=1:size(rule,1)
            rulebase(i).conse=rule(i,size(brbTree(brdTreeID).antRefval,1)+1:end);
            rulebase(i).ruleweight=1;
        end
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
            fprintf(fid_x1,' %s(',brbTree(brdTreeID).antecedent{antecedentID});
            fprintf(' %s(',brbTree(brdTreeID).antecedent{antecedentID});
            antcedentRefVal=cell2mat(brbTree(brdTreeID).antRefval(antecedentID,:));
            lbAU=horzcat(lbAU,ones(1,length(antcedentRefVal))*min(in(antecedentID,:)));
            ubAU=horzcat(ubAU,ones(1,length(antcedentRefVal))*max(in(antecedentID,:)));
            numOfAntecedentsRefVals=numOfAntecedentsRefVals+length(antcedentRefVal);
            fprintf(fid_x1,'%2.2f ',antcedentRefVal);
            fprintf(fid_x1,')');
            fprintf('%2.2f ',antcedentRefVal);
            fprintf(')');
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
        %initialVal{brdTreeID};
        fprintf(fid_x1,'\nIntial value\n');
        fprintf (fid_x1,'Attribute Weights\n');
        fprintf (fid_x1,'%d ', x0(1:numOfAttrWeight) );
        fprintf (fid_x1,'\nRuleWeights\n');
        fprintf (fid_x1,'%d ', x0(numOfAttrWeight+1:numOfAttrWeight+numOfRuleWeight) );
        fprintf (fid_x1,'\nBelief Degrees\n');
        z=x0(numOfAttrWeight+numOfRuleWeight+1:numOfVariables-numOfconRefval);
        fprintf (fid_x1,'%2.2f ',z);
        fprintf (fid_x1,'\nConsequent utlity values\n');
        fprintf (fid_x1,'%d ', x0(numOfVariables-numOfconRefval-numOfAntecedentsRefVals+1:numOfVariables-numOfAntecedentsRefVals) );
        fprintf (fid_x1,'\nAntecedent utlity values\n');
        fprintf (fid_x1,'%d ', x0(numOfVariables-numOfAntecedentsRefVals+1:numOfVariables) );

end