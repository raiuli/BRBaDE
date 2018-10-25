function f=applyPenaltyv2(x,numOfAttrWeight,numOfRuleWeight,numOfconRefval,numOfVariables,conseQuentRef,numOfAntecedentsRefVals,brbConfigdata,XVmin,XVmax)
        lam=0.5;
        constrainUnstaisfied=0;
        t=x;
        %fprintf ('Attribute Weights\n');
        %fprintf ('%2.2f ', t(1:numOfAttrWeight) );
        c=t(1:numOfAttrWeight);
        if sum(c< 0) >0
            constrainUnstaisfied=constrainUnstaisfied+lam*1^2;
        end
        if sum(c> 1) >0
            constrainUnstaisfied=constrainUnstaisfied+lam*1^2;
        end
        %fprintf ('\nRuleWeights\n');
        %fprintf ('%2.2f ', t(numOfAttrWeight+1:numOfAttrWeight+numOfRuleWeight) );
        c=t(numOfAttrWeight+1:numOfAttrWeight+numOfRuleWeight);
        if sum(c<0) >0
            constrainUnstaisfied=constrainUnstaisfied+lam*1^2;
        end
        if sum(c>1) >0
            constrainUnstaisfied=constrainUnstaisfied+lam*1^2;
        end
        %fprintf ('\nBelief Degrees\n');
        z=t(numOfAttrWeight+numOfRuleWeight+1:numOfVariables-numOfconRefval-numOfAntecedentsRefVals);
        j=1;
        for l=1:length(z)/length(conseQuentRef)
            %fprintf ('%2.2f ', z(j:j+length(conseQuentRef)-1));
            c=z(j:j+length(conseQuentRef)-1);
            if sum(c<0) >0
                constrainUnstaisfied=constrainUnstaisfied+lam*1^2;
            end
            if sum(c>1) >0
                constrainUnstaisfied=constrainUnstaisfied+lam*1^2;
            end
             if sum(c)~=1
                 constrainUnstaisfied=constrainUnstaisfied+lam*sum(c)^2;
             end    
            %fprintf ('\n');
            j=j+length(conseQuentRef);
        end
        consRefVal=t(numOfVariables-numOfconRefval-numOfAntecedentsRefVals+1:numOfVariables-numOfAntecedentsRefVals);
        a=sort(consRefVal,'descend');
        
        if ~isequal(a,consRefVal)
               constrainUnstaisfied=constrainUnstaisfied+lam*1^2;
        end
        %antecedent attribute utilty values
        antcedentsUVs=t(numOfVariables-numOfAntecedentsRefVals+1:numOfVariables);
        antcedentsUVsXVmin=XVmin(numOfVariables-numOfAntecedentsRefVals+1:numOfVariables);
        antcedentsUVsXVmax=XVmax(numOfVariables-numOfAntecedentsRefVals+1:numOfVariables);
        marker=1;
        for i=1:length(brbConfigdata.brbTree.antRefval)
            antecedentUVSize=size(brbConfigdata.brbTree.antRefval{i},2);
            antcedentUVs=antcedentsUVs(marker:marker+antecedentUVSize-1);
            antcedentUVsXVmin=antcedentsUVsXVmin(marker:marker+antecedentUVSize-1);
            antcedentUVsXVmax=antcedentsUVsXVmax(marker:marker+antecedentUVSize-1);
            
            if antcedentUVs(1)~=antcedentUVsXVmax(1)
                constrainUnstaisfied=constrainUnstaisfied+lam*1^2;
            end    
            if antcedentUVs(size(antcedentUVs,2))~=antcedentUVsXVmin(1)
                constrainUnstaisfied=constrainUnstaisfied+lam*1^2;
            end 
            a=sort(antcedentUVs,'descend');
        
            if ~isequal(a,antcedentUVs)
                   constrainUnstaisfied=constrainUnstaisfied+lam*1^2;
            end
            marker=marker+antecedentUVSize;
        end

        f=   constrainUnstaisfied;
        %fprintf ('Penalty=>%2.2f\n', f);
        %t(numOfVariables-numOfconRefval+1:numOfVariables) = sort(t(numOfVariables-numOfconRefval+1:numOfVariables),'descend');
end