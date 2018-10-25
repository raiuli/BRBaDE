function f=applyPenalty(x,numOfAttrWeight,numOfRuleWeight,numOfconRefval,numOfVariables,conseQuentRef)
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
        z=t(numOfAttrWeight+numOfRuleWeight+1:numOfVariables-numOfconRefval);
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
%             if sum(c)~=1
%                 constrainUnstaisfied=constrainUnstaisfied+lam*sum(c)^2;
%             end    
            %fprintf ('\n');
            j=j+length(conseQuentRef);
        end
        consRefVal=t(numOfVariables-numOfconRefval+1:numOfVariables);
        a=sort(consRefVal,'descend');
        
        if ~isequal(a,consRefVal)
               constrainUnstaisfied=constrainUnstaisfied+lam*1^2;
        end    
        
        f=   constrainUnstaisfied;
        %fprintf ('Penalty=>%2.2f\n', f);
        %t(numOfVariables-numOfconRefval+1:numOfVariables) = sort(t(numOfVariables-numOfconRefval+1:numOfVariables),'descend');
end