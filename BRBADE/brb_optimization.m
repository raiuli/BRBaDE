function [F, CR]=brb_optimization(d11,d12, d21,d22)
    formatOut = 'yyyy-mmm-dd_HH_MM_SS.FFF';
    dateString = datestr(datetime('now'),formatOut);
     %s = strcat('Log/Opti/BRBinsOpti_',dateString,'.txt');
     s = strcat('Log/Opti/BRBinsOpti_.txt');
    fid_x3 = fopen (s, 'w');
     
   % logger = logging.getLogger('log_brb_optimization', 'path', 'logger2.log');
    F=0;
    CR=0;

       
%        brbTree(1).antecedent=cellstr(['d11';'d12']);
%        brbTree(1).antRefval=[1 0.5 0;
%                       1 0.5 0;];
%        brbTree(1).consequent=cellstr('F');
%        brbTree(1).conRefval=[2 1 0.1];
%        
%        brbTree(2).antecedent=cellstr(['d21';'d22']);
%        brbTree(2).antRefval=[2 1 0;
%                       2 1 0;];
%        brbTree(2).consequent=cellstr('CR');
%        brbTree(2).conRefval=[ 1 0.5 0.1];
%        brbTree(1).antecedent=cellstr(['d11';'d12']);
%        brbTree(1).antRefval=[1 0.75 0.5 0;
%                       1  0.75 0.5 0;];
%        brbTree(1).consequent=cellstr('F');
%        brbTree(1).conRefval=[2 1.5 1 0.1];
%        
%        brbTree(2).antecedent=cellstr(['d21';'d22']);
%        brbTree(2).antRefval=[2 1.5 1 0;
%                       2 1.5 1 0;];
%        brbTree(2).consequent=cellstr('CR');
%        brbTree(2).conRefval=[ 1 0.75 0.5 0.1];

%        brbTree(1).antecedent=cellstr(['d11';'d12']);
%        brbTree(1).antRefval=[1 0.75 0.5 0;
%                       1  0.75 0.5 0;];
%        brbTree(1).consequent=cellstr('F');
%        brbTree(1).conRefval=[2 1.5 1 0.1];
%        
%        brbTree(2).antecedent=cellstr(['d21';'d22']);
%        brbTree(2).antRefval=[2 1.5 1 0;
%                       2 1.5 1 0;];
%        brbTree(2).consequent=cellstr('CR');
%        brbTree(2).conRefval=[ 1 0.75 0.5 0.1];
      
       brbTree(1).antecedent=cellstr(['d11';'d12']);
       brbTree(1).antRefval=[1 0.8 0.6 0.2 0;
                      1 0.8 0.6 0.2 0;];
       %brbTree(1).consequent=cellstr('F');
       brbTree(1).consequent=cellstr('CR');
       %brbTree(1).conRefval=[2 1.6 1.2 0.8 0.4 0.1];
       brbTree(1).conRefval=[2 1.7 1.4 1.1 0.8 0.5];
      
       brbTree(2).antecedent=cellstr(['d21';'d22']);
       brbTree(2).antRefval=[2 1.6 1.2 0.8 0.4 0;
                     2 1.6 1.2 0.8 0.4 0;];
       %brbTree(2).consequent=cellstr('CR');
       brbTree(2).consequent=cellstr('F');
       brbTree(2).conRefval=[ 1 0.8 0.6 0.2 0.1];
       %brbTree(2).conRefval=[ 1 0.96 0.92 0.88 0.8];
       
       
       line="d11,d21,d12,d22,F,CR";
        fprintf(fid_x3,'line');
        fprintf(fid_x3,'%f,%f,%f,%f',d11,d12,d21,d22);
       keySet=split(line,',');
       keySet=cellstr(keySet);
       line=d11+","+d21+","+d12+","+d22+","+F+","+CR;
       numberOfInputData=1;
       valueSet(numberOfInputData,:)=str2num(line); 
       valueSetcell=num2cell(valueSet,1);
       mapObj = containers.Map(keySet,valueSetcell);
      
       for brdTreeID=1:size(brbTree,2)
           conseQuentRef=brbTree(brdTreeID).conRefval;
           numOfAttrWeight=size(brbTree(brdTreeID).antRefval,1);
           numOfconRefval=size(brbTree(brdTreeID).conRefval,2);
           %read initial rule base for subrule base 1
           
           rule=calculateInitialRulebase(brbTree(brdTreeID).antRefval,brbTree(brdTreeID).conRefval);
           rulebase=struct;
           for i=1:size(rule,1)
                rulebase(i).conse=rule(i,size(brbTree(brdTreeID).antRefval,1)+1:end);
                rulebase(i).ruleweight=1;
           end    
           %rulebase.
           size(brbTree(brdTreeID).antecedent,1);
           observedOutput=mapObj(brbTree(brdTreeID).consequent{1});
           transformedRefVal={};
    
           fprintf(fid_x3,'\nAntecedents:');
           %logger.info('Antecedents:')
           for antecedentID=1:size(brbTree(brdTreeID).antecedent,1)
        
                 fprintf(fid_x3,' %s(',brbTree(brdTreeID).antecedent{antecedentID});
        
        %         logger.info(sprintf(' %s(',brbTree(brdTreeID).antecedent{antecedentID}));
                    in=mapObj(brbTree(brdTreeID).antecedent{antecedentID,1});
                    antcedentRefVal=brbTree(brdTreeID).antRefval(antecedentID,:);

                    fprintf(fid_x3,'%2.2f ',antcedentRefVal);
                    fprintf(fid_x3,')');
                     fprintf(fid_x3,'[%f]',in);
                    %logger.info(sprintf('%2.2f ',antcedentRefVal)+")");
                    tmp=inputTransform(in,antcedentRefVal,numberOfInputData);
                    transformedRefVal(antecedentID,:)={tmp};
              %      transformedRefVal(:,antecedentID)=tmp
                    attrWeight(antecedentID)=1;
           end
           fprintf(fid_x3,'=>%s (',brbTree(brdTreeID).consequent{1});
           fprintf(fid_x3,'%2.2f ',brbTree(brdTreeID).conRefval);
           fprintf(fid_x3,')\n');
           numOfRuleWeight=size(rulebase,2);
           numOfbeliefDegrees=numOfRuleWeight*numOfconRefval;
           numOfVariables=numOfconRefval+numOfAttrWeight+numOfRuleWeight+numOfbeliefDegrees;
           fprintf(fid_x3,'Number of Varaibles: %d=%d(CR)+%d(AW)+%d(RW)+%d(BD)\n',numOfVariables,numOfconRefval,numOfAttrWeight,numOfRuleWeight,numOfbeliefDegrees);
                
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
            fprintf(fid_x3,'\nIntial value\n');
            fprintf (fid_x3,'Attribute Weights\n');
            fprintf (fid_x3,'%d ', x0(1:numOfAttrWeight) );
            fprintf (fid_x3,'\nRuleWeights\n');
            fprintf (fid_x3,'%d ', x0(numOfAttrWeight+1:numOfAttrWeight+numOfRuleWeight) );
            fprintf (fid_x3,'\nBelief Degrees\n');
           z=x0(numOfAttrWeight+numOfRuleWeight+1:numOfVariables-numOfconRefval);
           j=1;
           for i=1:length(z)/length(conseQuentRef)
                 fprintf (fid_x3,'%2.2f ', z(j:j+length(conseQuentRef)-1));
                 fprintf (fid_x3,'\n');
                 j=j+length(conseQuentRef);
           end    
          fprintf (fid_x3,'\nConsequent\n');
          fprintf (fid_x3,'%d \n', x0(numOfVariables-numOfconRefval+1:numOfVariables) ); 
           save('dumpGlobalVariable.mat','fid_x3','transformedRefVal', 'conseQuentRef', 'rulebase', 'numberOfInputData',...
                'numOfVariables', 'numOfconRefval', 'numOfAttrWeight', 'numOfRuleWeight', 'numOfbeliefDegrees');
%         if brdTreeID==2
%             F=objFunAllParallel(x0);
%             fprintf (fid_x3,'\nCR: %f \n',F);
%         elseif brdTreeID==1
%            CR=objFunAllParallel(x0); 
%            fprintf (fid_x3,'\nF: %f \n',CR);
%         end 

        if brdTreeID==2
            F=objFunAllParallel(x0);
            fprintf (fid_x3,'\nF: %f \n',F);
        elseif brdTreeID==1
           CR=objFunAllParallel(x0); 
           fprintf (fid_x3,'\nCR: %f \n',CR);
        end 

      
       end
fclose(fid_x3);