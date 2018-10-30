

% Stat clear previous data

clear;
clc;
format compact;
fclose('all');
filename='SecDataset';
ext='.txt';
SamplePoint=100;
ffilename = strcat(filename,ext);

fid = fopen (ffilename, 'r');

numberOfInputData=0;

while ~feof(fid)
    numberOfInputData=numberOfInputData+1;
    line=fgetl(fid);
    if numberOfInputData==1
        keySet=split(line,',');
    else
        allvalueSet(numberOfInputData-1,:)=str2num(line);
    end
end
fclose(fid);

y = datasample(allvalueSet,SamplePoint);

indices = crossvalind('Kfold',y(:,3),5);
ffilename = strcat(filename,'_indices_',ext);
csvwrite(ffilename,indices);

ffilename = strcat(filename,num2str(SamplePoint),ext);
fid = fopen (ffilename, 'w');
for i=1:length(keySet)
   if i< length(keySet)
 fprintf(fid,'%s,',keySet{i});
   else
       fprintf(fid,'%s',keySet{i});
   end   
end
fprintf(fid,'\n');
fprintf(fid,[repmat('%2.2f,', 1, size(y, 2)) '\n'], y');
fclose(fid);
