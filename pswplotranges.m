function stop = pswplotranges(optimValues,state)
% formatOut = 'yyyy-mmm-dd';
% dateString = datestr(datetime('now'),formatOut);
% s = strcat('Log/BRBinsDE_',dateString,'.txt');
% fid_x1 = fopen (s, 'a');
stop = false; % This function does not stop the solver
ss=sprintf('pso_data.mat');
switch state
    case 'init'
        data(1,1)=0;
        data(1,2)=0;
        data(1,3)=0;
        data(1,4)=0;
        data(1,5)=0;
        data(1,6)=0;
        data(1,7)=0;
        data(1,8)=0;
        data(1,9)=0;
        data(1,10)=1;
        
        delete(ss)
        save(ss,'data')
    case 'iter'
        load(ss,'data');
        %fprintf(fid_x1,"\nPC:%4.2f,FC:%4.2f",optimValues.iteration,optimValues.bestfval);
        data(optimValues.iteration,9)=optimValues.bestfval;
        data(optimValues.iteration,10)=optimValues.iteration;
        save(ss,'data');
    case 'done'
        % No cleanup necessary
%  fclose(fid_x1);
end
