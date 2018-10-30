clear;
clc;
format compact;
fclose('all');

REPS = 1000; minTime = Inf; nsum = 10;
      tic
      for i=1:REPS
        tstart = tic;
        sum = 0; 
        for j=1:nsum, sum = sum + besselj(j,REPS); end
        telapsed = toc(tstart)
        minTime = min(telapsed,minTime)
      end
      averageTime = toc/REPS
      toc
      tic
      for i=1:REPS
        tstart = tic;
        sum = 0; 
        for j=1:nsum, sum = sum + besselj(j,REPS); end
        telapsed = toc(tstart)
        minTime = min(telapsed,minTime)
      end
      averageTime = toc/REPS
      aa=toc
      aa