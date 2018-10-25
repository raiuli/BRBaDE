function [bestmem,bestval,nfeval] = evalobjFunAll(fname,VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh,brbConfigdata,instance_num)
%function [bestmem,bestval,nfeval] = deEva(fname,VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh);
% load('dumpGlobalVariable.mat', 'numOfAntecedentsRefVals','observedOutput',...
%     'transformedRefVal', 'conseQuentRef', 'rulebase', 'sizeOfData',...
%     'numOfVariables', 'numOfconRefval', 'numOfAttrWeight', 'numOfRuleWeight', 'numOfbeliefDegrees','Aeq','beq');

% minimization of a user-supplied function with respect to x(1:D),
% using the differential evolution (DE) algorithm of Rainer Storn
% (http://www.icsi.berkeley.edu/~storn/code.html)
%
% Special thanks go to Ken Price (kprice@solano.community.net) and
% Arnold Neumaier (http://solon.cma.univie.ac.at/~neum/) for their
% valuable contributions to improve the code.
%
% Strategies with exponential crossover, further input variable
% tests, and arbitrary function name implemented by Jim Van Zandt
% <jrv@vanzandt.mv.com>, 12/97.
%
% Output arguments:
% ----------------
% bestmem        parameter vector with best solution
% bestval        best objective function value
% nfeval         number of function evaluations
%
% Input arguments:
% ---------------
%
% fname          string naming a function f(x,y) to minimize
% VTR            "Value To Reach". devec3 will stop its minimization
%                if either the maximum number of iterations "itermax"
%                is reached or the best parameter vector "bestmem"
%                has found a value f(bestmem,y) <= VTR.
% D              number of parameters of the objective function
% XVmin          vector of lower bounds XVmin(1) ... XVmin(D)
%                of initial population
%                *** note: these are not bound constraints!! ***
% XVmax          vector of upper bounds XVmax(1) ... XVmax(D)
%                of initial population
% y		        problem data vector (must remain fixed during the
%                minimization)
% NP             number of population members
% itermax        maximum number of iterations (generations)
% F              DE-stepsize F from interval [0, 2]
% CR             crossover probability constant from interval [0, 1]
% strategy       1 --> DE/best/1/exp           6 --> DE/best/1/bin
%                2 --> DE/rand/1/exp           7 --> DE/rand/1/bin
%                3 --> DE/rand-to-best/1/exp   8 --> DE/rand-to-best/1/bin
%                4 --> DE/best/2/exp           9 --> DE/best/2/bin
%                5 --> DE/rand/2/exp           else  DE/rand/2/bin
%                Experiments suggest that /bin likes to have a slightly
%                larger CR than /exp.
% refresh        intermediate output will be produced after "refresh"
%                iterations. No intermediate output will be produced
%                if refresh is < 1
%
%       The first four arguments are essential (though they have
%       default values, too). In particular, the algorithm seems to
%       work well only if [XVmin,XVmax] covers the region where the
%       global minimum is expected. DE is also somewhat sensitive to
%       the choice of the stepsize F. A good initial guess is to
%       choose F from interval [0.5, 1], e.g. 0.8. CR, the crossover
%       probability constant from interval [0, 1] helps to maintain
%       the diversity of the population and is rather uncritical. The
%       number of population members NP is also not very critical. A
%       good initial guess is 10*D. Depending on the difficulty of the
%       problem NP can be lower than 10*D or must be higher than 10*D
%       to achieve convergence.
%       If the parameters are correlated, high values of CR work better.
%       The reverse is true for no correlation.
%
% default values in case of missing input arguments:
% 	VTR = 1.e-6;
% 	D = 2;
% 	XVmin = [-2 -2];
% 	XVmax = [2 2];
%	y=[];
% 	NP = 10*D;
% 	itermax = 200;
% 	F = 0.8;
% 	CR = 0.5;
% 	strategy = 7;
% 	refresh = 10;
%
% Cost function:  	function result = f(x,y);
%                      	has to be defined by the user and is minimized
%			w.r. to  x(1:D).
%
% Example to find the minimum of the Rosenbrock saddle:
% ----------------------------------------------------
% Define f.m as:
%                    function result = f(x,y);
%                    result = 100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%                    end
% Then type:
%
% 	VTR = 1.e-6;
% 	D = 2;
% 	XVmin = [-2 -2];
% 	XVmax = [2 2];
% 	[bestmem,bestval,nfeval] = devec3("f",VTR,D,XVmin,XVmax);
%
% The same example with a more complete argument list is handled in
% run1.m
%
% About devec3.m
% --------------
% Differential Evolution for MATLAB
% Copyright (C) 1996, 1997 R. Storn
% International Computer Science Institute (ICSI)
% 1947 Center Street, Suite 600
% Berkeley, CA 94704
% E-mail: storn@icsi.berkeley.edu
% WWW:    http://http.icsi.berkeley.edu/~storn
%
% devec is a vectorized variant of DE which, however, has a
% propertiy which differs from the original version of DE:
% 1) The random selection of vectors is performed by shuffling the
%    population array. Hence a certain vector can't be chosen twice
%    in the same term of the perturbation expression.
%
% Due to the vectorized expressions devec3 executes fairly fast
% in MATLAB's interpreter environment.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 1, or (at your option)
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. A copy of the GNU
% General Public License can be obtained from the
% Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

%-----Check input variables---------------------------------------------
err=[];
if nargin<1, error('devec3 1st argument must be function name'); else
    if exist(fname)<1; err(1,length(err)+1)=1; end; end;
if nargin<2, VTR = 1.e-6; else
    if length(VTR)~=1; err(1,length(err)+1)=2; end; end;
if nargin<3, D = 2; else
    if length(D)~=1; err(1,length(err)+1)=3; end; end;
if nargin<4, XVmin = [-2 -2];else
    if length(XVmin)~=D; err(1,length(err)+1)=4; end; end;
if nargin<5, XVmax = [2 2]; else
    if length(XVmax)~=D; err(1,length(err)+1)=5; end; end;
if nargin<6, y=[]; end;
if nargin<7, NP = 10*D; else
    if length(NP)~=1; err(1,length(err)+1)=7; end; end;
if nargin<8, itermax = 200; else
    if length(itermax)~=1; err(1,length(err)+1)=8; end; end;
if nargin<9, F = 0.8; else
    if length(F)~=1; err(1,length(err)+1)=9; end; end;
if nargin<10, CR = 0.5; else
    if length(CR)~=1; err(1,length(err)+1)=10; end; end;
if nargin<11, strategy = 7; else
    if length(strategy)~=1; err(1,length(err)+1)=11; end; end;
if nargin<12, refresh = 10; else
    if length(refresh)~=1; err(1,length(err)+1)=12; end; end;
if length(err)>0
    fprintf(1,'error in parameter %d\n', err);
    %usage('devec3 (string,scalar,scalar,vector,vector,any,integer,integer,scalar,scalar,integer,integer)');
end

if (NP < 5)
    NP=5;
    fprintf(1,' NP increased to minimal value 5\n');
end
if ((CR < 0) | (CR > 1))
    CR=0.5;
    fprintf(1,'CR should be from interval [0,1]; set to default value 0.5\n');
end
if (itermax <= 0)
    itermax = 200;
    fprintf(1,'itermax should be > 0; set to default value 200\n');
end
refresh = floor(refresh);

formatOut = 'yyyy-mmm-dd_HH_MM_SS';
dateString = datestr(datetime('now'),formatOut);
s = strcat('Log/BRBinsDE_customeRWAW_',num2str(instance_num),'_',dateString,'.txt');
fid_x1 = fopen (s, 'w');
%s = strcat('Log/bestvalue_customeRWAW_',dateString,'.txt');
%fid_x2 = fopen (s, 'w');
s = strcat('Log/BRBinsDE_customeRWAW_',num2str(instance_num),' _IterationCount_',dateString,'.txt');
fid_x3 = fopen (s, 'w');

conseQuentRef=brbConfigdata.conseQuentRef;
numOfAttrWeight=    brbConfigdata.numOfAttrWeight;
numOfconRefval=    brbConfigdata.numOfconRefval;
%input= brbConfigdata.input;
%numOfAntecedentsRefVals=    brbConfigdata.numOfAntecedentsRefVals;
%outputOpti=  brbConfigdata.outputOpti;
%observedOutput=    brbConfigdata.observedOutput;
%    brbConfigdata.transformedRefVal=transformedRefVal;
%    brbConfigdata.rulebase=rulebase;
%sizeOfData=   brbConfigdata.sizeOfData;
numOfVariables=   brbConfigdata.numOfVariables;
numOfRuleWeight=   brbConfigdata.numOfRuleWeight;
numOfbeliefDegrees=    brbConfigdata.numOfbeliefDegrees;

%-----Initialize population and some arrays-------------------------------

pop = zeros(NP,D); %initialize pop to gain speed

%----pop is a matrix of size NPxD. It will be initialized-------------
%----with random values between the min and max values of the---------
%----parameters-------------------------------------------------------

for i=1:NP
    %pop(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
    t=zeros(1,D);
    t(1:numOfAttrWeight+numOfRuleWeight) = XVmin(1:numOfAttrWeight+numOfRuleWeight) + rand(1,numOfAttrWeight+numOfRuleWeight).*(XVmax(1:numOfAttrWeight+numOfRuleWeight) - XVmin(1:numOfAttrWeight+numOfRuleWeight));
    z=t(numOfAttrWeight+numOfRuleWeight+1:numOfVariables-numOfconRefval);
    k=1;
    for l=1:length(z)/length(conseQuentRef)
        r=rand(1,length(conseQuentRef));
        r=r/sum(r);
        z(k:k+length(conseQuentRef)-1)=r;
        k=k+length(conseQuentRef);
    end
    t(numOfAttrWeight+numOfRuleWeight+1:numOfVariables-numOfconRefval)=z;
    pop(i,:)=t;
end

popold    = zeros(size(pop));     % toggle population
val       = zeros(1,NP);          % create and reset the "cost array"
bestmem   = zeros(1,D);           % best population member ever
bestmemit = zeros(1,D);           % best population member in iteration
nfeval    = 0;                    % number of function evaluations

%------Evaluate the best member after initialization----------------------

ibest   = 1;                      % start with first population member
tic
val(1)  = feval(fname,pop(ibest,:),brbConfigdata);
toc
bestval = val(1);                 % best objective function value so far
nfeval  = nfeval + 1;
for i=2:NP % check the remaining members
    val(i) = feval(fname,pop(i,:),brbConfigdata);
end
for i=2:NP
    %val(i) = feval(fname,pop(i,:),brbConfigdata);
    nfeval  = nfeval + 1;
    if (val(i) < bestval)           % if member is better
        ibest   = i;                 % save its location
        bestval = val(i);
    end
end
bestmemit = pop(ibest,:);         % best member of current iteration
bestvalit = bestval;              % best value of current iteration

bestmem = bestmemit;              % best member ever

%------DE-Minimization---------------------------------------------
%------popold is the population which has to compete. It is--------
%------static through one iteration. pop is the newly--------------
%------emerging population.----------------------------------------

pm1 = zeros(NP,D);              % initialize population matrix 1
pm2 = zeros(NP,D);              % initialize population matrix 2
pm3 = zeros(NP,D);              % initialize population matrix 3
pm4 = zeros(NP,D);              % initialize population matrix 4
pm5 = zeros(NP,D);              % initialize population matrix 5
bm  = zeros(NP,D);              % initialize bestmember  matrix
ui  = zeros(NP,D);              % intermediate population of perturbed vectors
mui = zeros(NP,D);              % mask for intermediate population
mpo = zeros(NP,D);              % mask for old population
rot = (0:1:NP-1);               % rotating index array (size NP)
rotd= (0:1:D-1);                % rotating index array (size D)
rt  = zeros(NP);                % another rotating index array
rtd = zeros(D);                 % rotating index array for exponential crossover
a1  = zeros(NP);                % index array
a2  = zeros(NP);                % index array
a3  = zeros(NP);                % index array
a4  = zeros(NP);                % index array
a5  = zeros(NP);                % index array
ind = zeros(4);

iter = 1;
oldF=[];
oldCR=[];
oldBestval=[];
while ((iter < itermax) & (bestval > VTR))
    
    popold = pop;                   % save the old population
    %BRB inspired  part%
    popHistory(:,:,iter)=pop;

%        p = gcp();
        for i=1:NP
            f_history(i,iter) = feval(fname,pop(i,:),brbConfigdata);
        end
    %NP=10;
    % To request multiple evaluations, use a loop.
%     for idx = 1:NP
%         %f(idx) = parfeval(p,@magic,1,idx); % Square size determined by idx
%         f1(idx) = parfeval(p,@objFunAllParallel,1,pop(idx,:),i);
%     end
%     %telapsed1=toc(tstart)
%     % Collect the results as they become available.
%     magicResults = cell(1,NP);
%     for idx = 1:NP
%         % fetchNext blocks until next results are available.
%         [completedIdx,value] = fetchNext(f1);
%         magicResults{completedIdx} = value;
%         %fprintf('Got result with index: %d value %d.\n', completedIdx,value);
%         %fprintf('Got result with index: %d value.\n', completedIdx);
%     end
    %telapsed=toc(tstart)
   
    if (iter>1 )
        tmp=(popHistory(:,:,iter)-popHistory(:,:,iter-1)).^2;
        tmp= sum(sum(tmp,2))/NP;
        PC=sqrt(tmp);
        
        tmp=(f_history(:,iter)-f_history(:,iter-1)).^2;
        FC=sqrt(sum(tmp)/NP);
        d11=1-(1+PC)*exp(-PC);
        d12=1-(1+FC)*exp(-FC);
        
%         [fv n]=expstr(PC);
%         fpc=PC*10^n;
%         [fv n]=expstr(FC);
%         ffc=FC*10^n;
        d21=2*d11;
        d22=2*d12;
        %fprintf(fid_x1,"\nPC:%6.4f,FC:%6.4f,d11:%6.4f,d12:%6.4f,d21:%6.4f,d22:%6.4f",PC,FC,d11,d12,d21,d22);
        fprintf(fid_x1,"\nPC:%g,FC:%g,d11:%g,d12:%g,d21:%g,d22:%g",PC,FC,d11,d12,d21,d22);
        %if PC==0 || FC==0
        %    fprintf(fid_x1,",F:%4.2f,CR:%4.2f,Zero\n",F,CR);
        %else
        [newF,newCR]=brbesFrDEv2(d11,d12, d21,d22);
        oldBestval=horzcat(oldBestval,bestval);
        if (length(oldBestval)==10)&& (bestval > VTR)
            if (all(oldBestval-oldBestval(1)))
                newF=0.8;
                newCR=0.8;
            end
            %oldF=[];
            %oldCR=[];
        end     
        %if newF~=0||newCR~=0
        F=newF;
        CR=newCR;
        %end
        
        fprintf(fid_x1,",F:%6.4f,CR:%6.4f",F,CR);
        %end
        data(iter,1)=PC;
        data(iter,2)=FC;
        data(iter,3)=d11;
        data(iter,4)=d12;
        data(iter,5)=d21;
        data(iter,6)=d22;
        data(iter,7)=F;
        data(iter,8)=CR;
        data(iter,9)=bestval;
        data(iter,10)=iter;
        
    end
    
    ind = randperm(4);              % index pointer array
    
    a1  = randperm(NP);             % shuffle locations of vectors
    rt = rem(rot+ind(1),NP);        % rotate indices by ind(1) positions
    a2  = a1(rt+1);                 % rotate vector locations
    rt = rem(rot+ind(2),NP);
    a3  = a2(rt+1);
    rt = rem(rot+ind(3),NP);
    a4  = a3(rt+1);
    rt = rem(rot+ind(4),NP);
    a5  = a4(rt+1);
    
    pm1 = popold(a1,:);             % shuffled population 1
    pm2 = popold(a2,:);             % shuffled population 2
    pm3 = popold(a3,:);             % shuffled population 3
    pm4 = popold(a4,:);             % shuffled population 4
    pm5 = popold(a5,:);             % shuffled population 5
    
    for i=1:NP                      % population filled with the best member
        bm(i,:) = bestmemit;          % of the last iteration
    end
    
    mui = rand(NP,D) < CR;          % all random numbers < CR are 1, 0 otherwise
    
    if (strategy > 5)
        st = strategy-5;		  % binomial crossover
    else
        st = strategy;		  % exponential crossover
        mui=sort(mui');	          % transpose, collect 1's in each column
        for i=1:NP
            n=floor(rand*D);
            if n > 0
                rtd = rem(rotd+n,D);
                mui(:,i) = mui(rtd+1,i); %rotate column i by n
            end
        end
        mui = mui';			  % transpose back
    end
    mpo = mui < 0.5;                % inverse mask to mui
    
    if (st == 1)                      % DE/best/1
        ui = bm + F*(pm1 - pm2);        % differential variation
        ui = popold.*mpo + ui.*mui;     % crossover
    elseif (st == 2)                  % DE/rand/1
        ui = pm3 + F*(pm1 - pm2);       % differential variation
        ui = popold.*mpo + ui.*mui;     % crossover
    elseif (st == 3)                  % DE/rand-to-best/1
        ui = popold + F*(bm-popold) + F*(pm1 - pm2);
        ui = popold.*mpo + ui.*mui;     % crossover
    elseif (st == 4)                  % DE/best/2
        ui = bm + F*(pm1 - pm2 + pm3 - pm4);  % differential variation
        ui = popold.*mpo + ui.*mui;           % crossover
    elseif (st == 5)                  % DE/rand/2
        ui = pm5 + F*(pm1 - pm2 + pm3 - pm4);  % differential variation
        ui = popold.*mpo + ui.*mui;            % crossover
    end
    %-----Apply constraints--------------------------------------------------
    
    for i=1:NP
        constrainUnstaisfied=0;
        t=ui(i,:);
        %fprintf ('Attribute Weights\n');
        %fprintf ('%2.2f ', t(1:numOfAttrWeight) );
        c=t(1:numOfAttrWeight);
        if sum(c< 0) >0
            constrainUnstaisfied=1;
        end
        if sum(c> 1) >0
            constrainUnstaisfied=1;
        end
        %fprintf ('\nRuleWeights\n');
        %fprintf ('%2.2f ', t(numOfAttrWeight+1:numOfAttrWeight+numOfRuleWeight) );
        c=t(numOfAttrWeight+1:numOfAttrWeight+numOfRuleWeight);
        if sum(c<0) >0
            constrainUnstaisfied=1;
        end
        if sum(c>1) >0
            constrainUnstaisfied=1;
        end
        %fprintf ('\nBelief Degrees\n');
        z=t(numOfAttrWeight+numOfRuleWeight+1:numOfVariables-numOfconRefval);
        j=1;
        for l=1:length(z)/length(conseQuentRef)
            %fprintf ('%2.2f ', z(j:j+length(conseQuentRef)-1));
            c=z(j:j+length(conseQuentRef)-1);
            if sum(c<0) >0
                constrainUnstaisfied=1;
            end
            if sum(c>1) >0
                constrainUnstaisfied=1;
            end
            %fprintf ('\n');
            j=j+length(conseQuentRef);
        end
        %fprintf ('\nconstrainUnstaisfied =%2.2f \n',constrainUnstaisfied);
        if constrainUnstaisfied==1
            t(1:numOfAttrWeight+numOfRuleWeight) = XVmin(1:numOfAttrWeight+numOfRuleWeight) + rand(1,numOfAttrWeight+numOfRuleWeight).*(XVmax(1:numOfAttrWeight+numOfRuleWeight) - XVmin(1:numOfAttrWeight+numOfRuleWeight));
            z=t(numOfAttrWeight+numOfRuleWeight+1:numOfVariables-numOfconRefval);
            k=1;
            for l=1:length(z)/length(conseQuentRef)
                %fprintf ('%2.2f ', z(k:k+length(conseQuentRef)-1));
                r=rand(1,length(conseQuentRef));
                r=r/sum(r);
                z(k:k+length(conseQuentRef)-1)=r;
                %fprintf ('\n');
                k=k+length(conseQuentRef);
            end
            t(numOfAttrWeight+numOfRuleWeight+1:numOfVariables-numOfconRefval)=z;
        end
        ui(i,:)=t;
    end
    %-----Select which vectors are allowed to enter the new population------------
    tempval=zeros(1,NP);
    for i=1:NP
        tempval(i) = feval(fname,ui(i,:),brbConfigdata);
        nfeval  = nfeval + 1;
    end
    for i=1:NP
        %tempval = feval(fname,ui(i,:));   % check cost of competitor
        %nfeval  = nfeval + 1;
        if (tempval(i) <= val(i))  % if competitor is better than value in "cost array"
            pop(i,:) = ui(i,:);  % replace old vector with new one (for new iteration)
            val(i)   = tempval(i);  % save value in "cost array"
            
            %----we update bestval only in case of success to save time-----------
            if (tempval(i) < bestval)     % if competitor better than the best one ever
                bestval = tempval(i);      % new best value
                bestmem = ui(i,:);      % new best parameter vector ever
            end
        end
    end %---end for imember=1:NP
    
    bestmemit = bestmem;       % freeze the best member of this iteration for the coming
    % iteration. This is needed for some of the strategies.
    
    %----Output section----------------------------------------------------------
    
    if (refresh > 0)
        if (rem(iter,refresh) == 0)
            %xplt(NP,pop,0,0)
            fprintf(1,'Iteration: %d,  Best: %f,  F: %f,  CR: %f,  NP: %d\n',iter,bestval,F,CR,NP);
        end
        fprintf(fid_x3,'Iteration: %d,  Best: %f,  F: %f,  CR: %f,  NP: %d\n',iter,bestval,F,CR,NP); 
    end
    
    iter = iter + 1;
   
end %---end while ((iter < itermax) ...
fclose(fid_x3);
fclose(fid_x1);