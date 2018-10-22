function [bestmem,bestval,nfeval] = DiffEv(fname,VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)

   %-----Check input variables---------------------------------------------
% err=[];
% if nargin<1, error('devec3 1st argument must be function name'); else 
%   if exist(fname)<1; err(1,length(err)+1)=1; end; end;
% if nargin<2, VTR = 1.e-6; else 
%   if length(VTR)~=1; err(1,length(err)+1)=2; end; end;
% if nargin<3, D = 2; else
%   if length(D)~=1; err(1,length(err)+1)=3; end; end; 
% if nargin<4, XVmin = [-2 -2];else
%   if length(XVmin)~=D; err(1,length(err)+1)=4; end; end; 
% if nargin<5, XVmax = [2 2]; else
%   if length(XVmax)~=D; err(1,length(err)+1)=5; end; end; 
% if nargin<6, y=[]; end; 
% if nargin<7, NP = 10*D; else
%   if length(NP)~=1; err(1,length(err)+1)=7; end; end; 
% if nargin<8, itermax = 200; else
%   if length(itermax)~=1; err(1,length(err)+1)=8; end; end; 
% if nargin<9, F = 0.8; else
%   if length(F)~=1; err(1,length(err)+1)=9; end; end;
% if nargin<10, CR = 0.5; else
%   if length(CR)~=1; err(1,length(err)+1)=10; end; end; 
% if nargin<11, strategy = 7; else
%   if length(strategy)~=1; err(1,length(err)+1)=11; end; end;
% if nargin<12, refresh = 10; else
%   if length(refresh)~=1; err(1,length(err)+1)=12; end; end; 
% if length(err)>0
%   fprintf(1,'error in parameter %d\n', err);
%   fprintf('devec3 (string,scalar,scalar,vector,vector,any,integer,integer,scalar,scalar,integer,integer)');    	
% end

% if (NP < 5)
%    NP=5;
%    fprintf(1,' NP increased to minimal value 5\n');
% end
% if ((CR < 0) || (CR > 1))
%    CR=0.5;
%    fprintf(1,'CR should be from interval [0,1]; set to default value 0.5\n');
% end
% if (itermax <= 0)
%    itermax = 200;
%    fprintf(1,'itermax should be > 0; set to default value 200\n');
% end
refresh = floor(refresh);
formatOut = 'yyyy-mmm-dd_HH_MM_SS';
dateString = datestr(datetime('now'),formatOut);
s = strcat('Log/DE_',dateString,'.txt');
fid_x1 = fopen (s, 'w');
% s = strcat('Log/DE_bestvalue_',dateString,'.txt');
% fid_x2 = fopen (s, 'w');
%-----Initialize population and some arrays-------------------------------

pop = zeros(NP,D); %initialize pop to gain speed

%----pop is a matrix of size NPxD. It will be initialized-------------
%----with random values between the min and max values of the---------
%----parameters-------------------------------------------------------

for i=1:NP
   pop(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end

popold    = zeros(size(pop));     % toggle population
val       = zeros(1,NP);          % create and reset the "cost array"
bestmem   = zeros(1,D);           % best population member ever
bestmemit = zeros(1,D);           % best population member in iteration
nfeval    = 0;                    % number of function evaluations

%------Evaluate the best member after initialization----------------------

ibest   = 1;                      % start with first population member
val(1)  = feval(fname,pop(ibest,:),y);
bestval = val(1);                 % best objective function value so far
nfeval  = nfeval + 1;
for i=2:NP                        % check the remaining members
  val(i) = feval(fname,pop(i,:),y);
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
data(iter,1)=0;
     data(iter,2)=0;
     data(iter,3)=0;
     data(iter,4)=0;
     data(iter,5)=0;
     data(iter,6)=0;
     data(iter,7)=F;
     data(iter,8)=CR;
     data(iter,9)=bestval;
     data(iter,10)=iter;
%data=zeros(itermax,9);
while ((iter < itermax) & (bestval > VTR))
%while (iter < itermax)
  %fprintf(fid_x1,[repmat('%2.2f\t', 1, size(pop, 1)) '\n'], pop); 
  popold = pop;                   % save the old population
  
  
  %------------------%
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

%-----Select which vectors are allowed to enter the new population------------
  for i=1:NP
    tempval = feval(fname,ui(i,:),y);   % check cost of competitor
    nfeval  = nfeval + 1;
    if (tempval <= val(i))  % if competitor is better than value in "cost array"
       pop(i,:) = ui(i,:);  % replace old vector with new one (for new iteration)
       val(i)   = tempval;  % save value in "cost array"

       %----we update bestval only in case of success to save time-----------
       if (tempval < bestval)     % if competitor better than the best one ever
          bestval = tempval;      % new best value
          bestmem = ui(i,:);      % new best parameter vector ever
       end
    end
  end %---end for imember=1:NP

  bestmemit = bestmem;       % freeze the best member of this iteration for the coming 
                             % iteration. This is needed for some of the strategies.
  % fprintf(fid_x2,'%d,', iter); 
  % fprintf(fid_x2,[repmat('%2.2f,', 1, size(bestmem, 2)) '\n'], bestmem); 
%----Output section----------------------------------------------------------

  if (refresh > 0)
    if (rem(iter,refresh) == 0)
       fprintf(fid_x1,'Iteration: %d,  Best: %f,  F: %f,  CR: %f,  NP: %d',iter,bestval,F,CR,NP);
       fprintf('Iteration: %d,  Best: %f,  F: %f,  CR: %f,  NP: %d\n',iter,bestval,F,CR,NP);
    end
  end
  data(iter,1)=0;
  data(iter,2)=0;
  data(iter,3)=0;
  data(iter,4)=0;
  data(iter,5)=0;
  data(iter,6)=0;
  data(iter,7)=F;
  data(iter,8)=CR;
  data(iter,9)=bestval;
  data(iter,10)=iter;
  iter = iter + 1;
  
end %---end while ((iter < itermax) ...

ss=sprintf('%s_de.mat',fname)
delete(ss)
save(ss,'data')
% fnamm='ass.mat'