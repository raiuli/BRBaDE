function [bestmem,bestval,nfeval] = GADE(fname,VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,strategy,refresh)


%-----Initialize population and some arrays-------------------------------

pop = zeros(NP,D); %initialize pop to gain speed
CR_m=0.5;
F=0.5;
LP=20;
d1=0.01;
d2=0.01;
Zf=[F-d1,F,F+d2];
Zcr=[CR_m-d1,CR_m,CR_m+d2];
nfeval    = 0;                    % number of function evaluations
formatOut = 'yyyy-mmm-dd_HH_MM_SS';
dateString = datestr(datetime('now'),formatOut);
s = strcat('Log/GADE_',dateString,'.txt');
fid_x1 = fopen (s, 'w');
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
nfeval    = 0; 
iter      = 1;


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

MZf = containers.Map('KeyType','double','ValueType','double');
MZcr = containers.Map('KeyType','double','ValueType','double');
iter=1;
data(iter,1)=0;
  data(iter,2)=0;
  data(iter,3)=0;
  data(iter,4)=0;
  data(iter,5)=0;
  data(iter,6)=0;
  data(iter,7)=0;
  data(iter,8)=0;
  data(iter,9)=bestval;
  data(iter,10)=iter;
  iter=iter+1;
while ((iter < itermax) & (bestval > VTR))
  popold = pop;   
  %------------------%
  ind = randperm(4);              % index pointer array
  rot = (0:1:NP-1);
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
  
  
  F=Zf(randi([1,3],NP,1));
  muCR=Zcr(randi([1,3],NP,1));
  CR=normrnd(muCR, 0.2);
  
  mui = rand(NP,D) < CR';
  mpo = mui < 0.5; 
  ui = pm3 + F*(pm1 - pm2);       % differential variation
  ui(find(ui<XVmin))=XVmin(1);
  ui(find(ui>XVmax))=XVmax(1);
  ui = popold.*mpo + ui.*mui;     % crossover
  %-----Select which vectors are allowed to enter the new population------------
 
  RI_F=containers.Map('KeyType','double','ValueType','double');
  RI_CR=containers.Map('KeyType','double','ValueType','double');
  for i=1:NP
    tempval = feval(fname,ui(i,:),y);   % check cost of competitor
    %fprintf(1,"iter:%d,tempval:%d,val(i):%d \n",iter,tempval,val(i))
    nfeval  = nfeval + 1;
    if (tempval <= val(i))  % if competitor is better than value in "cost array"
       pop(i,:) = ui(i,:);  % replace old vector with new one (for new iteration)
       val(i)   = tempval;  % save value in "cost array"
        if isKey(RI_F,F(i))
            RI_F(F(i))=RI_F(F(i))+0;
        else
            RI_F(F(i))=0;
        end
        if isKey(RI_F,CR(i))
            
            RI_CR(CR(i))=RI_F(CR(i))+0;
        else    
            RI_CR(CR(i))=0;
        end
       
       %----we update bestval only in case of success to save time-----------
       if (tempval < bestval)     % if competitor better than the best one ever
          bestval = tempval;      % new best value
          bestmem = ui(i,:);      % new best parameter vector ever
       end
    else
        [fv n]=expstr(tempval);
        fx=val(i)*10^n;
        %fv=0.33;
        if fv>10 || fv<1
            fprintf(1,"error in fv:%g \n",fv)
        end
        %fx=11;
        %if fx>10 || fx<1
        %    fprintf(1,"error in fx:%g \n",fx)
        %end    
        if isKey(RI_F,F(i))
            RI_F(F(i))=RI_F(F(i))+(fx-fv);
        else
            RI_F(F(i))=(fx-fv);
        end
        if isKey(RI_F,CR(i))
            
            RI_CR(CR(i))=RI_F(CR(i))+(fx-fv);
        else    
            RI_CR(CR(i))=(fx-fv);
        end
    end
    
  end %---end for imember=1:NP
  if mod(iter,LP)==0
      for k=keys(RI_F)
        MZf(k{1})=RI_F(k{1})/NP;
      end
      for k=keys(RI_CR)
        MZcr(k{1})=RI_CR(k{1})/NP;
      end

      %fprintf(1,"iter:%d,length(MZf):%d,length(MZcr):%d \n",iter,length(MZf),length(MZcr)),
     v=[];
     k=[];
      for i=keys(RI_F)
        v=[v;MZf(i{1})];
        k=[k;i{1}];
      end
      %[sort_value, sort_index]=sort(v,'descend');
     [sort_value, sort_index]=sort(v,'descend');
      F=k(sort_index(1));
      Zf=[F-d1,F,F+d2];
      v=[];
      k=[];
      for i=keys(RI_CR)
        v=[v;MZcr(i{1})];
        k=[k;i{1}];
      end
      [sort_value, sort_index]=sort(v,'descend');
      CR=k(sort_index(1));
      Zcr=[CR-d1,CR,CR+d2];
  end
  bestmemit = bestmem;       % freeze the best member of this iteration for the coming    
  if (refresh > 0)
    if(iter==1)
       fprintf(fid_x1,'Iteration: %d,  Best: %f,  NP: %d',iter,bestval,NP);
       fprintf('Iteration: %d,  Best: %f,  NP: %d\n',iter,bestval,NP);
    end  
    if (rem(iter,refresh) == 0)
       fprintf(fid_x1,'Iteration: %d,  Best: %f,   NP: %d',iter,bestval,NP);
       fprintf('Iteration: %d,  Best: %f,    NP: %d\n',iter,bestval,NP);
       
%        for n=1:D
%          fprintf(fid_x1,'best(%d) = %f\n',n,bestmem(n));
%        end
       
    end
    if(iter==itermax)
      fprintf(fid_x1,'Iteration: %d,  Best: %f,   NP: %d',iter,bestval,NP);
       fprintf('Iteration: %d,  Best: %f,    NP: %d\n',iter,bestval,NP);
    end 
  end
  data(iter,1)=0;
  data(iter,2)=0;
  data(iter,3)=0;
  data(iter,4)=0;
  data(iter,5)=0;
  data(iter,6)=0;
  data(iter,7)=0;
  data(iter,8)=0;
  data(iter,9)=bestval;
  data(iter,10)=iter;
  iter=iter+1;
end    
ss=sprintf('%s_gade.mat',fname)
delete(ss)
save(ss,'data')