function newStructure=StructureOptimizationv1(brbStructs,brbParas,brbConfigdatas)

newStructure=[];
theta_new=[];%zeros(size(brbStructs(1).antRefval,1),size(brbStruct(1).antRefval,2));
ind=randperm(length(brbStructs));
if length(ind)==1
    k=ind(1);
    l=ind(1);
else
    k=ind(1);
    l=ind(2);
end

f_k=objFunAllParallelv6(cell2mat(brbParas(k)),brbConfigdatas(k));
f_l=objFunAllParallelv6(cell2mat(brbParas(l)),brbConfigdatas(l));
theta_k=brbStructs(k).antRefval;

if f_k>f_l
    tmp=brbStructs(k).antRefval;
    brbStructs(k).antRefval=brbStructs(l).antRefval;
    brbStructs(l).antRefval=tmp;
    %              for i=1:length(brbStructs(k).antRefval)
    %                  brbStructs(k).antRefval{i}=horzcat(brbStructs(k).antRefval{i},brbStructs(l).antRefval{i});
    %                  brbStructs(l).antRefval{i}=brbStructs(k).antRefval{i};
    %              end
end
%else
for i=1:size(brbStructs(k).antRefval,1)
    theta_k=cell2mat(brbStructs(k).antRefval(i,:));
    theta_l=cell2mat(brbStructs(l).antRefval(i,:));
    j_k=length(theta_k);
    j_l=length(theta_l);
    theta_new=[];
    if j_k>=j_l
        if length(theta_new)==0
            theta_new=theta_l;
        else
            theta_new=horzcat(theta_new,theta_l);
        end
        
        new_rf=ones(1,j_l)*99;
        rnd=rand(1,j_l)<0.25;
        theta_new=horzcat(theta_new,new_rf(rnd));
    else
        rnd=rand(1,j_l)<0.25;
        theta_new=horzcat(theta_new,theta_l(rnd));
        if length(theta_new)==0
            theta_new=theta_l(1);
        end
    end
    newStructure.antRefval(i,:)={theta_new};
end
newStructure.antecedent=brbStructs(1).antecedent;
%newStructure.antRefval=brbStructs(1).antRefval;
newStructure.consequent=brbStructs(1).consequent;
newStructure.conRefval=brbStructs(1).conRefval;
%end
end