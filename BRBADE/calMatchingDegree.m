
%function matchingDegree=calMatchingDegree(transformedRefValExt,attrWeightExt,noOfRules,ExtRef)
function matchingDegree=calMatchingDegree(transformedRefVal,attrWeight)
% for i=1:size(transformedRefVal,3)
%     if i==1
%         xb3=combvec(transformedRefVal(:,:,i));
%     else
%         xb3=combvec(xb3,transformedRefVal(:,:,i));
%     end
%     
% end
% xb3=xb3';
% matchingDegree = prod(xb3,2);
% %tic
% transformedRefVal=gpuArray(transformedRefVal1);
% attrWeight=gpuArray(attrWeight1);
% for i=1:size(transformedRefVal,3)
%     transformedRefVal(:,:,i)=power(transformedRefVal(:,:,i),attrWeight(i));
%     if i==1
%         xb3=combvec_gpu(transformedRefVal(:,:,i));
%     else
%         xb3=combvec_gpu(xb3,transformedRefVal(:,:,i));
%     end
%     
% end
% 
% 
% xb3=xb3';
% matchingDegree1 = prod(xb3,2);
% matchingDegree=gather(matchingDegree1);
% %toc

%transformedRefVal=gpuArray(transformedRefVal1);
%attrWeight=gpuArray(attrWeight1);
X=transformedRefVal;
Y=attrWeight;
% for i=1:size(transformedRefVal,3)
%     transformedRefVal(:,:,i)=power(transformedRefVal(:,:,i),attrWeight(i))
%     if i==1
%         xb3=combvec(transformedRefVal(:,:,i))
%     else
%         xb3=combvec(xb3,transformedRefVal(:,:,i))
%     end
%     
% end
% 
% 
% xb3=xb3'
% matchingDegree = prod(xb3,2);

% X = rand(310, 5, 5);
% Y = rand(1, 5);

sX = size(X, 2);
for i = 1:size(X, 3)
   X(:,:,i) = power(X(:,:,i), Y(i));
   if i == 1
     xb3 = X(:,:,i);  % combvec(X) is the same as X
   else
     s1  = size(xb3, 2);
     xb3 = [repmat(xb3, 1, sX); repelem(X(:,:,i), 1, s1)];
   end
end
xb3=xb3';
matchingDegree = prod(xb3,2);

return
end