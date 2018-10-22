function result= checkNumOFrefValue(bestStruct,newStruct)
    result=0;
    
    a=0;
    for i=1:length(bestStruct.antRefval)
        a=a+length(bestStruct.antRefval{i});
    end
    %b=prod(size(cell2mat(newStruct.antRefval)));
    b=0;
    for i=1:length(newStruct.antRefval)
        b=b+length(newStruct.antRefval{i});
    end    
    if b<a
        result=1;
    end
end