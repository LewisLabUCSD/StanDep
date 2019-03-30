function [A,index] = linearization_index(A,flag)

cnt = 0;
index = zeros(1,1);
if iscell(A{1,1})
    for i=1:length(A)
        k = length(A{i,1});
        for j=1:k
            cnt = cnt+1;
            index(cnt,1) = i;
            if strcmp(flag,'rows')
                B{cnt,1} = A{i,1}{1,j};
            else
                B{cnt,1} = A{i,1}{j,1};
            end
        end
    end
else
    for i=1:length(A)
        k = length(A{i,1});
        for j=1:k
            cnt = cnt+1;
            index(cnt,1) = i;
            if strcmp(flag,'rows')
                B(cnt,1) = A{i,1}(1,j);
            else
                B(cnt,1) = A{i,1}(j,1);
            end
        end
    end
end
A = B;