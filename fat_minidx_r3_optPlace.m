
condNum(condNum==0) = NaN;
% [val,idx] = min(condNum')

% Initialize a cell array to store the configurations
minConfigs = cell(size(condNum, 1), 1);

for  dof= 3:12

    [val,idx] = min(condNum')
    row = dof-2


    col=idx(row)
    % Determine the corresponding i, j, k configuration
    count = 0;
    for i = 1:dof
        for j = 1:dof
            if i ~= j
                for k = 1:dof
                    if k ~= i && k ~= j
                        count = count + 1;
                        if count == col
                            minConfigs{row} = [i, j, k];
                            break;
                        end
                    end
                end
            end
        end
    end
end



optPlace = minConfigs;
% Sort each sub-array in ascending order
optPlace = cellfun(@sort, optPlace, 'UniformOutput', false);

% Display the results
for row = 1:size(condNum, 1)
    fprintf('no. DOF %d: Minimum condition number configuration is [%d, %d, %d].\n', row+2, optPlace{row});
end
