
dof=24

% Preallocate for speed
pairList = [];
valueList = [];

for i = 1:dof
    for j = 1:dof
        if i ~= j && ~isnan(condNum(i,j))
            pairList = [pairList; i, j];
            valueList = [valueList; condNum(i,j)];
        end
    end
end

% Combine into a table for clarity
T = table(pairList(:,1), pairList(:,2), valueList, ...
    'VariableNames', {'DOF_i', 'DOF_j', 'ConditionNumber'});
disp(T)

[~,idx_min]=min(valueList);
 disp(pairList(idx_min,1))
 disp(pairList(idx_min,2))


% Sort the table by ascending condition number
T_sorted = sortrows(T, 'ConditionNumber');

[condList_sorted, sortIdx] = sortrows(T, 'ConditionNumber');
% T_sorted = table(pairList(sortIdx,1), pairList(sortIdx,2), condList_sorted, ...
%     'VariableNames', {'DOF_i', 'DOF_j', 'ConditionNumber'});


% Display the sorted table
disp('Sorted sensor pairs by condition number (ascending):')
disp(T_sorted)


figure
plot(T_sorted.ConditionNumber, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6)
grid on
xlabel('Sensor Pair Index (sorted)')
ylabel('Condition Number')
title('Sorted Condition Numbers for m=2, DOF=8')

% % Annotate the plot with sensor pairs
% hold on
% for idx = 1:height(T_sorted)
%     txt = sprintf('(%d,%d)', T_sorted.DOF_i(idx), T_sorted.DOF_j(idx));
%     text(idx, T_sorted.ConditionNumber(idx)*1.05, txt, ...
%         'HorizontalAlignment', 'center', 'FontSize', 8)
% end


figure
plot(T.ConditionNumber, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6)
grid on
xlabel('Sensor Pair Index (sorted)')
ylabel('Condition Number')
title('Sorted Condition Numbers for m=2, DOF=8')

% % Annotate the plot with sensor pairs
% hold on
% for idx = 1:height(T_sorted)
%     txt = sprintf('(%d,%d)', T_sorted.DOF_i(idx), T_sorted.DOF_j(idx));
%     text(idx, T_sorted.ConditionNumber(idx)*1.05, txt, ...
%         'HorizontalAlignment', 'center', 'FontSize', 8)
% end





