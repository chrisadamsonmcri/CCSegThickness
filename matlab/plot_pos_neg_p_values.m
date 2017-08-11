function [A, B, L] = plot_pos_neg_p_values(PValues, TValues, IDX, GroupNames)

if nargin < 3
	IDX = 1:numel(PValues);
end

I = (TValues > 0);
A = plot(IDX(I), PValues(I), 'r*');
hold on;
B = plot(IDX(~I), PValues(~I), 'b*');
ylim([0 1]);
XL = [min(IDX), max(IDX)];
xlim(XL);
line(XL, [0.05 0.05]);
if nargin > 3
	if(numel(GroupNames) ~= 2 || ~iscellstr(GroupNames))
		error('GroupNames must be a 2 element cell array of strings');
	end
	if(any(I) && any(~I))
		C = {[GroupNames{1} ' > ' GroupNames{2}], [GroupNames{1} ' < ' GroupNames{2}]};
	elseif(any(I))
		C = {[GroupNames{1} ' > ' GroupNames{2}]};
	else
		C = {[GroupNames{1} ' < ' GroupNames{2}]};
	end
	L = legend(C);
end
line([IDX(1), IDX(end)], [0.05, 0.05], 'Color', 'k');
xlabel('Node');
ylabel('\itp');