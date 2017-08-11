function [V, F, CData] = cc_seg_tube_load_cdata(Values, InterpMethod)

% [V, F, CData] = cc_seg_tube_load_cdata(Values)
%
% DESCRIPTIONS
%	Loads the ideal callosum and constructs the tube surface as vertices (V) and faces (F).
%	Also returns CData, which is Values interpolated for each vertex.
%	Utility script for cc_seg_tube_p_display and cc_seg_tube_scalar_display
%
% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
% See LICENSE for full license information.
%

if nargin < 2
	InterpMethod = 'linear';
end

% for the Witelson parcellations I think I need to make a separate script
% that parcellates the main image and then generates the value vector based
% on the locations of the streamlines
load ideal_callosum_interp_points_curved;

M = InterpPoints > 0;
InterpRange = [min(InterpPoints(M)); max(InterpPoints(M))];

% the 1 - is due to the nodes being the wrong way around on the
% interpPoints image
InterpPoints(M) = 1 - (InterpPoints(M) - InterpRange(1)) ./ (InterpRange(2) - InterpRange(1));
InterpPoints(M) = InterpPoints(M) .* (length(Values) - 1) + 1;

B = bwboundaries(InterpPoints > 0);
B{1} = B{1}(1:10:end, :);
I = sub2ind(size(InterpPoints), B{1}(:, 1), B{1}(:, 2));

% NaN values are handled by resetting to zero

T = Values;
if(any(isnan(T)))
	T(isnan(T)) = 0;
end

PAtB = interp1(1:length(Values), T, InterpPoints(I), InterpMethod);

[V, F, CData] = tube_surface([B{1}(:, 2), B{1}(:, 1), zeros(size(B{1}, 1), 1)], PAtB, 50);

V(:, 2) = -V(:, 2);
V = V(:, [1 3 2]);
