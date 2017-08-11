function [V, F, CData] = tube_surface_color_direction(XYZ, R, N, SurfaceClosed)

% [V, F, CData] = tube_surface(XYZ, R)
%
% DESCRIPTION
%	Returns a surface (patch) of vertices and faces that is a tube around the vertices XYZ
%	of radius R with the colordata being the direction of the streamline

% if(numel(Data) ~= size(XYZ, 1))
% 	error('The number of elements in Data must equal the number of vertices in XYZ');
% end

if nargin < 2
	R = 1;
end

if nargin < 3
	N = 20;
end

if nargin < 4
	SurfaceClosed = false;
end



% determine the unit normals, from streamtube
a = XYZ(2,:) - XYZ(1,:);
b = [0 0 1];
c = crossSimple(a,b);
if ~any(c)
   b = [1 0 0];
   c = crossSimple(a,b);
end
b = crossSimple(c,a);

normb = norm(b); if normb~=0, b = b/norm(b); end

unitnormals = zeros(size(XYZ));
unitnormals(1,:) = b;
    
for J = 1:size(XYZ, 1) - 1;
	a = XYZ(J+1,:) - XYZ(J,:);
	c = crossSimple(a,b);
	b = crossSimple(c,a);
	normb = norm(b); if normb~=0, b = b/norm(b); end
	%  b = b*R(j);
	unitnormals(J+1,:) = b;
end
%eyboard;

% determine the binormals, which is B = T x N
% T = tangent
% N = normal

if(SurfaceClosed)
	TangentVectors = XYZ - circshift(XYZ, [1 0]);
else
	TangentVectors = [XYZ(2, :) - XYZ(1, :); ...
	XYZ(3:end, :) - XYZ(2:end - 1, :); ...
	XYZ(end, :) - XYZ(end - 1, :)];
end

speed = sqrt(sum(TangentVectors .* TangentVectors, 2));

BiNormalVectors = bsxfun(@rdivide, cross(TangentVectors, unitnormals, 2), speed);
% hold off;
% plot3(XYZ(:, 1), XYZ(:, 2), XYZ(:, 3), '*-');
% hold on;
% h = quiver3(XYZ(:, 1), XYZ(:, 2), XYZ(:, 3), unitnormals(:, 1), unitnormals(:, 2), unitnormals(:, 3));
% set(h, 'Color', 'r');
% h = quiver3(XYZ(:, 1), XYZ(:, 2), XYZ(:, 3), BiNormalVectors(:, 1), BiNormalVectors(:, 2), BiNormalVectors(:, 3));
% set(h, 'Color', 'b');

%N = 20;
ANGLES = linspace(0, 2 * pi, N + 1);
ANGLES = ANGLES(1:N)';

VCell = cell(size(XYZ, 1), 1);
FCell = cell(size(XYZ, 1), 1);

COSANGLES = cos(ANGLES);
SINANGLES = sin(ANGLES);
for z = 1:size(XYZ, 1)
	VR = (bsxfun(@times, COSANGLES, unitnormals(z, :)) + bsxfun(@times, SINANGLES, BiNormalVectors(z, :))) * R / 2;
	VCell{z} = bsxfun(@plus, XYZ(z, :), VR);
	
	I = (1:N)';
	% for each vertex make faces
	MinusOneLevel = I - 1;
	MinusOneLevel(MinusOneLevel == 0) = N;
	%F = 
	FCell{z} = ...
		[I + (z - 1) * N, ... same level, same column
		 mod(I, N) + 1 + (z - 1) * N, ... one level up, same column
		 mod(I + z * N - 1, size(XYZ, 1) * N) + 1; ... same level, next column
		 I + (z - 1) * N, ... same level, same column
		 mod(I + z * N - 1, size(XYZ, 1) * N) + 1, ... same level, next column
		 mod(MinusOneLevel + z * N - 1, size(XYZ, 1) * N) + 1]; ... one level down, next column
	% if the surface is not closed then take out triangles that connect up
	% the ends
	if(~SurfaceClosed && z == size(XYZ, 1))
		FCell{z} = FCell{z}(all(FCell{z} > N, 2), :);
	end
	%FCell{z} = ...
	%	[I, mod((I + N) - 1, N * size(XYZ, 1)) + 1, mod(mod(I, N) + z * N, N * size(XYZ, 1)) + 1];
	%keyboard;
	% first row, I, I - 1, I + N
	% second row, I, I + N, I + N - 1
end

%streamline(VCell);
F = cat(1, FCell{:});
V = cat(1, VCell{:});

%II = [repmat(1, N, 1); repmat(2, N, 1)];
%CData = repmat(II(:)', size(XYZ, 1), 1);
%CData = repmat(II(:), 1, size(XYZ, 1));
% %CData = repmat(Data(:)', N, 1);
% GradientXYZ = [XYZ(2, :) - XYZ(1, :); ...
% 	XYZ(3:end, :) - XYZ(1:end - 2, :); ...
% 	XYZ(end, :) - XYZ(end - 1, :)];

GradientXYZ = imfilter(XYZ, [1 -8 8 -1]' / 12, 'same', 'conv', 'replicate');
GradientXYZ = imfilter(GradientXYZ, fspecial('gaussian', [15 1], 5), 'same', 'conv', 'replicate');
GradientXYZMAG = sqrt(sum(GradientXYZ .* GradientXYZ, 2));
CData = bsxfun(@rdivide, abs(GradientXYZ), GradientXYZMAG + double(GradientXYZMAG == 0));

IDX = repmat(1:size(XYZ, 1), N, 1);
IDX = IDX(:);

CData = CData(IDX, :);
% make an index array to repeat
%CData = repmat(CData, N, 1);
%CData = CData(:);

% patch('Faces', F, 'Vertices', V, 'FaceColor', 'r', 'EdgeColor', 'k');
% axis equal;
% keyboard;


function c=crossSimple(a,b)

c(1) = b(3)*a(2) - b(2)*a(3);
c(2) = b(1)*a(3) - b(3)*a(1);
c(3) = b(2)*a(1) - b(1)*a(2);
