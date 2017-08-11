function [V, F, CData] = tube_surface(XYZ, Data, R)

% [V, F, CData] = tube_surface(XYZ, Data, R)
%
% DESCRIPTION
%	Returns a surface (patch) of vertices and faces that is a tube around the vertices XYZ
%	of radius R with the colordata derived from Data, so that each vertex of the circle around XYZ(I, :) will
%	have CData(I)
% R = 1 if not given

if(numel(Data) ~= size(XYZ, 1))
	error('The number of elements in Data must equal the number of vertices in XYZ');
end

if nargin < 3
	R = 1;
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
%a
%b
%c

normb = norm(b); if normb~=0, b = b/norm(b); end
%keyboard;
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

TangentVectors = XYZ - circshift(XYZ, [1 0]);

speed = sqrt(sum(TangentVectors .* TangentVectors, 2));

BiNormalVectors = bsxfun(@rdivide, cross(TangentVectors, unitnormals, 2), speed);
% hold off;
% plot3(XYZ(:, 1), XYZ(:, 2), XYZ(:, 3), '*-');
% hold on;
% h = quiver3(XYZ(:, 1), XYZ(:, 2), XYZ(:, 3), unitnormals(:, 1), unitnormals(:, 2), unitnormals(:, 3));
% set(h, 'Color', 'r');
% h = quiver3(XYZ(:, 1), XYZ(:, 2), XYZ(:, 3), BiNormalVectors(:, 1), BiNormalVectors(:, 2), BiNormalVectors(:, 3));
% set(h, 'Color', 'b');

N = 20;
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
CData = repmat(Data(:)', N, 1);
CData = CData(:);

% patch('Faces', F, 'Vertices', V, 'FaceColor', 'r', 'EdgeColor', 'k');
% axis equal;
% keyboard;


function c=crossSimple(a,b)

c(1) = b(3)*a(2) - b(2)*a(3);
c(2) = b(1)*a(3) - b(3)*a(1);
c(3) = b(2)*a(1) - b(1)*a(2);
