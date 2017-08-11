clear;

BW = false(512, 512);

imshow(BW, []);

%[X, Y] = getpts;
XY = [222.0000  115.0000; ...
  177.0000  140.0000; ...
  142.0000  180.0000; ...
  113.0000  242.0000; ...
  126.0000  317.0000; ...
  178.0000  370.0000; ...
  267.0000  389.0000; ...
  325.0000  361.0000; ...
  371.0000  273.0000; ...
  291.0000  187.0000];

XY = round(XY);

I = sub2ind(size(BW), XY(:, 2), XY(:, 1));

BW(I) = 1;

SE = strel('disk', 10);

BW = imdilate(BW, SE);

[B] = bwboundaries(BW, 8);
Distances = zeros(length(B));
DistancesMinI = zeros(length(B));
DistancesMinJ = zeros(length(B));

for BoundaryI = 1:length(B) - 1
	for BoundaryJ = BoundaryI + 1:length(B)
		XC = bsxfun(@minus, B{BoundaryI}(:, 1), B{BoundaryJ}(:, 1)');
		YC = bsxfun(@minus, B{BoundaryI}(:, 2), B{BoundaryJ}(:, 2)');
		SZ = size(XC);
		XC = XC(:);
		YC = YC(:);
		[Distances(BoundaryI, BoundaryJ), I] = min(sqrt(XC .* XC + YC .* YC));
		[DistancesMinI(BoundaryI, BoundaryJ), DistancesMinJ(BoundaryI, BoundaryJ)] = ind2sub(SZ, I);
	end
end
Distances = Distances + Distances';
DistancesMinI = DistancesMinI + DistancesMinI';
DistancesMinJ = DistancesMinJ + DistancesMinJ';
if(length(B) > 2)
	[edges, cost] = prim(Distances);
else
	edges = [1, 2];
end

imshow(BW, []);
for z = 1:size(edges, 1)
	line([...
	B{edges(z, 1)}(DistancesMinI(edges(z, 1), edges(z, 2)), 2), ...
	B{edges(z, 2)}(DistancesMinJ(edges(z, 1), edges(z, 2)), 2)], ... 
	[B{edges(z, 1)}(DistancesMinI(edges(z, 1), edges(z, 2)), 1), ...
	B{edges(z, 2)}(DistancesMinJ(edges(z, 1), edges(z, 2)), 1)]);
end