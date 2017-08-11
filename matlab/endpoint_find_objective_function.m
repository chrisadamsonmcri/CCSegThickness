function [ArcLength, MaskClosed] = endpoint_find_objective_function(MainContour, LeftOffset, RightOffset, RightStartIDX, MaskClosed)

BT = circshift(MainContour, [-LeftOffset, 0]);
%RightStartIDX + Z + RightOffset
Contours.xi = BT(1:RightStartIDX - LeftOffset + RightOffset, 2);
Contours.yi = BT(1:RightStartIDX - LeftOffset + RightOffset, 1);
Contours.xo = BT(RightStartIDX - LeftOffset + RightOffset + 1:size(BT, 1), 2);
Contours.yo = BT(RightStartIDX - LeftOffset + RightOffset + 1:size(BT, 1), 1);

Contours.xo = flipud(Contours.xo);
Contours.yo = flipud(Contours.yo);

%NumContours = 250;
%clf;
%plot(Contours.xo, Contours.yo, 'r', Contours.xi, Contours.yi, 'b');

%plot(Contours.xi, Contours.yi, Contours.xo(end:-1:1), Contours.yo(end:-1:1));
%keyboard;
[~, ~, ~, SolvedImage, ~, ~, ~, ~, ~, ~, ~, MaskClosed] = laplace_get_points_2d_auto_mw(Contours.xi, Contours.yi, Contours.xo(end:-1:1), Contours.yo(end:-1:1), 1, 1, 1, 0, MaskClosed);

%imshow(SolvedImage, []);
%keyboard;
C = contourc(SolvedImage, [0.5, 0.5]);
if(isempty(C))
	ArcLength = 0;
else
	CT = contour2cell(C);
	SZ = cellfun('size', CT, 1);
	[~, I] = max(SZ);
	C = CT{I};
	%disp(num2str(length(CT)));
	ArcLength = arc_length(C);
end