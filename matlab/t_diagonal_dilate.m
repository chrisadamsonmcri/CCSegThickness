clear;

BW = zeros(100, 100);
BW(50, 50) = 1;
F = diagonal_dilate(BW, [90, 30], -3 / 9);

subplot 121;
imshow(BW, []);
subplot 122;
imshow(F, []);
size(BW)
size(F)