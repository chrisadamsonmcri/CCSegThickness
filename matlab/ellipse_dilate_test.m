clear;

D = false(101, 101);
D(50, 50) = 1;

V = radial_dilate(D, [7, 11]);

subplot 121;
imshow(D, []);
subplot 122;
imshow(V, []);