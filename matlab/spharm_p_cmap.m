function [RGB] = cranio_spharm_p_cmap(P)

H = cat(2, linspace(0, 2/3, 255), 2/3);

PValueCMAPX = cat(2, linspace(0, 0.05, 255), 1);

IH = interp1(PValueCMAPX, H, P);
RGB = hsv2rgb(cat(2, IH(:), ones(numel(IH), 1), ones(numel(IH), 1)));