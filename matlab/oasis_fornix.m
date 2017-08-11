clear;

D = dir(fullfile('oasis_database', 'flipped', 'OAS1*.nii.gz'));

[FNames{1:length(D)}] = deal(D.name);

[match, tokens] = regexp(FNames, '^(OAS1_\d+).*\.nii\.gz$', 'match', 'tokens');

Codes = cell(1, length(D));

for z = 1:length(FNames)
	Codes{z} = tokens{z}{1}{1};
end

Codes = unique(Codes);
%R = [176.7955  182.1992  174.2691  101.3193];

for z = 1:length(Codes)
	disp(Codes{z});
	[NIISeg, AVWSeg] = load_nii(fullfile('oasis_database', 'flipped', [Codes{z} '_cc.nii.gz']));
	[NII, AVW] = load_nii(fullfile('oasis_database', 'flipped', [Codes{z} '_msp.nii.gz']));
	OriginalAVWSegClass = class(AVWSeg);
	AVWSeg = (AVWSeg > 0);
	CroppingBox = imdilate(AVWSeg, ones(25, 25));
	[I, J] = find(CroppingBox);
	CroppingI = min(I):max(I);
	CroppingJ = min(J):max(J);
	AVWCropped = AVW(CroppingI, CroppingJ);
	AVWSegCropped = AVWSeg(CroppingI, CroppingJ);
	AVWCropped = double(AVWCropped);
	
	MeanIntensity = mean(double(AVWCropped(AVWSegCropped)));
	VarIntensity = var(double(AVWCropped(AVWSegCropped)));
	StdIntensity = sqrt(VarIntensity);
	
	XC = double(AVWCropped) - MeanIntensity;
	N = 1 ./ (sqrt(2 * pi) * StdIntensity);
	P = N .* exp(-(XC .* XC) ./ 2 ./ VarIntensity);
	P(XC > 0) = N;
	L = bwlabel((P > N / 5) & ~imdilate(AVWSegCropped, ones(3, 3)));
	R = regionprops(L, 'Area', 'PixelIdxList', 'PixelList', 'Centroid');
	CCR = regionprops(AVWSegCropped, 'Centroid');
	XC = bsxfun(@minus, CCR.Centroid, cat(1, R.Centroid));
	
	M = max(abs(XC), [], 2);
	%if(isequal(Codes{z}, 'OAS1_0085'))
	%	keyboard;
	%end
	%keyboard;
	R(M > 20) = [];
	% redo L
	OutIMG = (AVWCropped - min(AVWCropped(:))) ./ (max(AVWCropped(:)) - min(AVWCropped(:)));
	OutIMG = round(OutIMG * 255);
	OutIMG = repmat(uint8(OutIMG), [1, 1, 3]);

	if(~isempty(R))
		[~, I] = max([R.Area]);
		%[ID, JD] = find(L == I);
		[ID, JD] = deal(R(I).PixelList(:, 2), R(I).PixelList(:, 1));
		%OutIMG(:, :, 3) = 0;
		for k = 1:length(ID)
			OutIMG(ID(k), JD(k), :) = [0, 0, 255];
		end
		AVWFornixSeg = zeros(size(AVWSeg), OriginalAVWSegClass);
		T = zeros(size(AVWSegCropped), OriginalAVWSegClass);
		T(R(I).PixelIdxList) = 1;
		AVWFornixSeg(CroppingI, CroppingJ) = T;
		NIIFornix = NIISeg;
		NIIFornix.img = permute(flipdim(AVWFornixSeg, 1), [2, 1]);
		save_nii(NIIFornix, fullfile('oasis_database', 'flipped', [Codes{z} '_fornix.nii.gz']));
	end
	
	[~, ~, ~] = mkdir('fornix_figures');
	imwrite(OutIMG, fullfile('fornix_figures', [Codes{z} '_fornix.png']));
end