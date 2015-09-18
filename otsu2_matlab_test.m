clear;

inputFile = '/mnt/addo/data/cc_seg/louisee/RawT1ReorientCropped/AMY0024_T1.nii.gz';
if exist(inputFile, 'file') ~= 2
	inputFile = '/data/addo/louisee/RawT1ReorientCropped/AMY0024_T1.nii.gz';
end
NII = load_nii(inputFile);

NIIData = NII.img;

X = NIIData(NIIData > 0);
clear NIIData NII;

X = double(X);
X = (X - min(X)) / (max(X) - min(X));
X = uint8(round(X * 255));

%THRESH1, OmegaZeros1, MuZeros1, SigmaZeros1, OmegaOnes1, MuOnes1, SigmaOnes1, OmegaTwos1, MuTwos1, SigmaTwos1, MuOneHistArray1, SigmaB1, SigmaB2 = CCSegPipe.otsu2(X, returnWorkingValues = True)
[THRESH, METRIC] = multithresh(X, 2);