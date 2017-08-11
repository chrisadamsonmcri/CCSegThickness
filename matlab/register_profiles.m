function [RegisteredProfiles, MeanPeaks, AllPeaks, AllReparamX] = register_profiles(P)

RegisteredProfiles = zeros(size(P));

NumNodes = size(P, 1);

%MeanProfile = nanmean(P, 2);

SE = strel('arbitrary', ones(floor(NumNodes / 5), 1));

%Peaks = find(MeanProfile == imdilate(MeanProfile, SE));

%LeftPeak = Peaks(1);
%RightPeak = Peaks(end);

% hardcoded
LeftPeak = round(NumNodes * 0.12);
RightPeak = round(NumNodes * 0.90);
MeanPeaks = [LeftPeak, RightPeak];


%[~, MeanProfileFirstPeak] = max(MeanProfile(1:floor(NumNodes / 2)));
%[~, MeanProfileSecondPeak] = max(MeanProfile(ceil(NumNodes / 2):end));

%plot(MeanProfile);

%hold on;
%plot(MeanProfileFirstPeak, MeanProfile(MeanProfileFirstPeak), '*');
%plot(MeanProfileSecondPeak, MeanProfile(MeanProfileSecondPeak), '*');
%plot(Peaks, MeanProfile(Peaks), '*');

AllPeaks = zeros(size(P, 2), 2);
AllReparamX = zeros(size(P));
for z = 1:size(P, 2)
	T = P(:, z);
	T(isnan(T)) = 0;
	CurPeaks = find(T == imdilate(T, SE));
	
	
	%plot(T);
	%hold on;
	%keyboard;
	if(~isempty(CurPeaks))
		CurLeftPeaks = CurPeaks(CurPeaks < NumNodes / 2);
		CurRightPeaks = CurPeaks(CurPeaks > NumNodes / 2);

		[~, I] = max(T(CurLeftPeaks));
		CurLeftPeak = CurLeftPeaks(I);
		[~, I] = max(T(CurRightPeaks));
		CurRightPeak = CurRightPeaks(I);
		%plot(CurPeaks, T(CurPeaks), '*');
		%CurLeftPeak = CurPeaks(1);
		%CurRightPeak = CurPeaks(end);
		AllPeaks(z, 1) = CurLeftPeak;
		AllPeaks(z, 2) = CurRightPeak;
		%if(CurLeftPeak ~= LeftPeak && CurRightPeak ~= RightPeak)

			ReparamX = zeros(NumNodes, 1);
			ReparamX(1:CurLeftPeak) = linspace(1, LeftPeak, CurLeftPeak);
			ReparamX(CurLeftPeak:CurRightPeak) = linspace(LeftPeak, RightPeak, CurRightPeak - CurLeftPeak + 1);
			ReparamX(CurRightPeak:end) = linspace(RightPeak, NumNodes, NumNodes - CurRightPeak + 1);
			RegisteredProfiles(:, z) = interp1(ReparamX, T, 1:NumNodes);
			AllReparamX(:, z) = ReparamX;
		%end
		%keyboard;
	end
	%plot
	%keyboard;
end