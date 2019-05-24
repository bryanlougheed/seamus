function seamus_run(simstart, siminc, simend, btinc, fpcm, realD14C, blankbg, adpoints, bdpoints, savename, varargin)
%seamus_run(simstart, siminc, simend, btinc, fpcm, realD14C, blankbg, adpoints, bdpoints, savename)
%
%   Sediment simulation module of D14C-enabled SEdiment AccuMUlation Simulator (SEAMUS)
%   B.C. Lougheed, 2019
%   bryan.lougheed@geo.uu.se
%
%   This module runs the actual sediment simulation and outputs variables to a .mat file.
%
%   Required input:
%   ====================
%   simstart = simulation start calendar year BP (whole year, integer)
%   siminc   = simulation timestep (whole year, integer)
%   simend   = simulation end calendar year BP (whole year, integer)
%   btinc    = bioturbate every x years (integer, multiple of or divisible by siminc)
%   fpcm     = specimens per cm (archive capacity, higher = better stats but more memory. integer)
%   realD14C = 0 or 1. Include IntCal D14C (readld14c = 1) or assume constant D14C (realD14C = 0)
%   blankbg  = Blank value (14C years) to assign to single specimens.
%   adpoints = age-depth points. n x 2 matrix. Col 1 = age ka, Col 2 = core depth (cm).
%   bdpoints = bioturbation depth. n x 2 matrix. Col 1 = age ka, Col 2 = bioturbation depth (cm)
%   savename = name to give output .mat file (string, prefix with directory if necessary)
%
%   Optional input:
%   ====================
%   'calcurve' = 14C calibration curve with which to assign 14C activities to single specimens.
%                Refers to corresponding .14c file in matlab search path. Default = 'Marine13'
%                e.g., 'calcurve','Marine13'
%
%   'speciesA' = Temporal abundance of Species A. n by 2 matrix.
%                Col 1 = age ka, Col 2 = abundance of total sediment flux, normalised between 0 and 1
%                Default is constant abundance of 1.
%                e.g., 'speciesA', MATRIX
%
%   'speciesB' = Temporal abundance of Species B. n by 2 matrix.
%                Col 1 = age ka, Col 2 = abundance of total sediment flux, normalised between 0 and 1
%                Default is constant abundance of 0.
%                e.g., 'speciesB', MATRIX
%
%                Note: SpeciesA and SpeciesB must never sum to more than 1.
%
%   'resageA'  = Temporal reservoir age for Species A. n by 2 matrix.
%                Col 1 = age ka, Col 2 = res age in 14C years.
%                Default is constant res age of 0.
%                e.g., 'resageA', MATRIX
%
%   'offsetA'  = Temporal cal age offset from cal curve for Species A. n by 2 matrix.
%                Col 1 = age ka, Col 2 = age offset in years.
%                e.g., 'offsetA', MATRIX
%                Default is constant offset of 0.
%
%   'resageB'  = As for resageA but for Species B
%
%   'offsetB'  = As for offsetA but for Species B.
%
%   'carrierA' = Desired non-radioactive carrier signal(s) (e.g. d18O, etc) for Species A. n by x+1 matrix.
%                Col 1 = age ka, Col 2 = carrier 1, Col 3 = carrier 2, etc.
%                Default is no carrier signals.
%                e.g. 'carrierA', MATRIX
%
%   'carrierB' = Desired non-radioactive carrier signal(s) (e.g. d18O, etc) for Species B. n by x+1 matrix,
%                where x is number of carriers.
%                Col 1 = age ka, Col 2 = carrier 1, Col 3 = carrier 2, etc.
%                Default is no carrier signals.
%                e.g. 'carrierB', MATRIX
%
%   'do32bit'  = Reduce memory requirements in half by saving variables as 32 bit.
%                0 = 64 bit (default), 1 = 32 bit
%                e.g. 'do32bit',1
%
%   Output:
%   ====================
%   A .mat file called savename.mat, containing the following output variables
%   of shared dimensions, whereby the nth position in each variable corresponds
%   to the same specimen:
%
%   depths          = m by 1 matrix of specimen final core depths (cm)
%   depths_original = m by 1 matrix of specimen original deposition depths (cm)
%   cycles          = m by 1 matrix of specimen bioturbation cycles (# of cycles)
%   types           = m by 1 matrix of specimen type identifier (speciesA = 0, speciesB = 1)
%   ages            = m by 1 matrix of specimen ages (age in calendar years before 1950)
%   foram14c        = m by 1 matrix of specimen 14C ages (14C age BP)
%   foramfmc        = m by 1 matrix of specimen 14C activity (fMC)
%   carrierA        = m by x matrix of carrierA signal(s) (x = # of signals)
%   carrierB        = m by x matrix of carrierB signal(s) (x = # of signals)

% Optional parameters input parser (parse varargin)
p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = false;
p.FunctionName='seamus_run';

defaultspeciesA =      [ 0    1
	                     100  1 ];
defaultresageA =       [ 0    0
	                     100  0 ];
defaultoffsetA =       [ 0    0
	                     100  0 ];
defaultspeciesB =      [ 0    0
	                     100  0 ];
defaultresageB =       [ 0    0
	                     100  0 ];
defaultoffsetB =       [ 0    0
	                     100  0 ];
defaultcarrierA = [];
defaultcarrierB = [];
defaultcalcurve = 'Marine13';
defaultdo32bit = 0;

addParameter(p,'speciesA',defaultspeciesA,@isnumeric);
addParameter(p,'resageA',defaultresageA,@isnumeric);
addParameter(p,'speciesB',defaultspeciesB,@isnumeric);
addParameter(p,'resageB',defaultresageB,@isnumeric);
addParameter(p,'offsetA',defaultoffsetA,@isnumeric);
addParameter(p,'offsetB',defaultoffsetB,@isnumeric);
addParameter(p,'carrierA',defaultcarrierA,@isnumeric);
addParameter(p,'carrierB',defaultcarrierB,@isnumeric);
addParameter(p,'calcurve',defaultcalcurve,@isstr);
addParameter(p,'do32bit',defaultdo32bit,@isnumeric);

parse(p,varargin{:});
speciesA = p.Results.speciesA;
resageA = p.Results.resageA;
speciesB = p.Results.speciesB;
resageB = p.Results.resageB;
offsetA = p.Results.offsetA;
offsetB = p.Results.offsetB;
carrierA = p.Results.carrierA;
carrierB = p.Results.carrierB;
calcurve = p.Results.calcurve;
do32bit = p.Results.do32bit;
% /input parser


% calc simulation values for every sim step
disp('Prepping simulation variables')
simstart = round(simstart);
simend = round(simend);
simsteps = simstart:-siminc:simend;
btincs = simstart:-btinc:simend;
fspecA = interp1(speciesA(:,1)*1000, speciesA(:,2), simsteps,'linear','extrap');
fspecB = interp1(speciesB(:,1)*1000, speciesB(:,2), simsteps,'linear','extrap');
Aoffsets = interp1(offsetA(:,1)*1000, offsetA(:,2), simsteps,'linear','extrap');
Boffsets = interp1(offsetB(:,1)*1000, offsetB(:,2), simsteps,'linear','extrap');
resagesA = interp1(resageA(:,1)*1000, resageA(:,2), simsteps,'linear','extrap');
resagesB = interp1(resageB(:,1)*1000, resageB(:,2), simsteps,'linear','extrap');
addepths = interp1(adpoints(:,1)*1000, adpoints(:,2), simsteps,'linear','extrap');
biodepths =  interp1(bdpoints(:,1)*1000, bdpoints(:,2), simsteps,'linear','extrap');
sedrates = diff(addepths) ./ diff(simsteps);
sedrates = [sedrates(1) sedrates];
nforams = round(fpcm ./ ((1./sedrates./siminc))); % nforams per timestep (scaled to sedrate and core capacity)
clear sedrates % no longer needed, free up memory
stinds = cumsum(nforams) - nforams + 1; % startindex for each timestep
if isempty(carrierA) == 0
	carrierAin = NaN(numel(simsteps),size(carrierA,2)-1);
	for i = 2:size(carrierA,2)
		ind = ~isnan(carrierA(:,i));
		carrierAin(:,i-1) = interp1(carrierA(ind,1)*1000, carrierA(ind,i), simsteps,'linear'); % don't extrap
	end
	clear carrierA
end
if isempty(carrierB) == 0
	carrierBin = NaN(numel(simsteps),size(carrierB,2)-1);
	for i = 2:size(carrierB,2)
		ind = ~isnan(carrierB(:,i));
		carrierBin(:,i-1) = interp1(carrierB(ind,1)*1000, carrierB(ind,i), simsteps,'linear'); % don't extrap
	end
	clear carrierB
end

% prep core simualation output variables
if do32bit == 0
	ages = repelem(simsteps,nforams);
	depths = repelem(addepths,nforams);
	depths_original = depths;
	resagesA = repelem(resagesA,nforams);
	resagesB = repelem(resagesB,nforams);
	Aoffsets = repelem(Aoffsets,nforams);
	Boffsets = repelem(Boffsets,nforams);
	types = NaN(size(depths));
	cycles = zeros(size(depths));
elseif do32bit == 1
	ages = single(repelem(simsteps,nforams));
	depths = single(repelem(addepths,nforams));
	depths_original = depths;
	resagesA = single(repelem(resagesA,nforams));
	resagesB = single(repelem(resagesB,nforams));
	Aoffsets = single(repelem(Aoffsets,nforams));
	Boffsets = single(repelem(Boffsets,nforams));
	types = single(NaN(size(depths)));
	cycles = single(zeros(size(depths)));
end

% populate species types for each simstep according to abundance
% quite fast in loop, but can probably be vectorised in future version
for i = 1:length(simsteps)
	types(stinds(i) : stinds(i)+round(fspecA(i)*nforams(i))-1) = 0;
	types(stinds(i)+round(fspecA(i)*nforams(i)) : stinds(i)+round(fspecA(i)*nforams)+round(fspecB(i)*nforams(i))-1) = 1;
end

% RUN DEPTH BIOTURBATION SIMULATION
% can't be vectorised, must run in loop
% because each step dependent on previous step
disp('Running sedimentation simulation')
updateshow = 1:(length(btincs)-1)/10:length(btincs);
updateshow = updateshow(2:end);
updatedisp = 0;
for i = 1:length(btincs)
	% get corresponding simstep index
	s = find(simsteps == btincs(i));
	% find all forams in the contemporaneous bioturbation depth
	ind = find(depths >= addepths(s) & depths < addepths(s) + biodepths(s));
	% bioturbate the depth values with uniform mixing in bioturbation depth
	depths(ind) = rand(length(ind),1)*biodepths(s) + addepths(s);
	% keep track of bioturbation cycles per foram
	cycles(ind) = cycles(ind) + 1;
	if ismember(i,updateshow) == 1
		updatedisp = updatedisp + 10;
		disp([num2str(round(updatedisp)),'%']);
	end
end
clear addepths biodepths simsteps % free up some memory
% --/ END DEPTH BIOTURBATION SIMULATION

% ASSIGN A 14C ACTIVITY TO EACH FORAM
disp('Calculating 14C ages of single specimens')
if realD14C == 0 % if assuming constant d14C, convert real years directly to activity (use real half-life)
	foramfmc(types == 0) = exp([ages(types == 0)+Aoffsets(types == 0)-resagesA(types == 0)]/-8267);
	foramfmc(types == 1) = exp([ages(types == 1)+Boffsets(types == 1)-resagesB(types == 1)]/-8267);
elseif realD14C == 1 % if including cal curve d14C
	File = fopen([calcurve,'.14c']);
	Contents = textscan(File,'%f %f %f %f %f','headerlines',11,'delimiter',',');
	fclose(File);
	curvecal = flipud(Contents{1});
	curve14c = flipud(Contents{2});
	%curved14c = flipud(Contents{4}); % not required
	if do32bit == 0
		foram14c = 70000*ones(size(ages)); % super blank, replaced with chosen blank later
	elseif do32bit == 1
		foram14c = single(70000*ones(size(ages)));
	end
	hicurvecal = curvecal(1):siminc:curvecal(end);
	hicurve14c = interp1(curvecal, curve14c, hicurvecal);
	uages = unique(ages);
	updateshow = 1:(length(uages)-1)/10:length(uages);
	updateshow = updateshow(2:end);
	updatedisp = 0;
	for i = 1:length(uages)
		% Species A
		ind = find(ages == uages(i) & types == 0);
		if isempty(ind) ~=1 && uages(i) <= max(hicurvecal) - Aoffsets(ind(1)) % All Aoffsets of timestep are the same, so use first one
			foram14c(ind) = hicurve14c(hicurvecal == uages(i) + Aoffsets(ind(1)) ) + resagesA(ind(1)); 
		end
		% Species B
		ind = find(ages == uages(i) & types == 1);
		if isempty(ind) ~=1 && uages(i) <= max(hicurvecal) - Boffsets(ind(1)) % ditto
			foram14c(ind) = hicurve14c(hicurvecal == uages(i) + Boffsets(ind(1)) ) + resagesB(ind(1));
		end
		
		if uages(i) > max(hicurvecal)
			disp('100%')
			break
		end
		
		if ismember(i,updateshow) == 1
			updatedisp = updatedisp + 10;
			disp([num2str(round(updatedisp)),'%']);
		end
	end
	foramfmc = exp(foram14c/-8033); % convert IntCal 14C age to activity (use Libby half-life)
end
clear hicurvecal hicurve14c % free up some memory
blankbgfmc = exp(blankbg/-8033); % convert chosen 14C age blank to activity blank (use Libby half-life)
foramfmc(foramfmc<blankbgfmc) = blankbgfmc; % include background on all forams
foram14c = -8033*log(foramfmc); % calculate 14C age for each foram (Libby half-life)
% / ASSIGNING OF 14C ACTIVITY

if strcmpi(savename(end-3:end),'.mat') ~= 1
	savename = [savename,'.mat'];
end

% carrier outputs (do last for memory efficiency)
if exist('carrierAin','var') == 1
	carrierA = repelem(carrierAin,nforams,1);
	clear carrierAin
end
if exist('carrierBin','var') == 1
	carrierB = repelem(carrierBin,nforams,1);
	clear carrierBin
end
clear nforams

save(savename,'-v7.3',...%'-nocompression',...
		'depths',...
		'depths_original',...
		'cycles',...
		'types',...
		'ages',...
		'foram14c',...
		'foramfmc',...
		'blankbg',...
	    'carrierA',...
		'carrierB');
	
	
end % end seamus_run