function seamus_pick(matfilein, matfileout, calcurve, pickint, Apickfordate, Bpickfordate, varargin)
% seamus_pick(matfilein, matfileout, calcurve, pickint, Apickfordate, Bpickfordate)
% 
% Picking and dating module of D14C-enabled SEdiment AccuMUlation Simulator (SEAMUS)
% Version 1.20 (2019-11-29)
% B.C. Lougheed, 2019
% bryan.lougheed@geo.uu.se
% 
% This module picks specimen samples from a simulated sediment core in discrete X cm intervals,
% and calculates the 14C AMS ages, calibrated 14C ages and true ages for each cm.
% 
% Required input:
% ===================
% matfilein      = string directing to the mat file created by seamus_run()
% matfileout     = string with the mat file name to save output variables to
% Apickfordate   = number of whole specimens per sample to pick from species A (-1 = pick all)
% Bpickfordate   = number of whole specimens per sample to pick from species B (-1 = pick all)
% calcurve       = cal curve to use for calibrating 14C dates (e.g. 'Marine13')
% pickint        = discrete depth slice size (cm)
% 
% Optional input
% ===================
%
% Aresage = m by 3 matrix of res age to apply to downcore dates of Species A
% col 1 = depth (cm), col 2 = resage (14C yr), col 3 = resageerr (14Cyr))
% Default is res age of zero p/m zero for all samples.
% e.g., 'Aresage',matrix
%
% Bresage = As for Aresage, but for Species B.
%
% Abroken = m by 2 matrix of fraction broken speicmens (won't be picked) for Species A
% (col 1 = depth (cm), col 2 = fraction broken (between 0 and 1)
% Default is no broken forams for all depths (e.g. fraction broken 0).
% e.g., 'Abroken',matrix
%
% Bbroken = As for Abroken, but for Species B.
%
% AMSerror = 1 by 2 matrix. Specify AMS measurement error for a modern (1.0 F14C) measurement, 
% and a measurement one year younger than the measurement blank. Enter in 14C yrs. E.g. [25 400]
% for 25 14C yr and 400 14C yr errors for 1.0 F14C and the blank, respectively.
% Default is [30 200].
% 		
% Output:
% ===================
% mat file with name matfileout, containing the following variables, which are all
% the same dimension.
% 		
% The row position in each output corresponds to the same picked discrete depth sample (DDS).
% 	
% discdepth      = The centre of the DDS (cm).
% Adiscagemed    = Species A median age for the DDS (years).
% Adiscagemean   = Species A mean age of the DDS (years).
% AdiscAMSage    = Species A predicted lab 14C age for the DDS (14C years).
% AdiscAMSerr    = Species A predicted lab 14C error for the DDS (14C years).
% Adisc14Cage    = Species A mean 14C age for the DDS (14C years).
% Adisccalagemed = Species A calibrated age for the DDS (cal years).
% Adiscblank     = Species A number of 14C blank forams (n specimens).
% Adiscwhole	 = Species A specimens in this DDS with cycles number < than this number are considered
% 				   as still whole and were picked for sample.
% Adiscnforam    = Number of whole Species A specimens picked in DDS sample. (n specimens)
% Adisccarmean   = Species A DDS carrier signal mean. Each column corresponds to a signal.
% Adisccarsdev   = Species A DDS carrier signal stdev. Each column corresponds to a signal.
% Bdiscagemed    = Species B median age for the DDS (years).
% Bdiscagemean   = Species B mean age of the DDS (years).
% BdiscAMSage    = Species B predicted lab 14C age for the DDS (14C years).
% BdiscAMSerr    = Species B predicted lab 14C error for the DDS (14C years).
% Bdisc14Cage    = Species B mean 14C age for the DDS (14C years).
% Bdisccalagemed = Species B calibrated age for the DDS (cal years).
% Bdiscblank     = Species A number of 14C blank forams in DDS(n specimens).
% Bdiscwhole     = Species B specimens in this DDS with cycles number < than this number are considered
% 	        	   as still whole and were picked for sample.
% Bdiscnforam    = Number of whole Species B specimens picked in DDS sample. (n specimens)
% Bdisccarmean   = Species B DDS carrier signal mean. Each column corresponds to a signal.
% Bdisccarsdev   = Species B DDS carrier signal stdev. Each column corresponds to a signal.

% Optional parameters input parser (parse varargin)
p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = false;
p.FunctionName='seamus_pick';

defaultAMSerror =     [30 200];

defaultAresage =      [ 0    0  0
	                    100  0  0];
defaultBresage =      [ 0    0  0
	                    100  0  0];

defaultAbroken =      [ 0    0
	                    100  0];
defaultBbroken =      [ 0    0
	                    100  0];

addParameter(p,'AMSerror',defaultAMSerror,@isnumeric);
addParameter(p,'Aresage',defaultAresage,@isnumeric);
addParameter(p,'Bresage',defaultBresage,@isnumeric);
addParameter(p,'Abroken',defaultAbroken,@isnumeric);
addParameter(p,'Bbroken',defaultBbroken,@isnumeric)

parse(p,varargin{:});
Aresage = p.Results.Aresage;
Bresage = p.Results.Bresage;
Abroken = p.Results.Abroken;
Bbroken = p.Results.Bbroken;
AMSerror = p.Results.AMSerror;

% load needed variables from simulation matfile
% m = matfile(matfilein);  % matfile not yet octave compatible
% depths = m.depths;
% cycles = m.cycles;
% types = m.types;
% foram14c = m.foram14c;
% foramfmc = m.foramfmc;
% ages = m.ages;
% blankbg = m.blankbg;
% carrierA = m.carrierA;
% carrierB = m.carrierB;
load(matfilein, 'depths','cycles','types','foram14c','foramfmc','ages','blankbg','carrierA','carrierB') %octave compatible

fromdepth = (round(min(depths/10))*10);
todepth = (round(max(depths/10))*10)-10;
depths2do = [fromdepth:pickint:todepth]';
discdepth = (depths2do + depths2do+pickint)/2;

Areserr =  interp1(Aresage(:,1), Aresage(:,3), discdepth,'linear','extrap'); % has to come first, because reusing var name
Aresage =  interp1(Aresage(:,1), Aresage(:,2), discdepth,'linear','extrap');
Abroken =  interp1(Abroken(:,1), Abroken(:,2), discdepth,'linear','extrap');
Adiscagemed = NaN(numel(discdepth),1);
Adiscagemean = NaN(numel(discdepth),1);
Adisc95_4 = NaN(numel(discdepth),2);
AdiscAMSage = NaN(numel(discdepth),1);
AdiscAMSerr = NaN(numel(discdepth),1);
Adisc14Cage = NaN(numel(discdepth),1);
Adisccalagemed = NaN(numel(discdepth),1);
%Adisc14CPDF = cell(numel(depths2do),iter);
Adisccal68_2 = cell(numel(depths2do),1);
Adisccal95_4 = cell(numel(depths2do),1);
Adiscblank = NaN(numel(discdepth),1);
Adiscwhole = NaN(numel(discdepth),1);
Adisccarmean = NaN(numel(discdepth),size(carrierA,2));
Adisccarsdev = NaN(numel(discdepth),size(carrierA,2));
Adiscnforam = zeros(numel(discdepth),1);

Breserr =  interp1(Bresage(:,1), Bresage(:,3), discdepth,'linear','extrap'); % has to come first, because reusing var name
Bresage =  interp1(Bresage(:,1), Bresage(:,2), discdepth,'linear','extrap');
Bbroken =  interp1(Bbroken(:,1), Bbroken(:,2), discdepth,'linear','extrap');
Bdiscagemed = NaN(numel(discdepth),1);
Bdiscagemean = NaN(numel(discdepth),1);
Bdisc95_4 = NaN(numel(discdepth),2);
BdiscAMSage = NaN(numel(discdepth),1);
BdiscAMSerr = NaN(numel(discdepth),1);
Bdisc14Cage = NaN(numel(discdepth),1);
Bdisccalagemed = NaN(numel(discdepth),1);
%Bdisc14CPDF = cell(numel(depths2do),iter);
Bdisccal68_2 = cell(numel(depths2do),1);
Bdisccal95_4 = cell(numel(depths2do),1);
Bdiscblank = NaN(numel(discdepth),1);
Bdiscwhole = NaN(numel(discdepth),1);
Bdisccarmean = NaN(numel(discdepth),size(carrierB,2));
Bdisccarsdev = NaN(numel(discdepth),size(carrierB,2));
Bdiscnforam = zeros(numel(discdepth),1);

warning('off')

havebox = license('test','statistics_toolbox');
isoctave = exist('OCTAVE_VERSION', 'builtin');

% ---------SPECIES A----------  % VECTORISE THIS STUFF SOME DAY
if isempty(find(types == 0,1)) ~= 1
	for i = 1:length(discdepth) % for each core depth slice
		
		d1 = depths2do(i);
		d2 = d1+pickint;
		
		if isempty(find(depths >= d1 & depths < d2 & types == 0 , 1)) == 1
			continue
		end
		
		if Abroken(i) > 0 % get index of whole forams in this depth interval
			if havebox == 1
				Adiscwhole(i) = prctile(cycles(depths >= d1 & depths < d2 & types == 0),100-(100*Abroken(i)));
			elseif havebox == 0
				Adiscwhole(i) = octprctile(cycles(depths >= d1 & depths < d2 & types == 0),100-(100*Abroken(i)));			
			end
			foramsnow = find(depths >= d1 & depths < d2 & types == 0 & cycles <= Adiscwhole(i));
		else % or assume all forams are whole (index of all forams in this depth)
			Adiscwhole(i) = max(cycles(depths >= d1 & depths < d2 & types == 0));
			foramsnow = find(depths >= d1 & depths < d2 & types == 0);
		end
		
		if Apickfordate > -1
			if  Apickfordate <= numel(foramsnow)
				% ind = randsample(foramsnow, Apickfordate); % randomly pick X number of forams
				ind = foramsnow(randperm(numel(foramsnow),Apickfordate)); % faster and doesn't require stats toolbox
			else
				continue % not enough available for desired sample size, skip to next discrete depth
			end
		elseif Apickfordate == -1
			ind = foramsnow;
			if numel(ind) == 0
				continue % none available, skip to next discrete depth
			end
		end
		
		% number of forams
		Adiscnforam(i) = numel(ind);
		% number of 14C blank forams in the sample
		Adiscblank(i) = numel(find(foram14c(ind) >= blankbg));
		% True mean age
		Adiscagemean(i) = mean(ages(ind));
		% True median age
		Adiscagemed(i) = median(ages(ind));
        % True 95.45 interval
        Adisc95_4(i,1:2) = prctile(ages(ind),[100*(1-erf(2/sqrt(2)))/2 100-100*(1-erf(2/sqrt(2)))/2]);
		% True mean 14C age
		Adisc14Cage(i) = mean(foram14c(ind));
		% AMS 14C age
		meanfmc = mean(foramfmc(ind));
		AdiscAMSage(i) = -8033*log(meanfmc); % convert to 14C years; Libby half-life
		AdiscAMSerr(i) = interp1([1.0 exp((blankbg+1)/-8033)],[AMSerror(1) AMSerror(2)],meanfmc); % typical AMS error, scaled to fmc
		% Calibrated 14C age
		[Adisccal95_4{i}, Adisccal68_2{i}, ~, Adisccalagemed(i)] = seamus_matcal(AdiscAMSage(i), AdiscAMSerr(i), calcurve, 'CalBP','plot',0,'resage',Aresage(i),'reserr',Areserr(i));
		% carrier signals
		if isempty(Adisccarmean) ~= 1
			if isoctave == 0
				Adisccarmean(i,:) = mean(carrierA(ind,:),1,'omitnan');
				Adisccarsdev(i,:) = std(carrierA(ind,:),0,1,'omitnan');
			else
				Adisccarmean(i,:) = nanmean(carrierA(ind,:),1);
				Adisccarsdev(i,:) = nanstd(carrierA(ind,:),0,1);
			end
		end
	end
end

% ---------SPECIES B---------- % VECTORISE THIS STUFF SOME DAY
if isempty(find(types == 1,1)) ~= 1
	for i = 1:length(discdepth) % for each core depth slice
		
		d1 = depths2do(i);
		d2 = d1+pickint;
		
		if isempty(find(depths >= d1 & depths < d2 & types == 1 , 1)) == 1
			continue
		end
		
		if Bbroken(i) > 0 % get index of whole forams in this depth interval
			if havebox == 1
				Bdiscwhole(i) = prctile(cycles(depths >= d1 & depths < d2 & types == 0),100-(100*Bbroken(i)));
			elseif havebox == 0
				Bdiscwhole(i) = octprctile(cycles(depths >= d1 & depths < d2 & types == 0),100-(100*Bbroken(i)));
			end
			foramsnow = find(depths >= d1 & depths < d2 & types == 1 & cycles <= Bdiscwhole(i));
		else % or assume all forams are whole (index of all forams in this depth)
			Bdiscwhole(i) = max(cycles(depths >= d1 & depths < d2 & types == 1));
			foramsnow = find(depths >= d1 & depths < d2 & types == 1);
		end
		
		if Bpickfordate > -1
			if  Bpickfordate <= numel(foramsnow)
				% ind = randsample(foramsnow, Bpickfordate); % randomly pick X number of forams
				ind = foramsnow(randperm(numel(foramsnow),Bpickfordate)); % faster and doesn't require stats toolbox
			else
				continue % not enough available for desired sample size
			end
		elseif Bpickfordate == -1
			ind = foramsnow;
			if numel(ind) == 0
				continue % none available
			end
		end
		
		% number of forams
		Bdiscnforam(i) = numel(ind);
		% number of 14C blank forams in the sample
		Bdiscblank(i) = numel(find(foram14c(ind) >= blankbg));
		% True mean age
		Bdiscagemean(i) = mean(ages(ind));
		% True median age
		Bdiscagemed(i) = median(ages(ind));
		% True 95.45 interval
        Bdisc95_4(i,1:2) = prctile(ages(ind),[100*(1-erf(2/sqrt(2)))/2 100-100*(1-erf(2/sqrt(2)))/2]);
		% True mean 14C age
		Bdisc14Cage(i) = mean(foram14c(ind));
		% AMS 14C age
		meanfmc = mean(foramfmc(ind));
		BdiscAMSage(i) = -8033*log(meanfmc); % convert to 14C years; Libby half-life
		BdiscAMSerr(i) = interp1([1.0 exp((blankbg-1)/-8033)],[AMSerror(1) AMSerror(2)],meanfmc); % typical AMS error, scaled to fmc
		% Calibrated 14C age
		[Bdisccal95_4{i}, Bdisccal68_2{i}, ~, Bdisccalagemed(i)] = seamus_matcal(BdiscAMSage(i), BdiscAMSerr(i), calcurve, 'CalBP','plot',0,'resage',Bresage(i),'reserr',Breserr(i));
		% carrier signals
		if isempty(Adisccarmean) ~= 1
			if isoctave == 0
				Bdisccarmean(i,:) = mean(carrierB(ind,:),1,'omitnan');
				Bdisccarsdev(i,:) = std(carrierB(ind,:),0,1,'omitnan');
			else
				Bdisccarmean(i,:) = nanmean(carrierB(ind,:),1);
				Bdisccarsdev(i,:) = nanstd(carrierB(ind,:),0,1);
			end
		end
	end
end
warning('on')

if strcmpi(matfileout(end-3:end),'.mat') ~= 1
	matfileout = [matfileout,'.mat'];
end

if isoctave == 0
	save(matfileout,'-v7.3', '-nocompression',...
		'discdepth',...
		'Adiscagemed',...
		'Adiscagemean',...
		'Adisc95_4',...
		'AdiscAMSage',...
		'AdiscAMSerr',...
		'Adisc14Cage',...
		'Adisccalagemed',...
		'Adisccal95_4',...
		'Adisccal68_2',...
		'Adiscblank',...
		'Adiscwhole',...
		'Adisccarmean',...
		'Adisccarsdev',...
		'Adiscnforam',...
		'Bdiscagemed',...
		'Bdiscagemean',...
		'Bdisc95_4',...
		'BdiscAMSage',...
		'BdiscAMSerr',...
		'Bdisc14Cage',...
		'Bdisccalagemed',...
		'Bdisccal95_4',...
		'Bdisccal68_2',...
		'Bdiscblank',...
		'Bdiscwhole',...
		'Bdisccarmean',...
		'Bdisccarsdev',...
		'Bdiscnforam');
else
	save(matfileout,...
		'discdepth',...
		'Adiscagemed',...
		'Adiscagemean',...
		'Adisc95_4',...
		'AdiscAMSage',...
		'AdiscAMSerr',...
		'Adisc14Cage',...
		'Adisccalagemed',...
		'Adisccal95_4',...
		'Adisccal68_2',...
		'Adiscblank',...
		'Adiscwhole',...
		'Adisccarmean',...
		'Adisccarsdev',...
		'Adiscnforam',...
		'Bdiscagemed',...
		'Bdiscagemean',...
		'Bdisc95_4',...
		'BdiscAMSage',...
		'BdiscAMSerr',...
		'Bdisc14Cage',...
		'Bdisccalagemed',...
		'Bdisccal95_4',...
		'Bdisccal68_2',...
		'Bdiscblank',...
		'Bdiscwhole',...
		'Bdisccarmean',...
		'Bdisccarsdev',...
		'Bdiscnforam');
end



% Nested function octprctile
function a = octprctile(x, p)
	y = sort(x);
	if size (y,1) == 1
		y = y(:);
	end
	trim = 1 + (size(y,1)-1)*p(:)*0.01;
	delta = (trim - floor(trim))*ones(1,size(y,2));
	a = y(floor(trim), :) .* (1-delta) + y(ceil(trim), :) .* delta;
end



end % end seamus_pick
