%% Welcome to the SEAMUS tutorial, a step-by-step hands-on guide to the simulation and plotting of results
% The tutorial is divided into code blocks, so that you can run each code block at a time.

%% (1): Let's load some kind of interesting climate signal

% (1a): load the NGRIP d18O record vs GICC05 time ----------
% Of course, e.g. forams in the ocean don't directly record the NGRIP
% values, but we just use it here as a nice signal to play with.
d = load('Greenland on GICC05 50 yr mean (repeats removed).txt');
ngripage = (d(:,1) - 50)/1000; % Change b2k yrs to b1950 ka
ngrip1 = d(:,3);   % d18O of NGRIP1 (holocene)
ngrip2 = d(:,5);   % d18O of NGRIP2 (older parts)
ngrip1 = ngrip1(isnan(ngrip2)); % Remove NGRIP1 and NGRIP2 overlap
ngrip2 = ngrip2(~isnan(ngrip2)); % Remove NGRIP1 and NGRIP2 overlap
ngrip18O = [ngrip1; ngrip2]; % merge NGRIP1 and NGRIP2
ngrip18O = ngrip18O(ngripage<=50); % limit to past 50 ka
ngripage = ngripage(ngripage<=50); % limit to past 50 ka


%% (2): Let's prepare the input variables for seamus_run sediment simulation

% (2a): General simulation settings for seamus_run --------
simstart = 70000;          % start year (yrs BP)
siminc = 10;               % simulation timestep (yrs)
simend = 0;                % end year (yrs BP)
btinc = 10;                % bioturbation timestep (yrs).
realD14C = 1;              % Assign 'real' D14C to single specimens (1 = yes, 0 = no)
calcurve = 'Marine13';     % Calibration curve to use for assigning 'real' D14C
blankbg = 46806;           % the 14C blank value to assign to single specimens (46806 = higest value in Marine13 = practical blank)
fpcm = 10^3;               % Core capacity in specimens per cm. Let's use 10^3 for now to run quickly (you can change later).
runfile = 'test_run.mat';  % Name of output file from seamus_run

% (2b): Dynamic inputs for seamus_run --------------------

% bioturbation depths   ka   BD (cm)
bdpoints =             [0    10
	                    50   10
	                    100  10];
% The above will result in a constant bioturbation depth
% of 10 cm. You can enter dynamic values by adding extra rows
% to the matrix and the simulation will interpolate/extrapolate
% where needed. Or you can simply load data from external file.


% age-depth points      ka  depth (cm)
adpoints =              [0   0
	                    25  250
	                    50  500]; 
% The above will result in constant SAR of 10 cm/ka.
% You can enter as many age depth points as you like, and
% the simulation will interpolate/extrapolate where needed.
% You could also load an external file instead.

% Species A carrier signal
%                       ka        signal  
carrierA =              [ngripage ngrip18O];
% Here we simply use the NGRIP data that we prepared earlier.
% For carrier signals, the simulation will interpolate, but
% not extrapolate.

% Species A abundance   ka     fraction of specimen flux
speciesA =             [0      1
	                    25     1
						50     1];

%% (3): Run the seamus_run sediment simulation (can take some time)			

seamus_run(simstart, siminc, simend, btinc,	fpcm, realD14C,... 
	blankbg, adpoints, bdpoints, runfile,...
	'calcurve', calcurve,'carrierA',carrierA,'speciesA',speciesA);
% you can add more optional inputs, see documentation

%% (4): Now prepare the input variables for the seamus_pick picking and 14C dating simulaiton

% (4a): General simulation settings for seamus_pick
pickint = 1;                  % thickness of the core slices (cm)
Apickfordate = -1;            % number of Species A specimens to pick per sample (enter -1 for all available)
Bpickfordate = -1;            % number of Species A specimens to pick per sample (enter -1 for all available)
calcurve = 'Marine13';        % calibration curve to use when calibrating dates
pickfile = 'test_pick.mat';   % file to save the picking results to
matfilein = runfile;          % The name of the sediment simulation to analyse (the file saved by seamus_run)


%% (5): Now run the seamus_pick picking and 14C dating simulation (can take some time)

seamus_pick(matfilein, pickfile, calcurve, pickint, Apickfordate, Bpickfordate)
% you can add more optional inputs, see documentation

%% (6) Now we can make some figures. Let's start with a simple downcore figure showing depth vs age

% ---- (6a): Get some stuff to plot --------------------

% tell matlab where the sediment simulation and picking simulation files are
r = matfile(runfile);
p = matfile(pickfile);

% Load single specimen data to plot:
depths = r.depths;
ages = r.ages; % <--- hint: change r.ages to r.foram14c to get the raw 14C ages for the single specimens
types = r.types;
% get only Species A single specimens (types value of 0)
depths = depths(types == 0);
ages = ages(types == 0);

% load the downcore discrete-depth median age from the pickfile
Adiscage = p.Adiscagemed;
discdepth = p.discdepth; % and associated depth slice values
% hint: change p.Adiscagemed to p.Adisccalagemed and you will get the calibrated 14C age, or p.discAMSage for the uncalibrated lab 14C age.

% ---- (6b): Now for the plotting
figure(1)
clf
% plot a cloud of single specimen values
plot(ages/1000,depths,'b.','markersize',1,'color',[0.8 0.8 0.8]) % divide age by 1000 to make ka (easier for the plot labels)
% plot the median age on top
hold on
plot(Adiscage/1000,discdepth,'k-','linewidth',2);
% make the figure nicer
set(gca,'ydir','reverse')
xlabel('ka')
ylabel('Depth (cm)')
grid on

%% (7) Figure of carrier signal vs simulated core depth

% ---- (7a): Get some stuff to plot --------------------

% tell matlab where the sediment simulation and picking simulation files are
r = matfile(runfile);
p = matfile(pickfile);

% Load single specimen data to plot:
depths = r.depths;
carrierA = r.carrierA; % get the carrier signal(s) for Species A
carrierA = carrierA(:,1); % get the first carrierA signal for Species A (in this example we only have one anyway: d18O)
types = r.types;

% get only Species A single specimens (types value of 0)
depths = depths(types == 0);
carrierA = carrierA(types == 0);

% load the downcore discrete-depth carrier signal
Adisccarmean = p.Adisccarmean;
Adisccarmean = Adisccarmean(:,1); % get the first carrierA signal
discdepth = p.discdepth;

% ---- (7b): Now for the plotting
figure(2)
clf
% plot a cloud of single specimen values
plot(depths,carrierA,'b.','markersize',1,'color',[0.8 0.8 0.8])
% plot the median age on top
hold on
plot(discdepth,Adisccarmean,'k-','linewidth',2);
% make the figure nicer
ylabel('\delta^1^8O')
xlabel('Depth (cm)')
grid on


%% (8) A histogram of single specimens from a particular depth interval

% ---- (8a): Get some stuff to plot --------------------

% tell matlab where the sediment simulation file is
r = matfile(runfile);

% which depth interval to do? (interval >=depth1 and <depth2)
depth1 = 100;
depth2 = 101;

% Load single specimen data from the simulation
depths = r.depths;
ages = r.ages;
carrierA = r.carrierA;

% get only the depth interval and species type we are interested in
ind = find(types == 0 & depths >= depth1 & depths < depth2);
ages = ages(ind);
carrierA = carrierA(ind,1); % the first carrierA signal

% ---- (8b): Now for the plotting. Let's make two histograms, one for age and one for d18O

figure(3)
clf

subplot(1,2,1) % subplot 1
histogram(ages/1000)
xlabel('ka')
ylabel('n specimens')

subplot(1,2,2) % subplot 2
histogram(carrierA)
xlabel('\delta^1^8O')
ylabel('n specimens')
set(gca,'xdir','reverse')


%% (9) A calibration plot (e.g. similar to Fig. 3 in the manuscript)
% This figure will show histograms of 14C age and age from the sediment simulation for
% for a particular discrete depth, along with their expected laboratory 14C age determination
% and resulting calibrated 14C age distribution.


% ---- (9a): First tell matlab where the sediment simulation and picking simulation files are
r = matfile(runfile);
p = matfile(pickfile);

% ---- (9b): Now specify some settings for the plot
age2do = 12800; % Find and plot the discrete depth slice that has true median age closest to this value (years)
curve2use = 'Marine13'; % cal curve to use for calibration
resage2use = 0; % reservoir age to apply during calibration
reserr2use = 0; % reservoir age error to apply during calibration
binwidth = 100; % bin width for histograms (years or 14C years)
manxlim = []; % manual x limits on plot (leave empty for automatic, or e.g. [3 5] for 3 ka to 5 ka)
manylim = []; % manual y limits on plot (leave empty for automatic, or e.g. [3 5] for 3 ka to 5 ka)
% ------

% Now run the code block to produce the figure

figure(4)
clf

[~, di] = min(abs(p.Adiscagemed-age2do)); % finds closest picked sample to desired age2do
sliceint = mean(diff(p.discdepth)); % find the slice interval of the core
depthinterval = [p.discdepth(di,1)-sliceint/2 p.discdepth(di,1)+sliceint/2]; % depth interval of the sample
forind = find(r.depths >= depthinterval(1) & r.depths < depthinterval(2) & r.cycles <= p.Adiscwhole(di,1)); % find all whole forams in the depth interval
AMSage = p.AdiscAMSage(di,1) - resage2use;
AMSerr = (p.AdiscAMSerr(di,1)^2 + reserr2use^2)^0.5;
thisAMSpdf = normpdf([AMSage-3.5*AMSerr:AMSage+3.5*AMSerr], AMSage, AMSerr); % Gaussian distribution of AMS age + AMS error with resage and reserr
[p95_4, p68_2, calprob, medage] = matcal(p.AdiscAMSage(di,1), p.AdiscAMSerr(di,1), 'marine13', 'calbp', 'resage', resage2use, 'reserr', reserr2use, 'plot', 0); % take AMS age and AMS error and make calibrated age distribution using MatCal + Marine13
dgreen = [0 158 115]/255;
lorange = [219 159 114]/255;

% plot true/calibrated ages on X axis
ages = r.ages;
ages = ages(forind);
axcalage = axes;
binedges = min(ages):binwidth:max(ages);
hagehist = histogram(ages/1000,binedges/1000, 'facecolor', dgreen, 'edgecolor', 'none'); % histogram of true ages from SEAMUS_run in green
hold on
scaling = max(hagehist.Values) / max(calprob(:,2));
area(calprob(calprob(:,2)>0.001*max(calprob(:,2)),1)/1000, calprob(calprob(:,2)>0.001*max(calprob(:,2)),2)*scaling,'facecolor', lorange, 'edgecolor', lorange) % calibrated 14C age distribution from MatCal in orange
uistack(hagehist,'top')
X = hagehist.BinEdges;
xlim([min(X)-(max(X)-min(X)) max(X)])
ylim([0 max(hagehist.Values)*4]);
plot(median(ages)/1000,0.28*max(ylim), 'bd', 'markerfacecolor',dgreen, 'markeredgecolor',dgreen) % plot median true age
plot(medage/1000,0.28*max(ylim), 'bs', 'markerfacecolor',lorange, 'markeredgecolor',lorange) % plot median calibrated age
truemedian = round(median(ages));
calmedian = round(medage);
agesig2(1) = prctile(ages,100*(1-erf(2/sqrt(2)))/2);
agesig2(2) = prctile(ages,100-100*(1-erf(2/sqrt(2)))/2);

set(gca,'yticklabel',[])
set(gca,'ytick',[])
set(gca,'xticklabel',[])
set(gca,'xtick',[])
set(gca, 'color', 'none')
if isempty(manxlim) == 0
	xlim(manxlim)
end
xlims = xlim;
if xlims(1) < 0;
	xlims(1) = 0;
	xlim(xlims);
end
axis off

% plot 14C ages in Y axis
plotscale = 1000;
foram14c = r.foram14c;
foram14c = foram14c(forind);
ax14cams = axes;
binedges = min(foram14c):binwidth:max(foram14c);
h14chist = histogram(foram14c/plotscale,binedges/plotscale, 'facecolor', dgreen, 'edgecolor', 'none','orientation','horizontal');
hold on
scaling = max(h14chist.Values) / max(thisAMSpdf);
hamsn = patch(thisAMSpdf*scaling, [AMSage-3.5*AMSerr:AMSage+3.5*AMSerr]/1000,lorange,'edgecolor',lorange);
X = h14chist.BinEdges;
ylim([min(X)-(max(X)-min(X)) max(X)])
xlim([0 max(h14chist.Values)*4])
plot(0.28*max(xlim), mean(foram14c)/1000,'bd', 'markerfacecolor',dgreen, 'markeredgecolor',dgreen) % plot true mean 14C age
plot(0.28*max(xlim), AMSage/1000, 'bs', 'markerfacecolor',lorange, 'markeredgecolor',lorange) % plot AMS mean age
set(gca,'yticklabel',[])
set(gca,'ytick',[])
set(gca,'xticklabel',[])
set(gca,'xtick',[])
set(gca, 'color', 'none')
uistack(h14chist,'top')
if isempty(manylim) == 0
	ylim(manylim)
end
ylims = ylim;
if ylims(1) < 0;
	ylims(1) = 0;
	ylim(ylims);
end
axis off

% plot cal curve
axintcal = axes;
File = fopen(['private/',curve2use,'.14c']);
Contents = textscan(File,'%f %f %f %f %f','headerlines',11,'delimiter',',');
fclose(File);
curvecal = flipud(Contents{1});
curve14c = flipud(Contents{2});
curve14cerr = flipud(Contents{3});
xdata = curvecal/1000;
ydata = curve14c/1000;
onesig = curve14cerr/1000;
fill([xdata' fliplr(xdata')],[ydata'+2*onesig' fliplr(ydata'-2*onesig')],[0.8 0.8 0.8],'edgecolor','none'); % 2 sigma range
hold on
fill([xdata' fliplr(xdata')],[ydata'+onesig' fliplr(ydata'-onesig')],[0.6 0.6 0.6],'edgecolor','none'); % 1 sigma range
xlim(xlims);
ylim(ylims); % match to xlims and ylims of other plots
xlabel('Cal ka BP')
ylabel('^1^4C age (^1^4C ka BP)')
set(gca, 'color', 'none')

% text labels
tboxstr = ['Discrete depth interval: ', num2str(round(depthinterval(1))),'-',num2str(round(depthinterval(2))), ' cm'];
tboxf2 = annotation('textbox',get(gca,'position'),'String',tboxstr);
set(tboxf2, 'linestyle','none')
set(tboxf2, 'horizontalalignment','left')

% clean up order of plotting
axes(ax14cams)
uistack(h14chist,'top')
ylim(ylims)
axes(axcalage)
xlim(xlims)
uistack(hagehist,'top')
uistack(axintcal, 'top')
