%% Welcome to the SEAMUS tutorial, a step-by-step hands-on guide to the simulation and plotting of results
% This particular tutorial is designed for the Octave environment.

pkg load statistics % load the Octave statistics package

% (1): Let's load some kind of interesting climate signal
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


% (2): Let's prepare the input variables for seamus_run sediment simulation
% (2a): General simulation settings for seamus_run --------
simstart = 70000;          % start year (yrs BP)
siminc = 10;               % simulation timestep (yrs)
simend = 0;                % end year (yrs BP)
btinc = 10;                % bioturbation timestep (yrs).
realD14C = 1;              % Assign 'real' D14C to single specimens (1 = yes, 0 = no)
calcurve = 'Marine13';     % Calibration curve to use for assigning 'real' D14C
blankbg = 46806;           % the 14C blank value to assign to single specimens (46806 = higest value in Marine13 = practical blank)
fpcm = 10^2;               % Core capacity in specimens per cm. Let's use 10^2 for now to run quickly (you can change later). (Plotting large amounts of data points difficult in Octave, though)
runfile = 'test_run.mat';  % Name of output file from seamus_run

% (2b): Dynamic inputs for seamus_run --------------------

% bioturbation depths   ka   BD (cm)
bdpoints =				[0    10
						50   10
						100  10];
% The above will result in a constant bioturbation depth
% of 10 cm. You can enter dynamic values by adding extra rows
% to the matrix and the simulation will interpolate/extrapolate
% where needed. Or you can simply load data from external file.


% age-depth points      ka  depth (cm)
adpoints =				[0   0
						25   250
						50   500];
% The above will result in constant SAR of 10 cm/ka.
% You can enter as many age depth points as you like, and
% the simulation will interpolate/extrapolate where needed.
% You could also load an external file instead.

% Species A carrier signal(s)
%                       ka        signal
carrierA =              [ngripage ngrip18O];
% Here we simply use the NGRIP data that we prepared earlier.
% For carrier signals, the simulation will interpolate, but
% not extrapolate.

% Species A abundance   ka     fraction of specimen flux
speciesA =				[0      1
						25     1
						50     1];

            
            
% (3): Run the seamus_run sediment simulation (can take some time)
tic
seamus_run(simstart, siminc, simend, btinc,	fpcm, realD14C,...
	blankbg, adpoints, bdpoints, runfile,...
	'calcurve', calcurve,'carrierA',carrierA,'speciesA',speciesA);
toc
% you can add more optional inputs, see documentation
% For full options type 'help seamus_run' into command window or consult
% supplementary tables in manuscript.



% (4): Now prepare the input variables for the seamus_pick picking and 14C dating simulaiton
% (4a): General simulation settings for seamus_pick
pickint = 1;                  % thickness of the core slices (cm)
Apickfordate = -1;            % number of Species A specimens to pick per sample (enter -1 for all available)
Bpickfordate = -1;            % number of Species A specimens to pick per sample (enter -1 for all available)
calcurve = 'Marine13';        % calibration curve to use when calibrating dates
pickfile = 'test_pick.mat';   % file to save the picking results to
matfilein = runfile;          % The name of the sediment simulation to analyse (the file saved by seamus_run)



% (5): Now run the seamus_pick picking and 14C dating simulation (can take some time)
tic
seamus_pick(matfilein, pickfile, calcurve, pickint, Apickfordate, Bpickfordate)
toc
% you can add more optional inputs, see documentation
% For full options type 'help seamus_pick' into command window or consult
% supplementary tables in manuscript.



% (6) Now we can analyse the output data and make some figures.
% Let's start with a simple downcore figure showing depth vs age
% ---- (6a): Get some stuff to plot --------------------
% Load single specimen data to plot:
load(runfile, 'depths', 'ages', 'types');
% get only Species A single specimens (types value of 0)
depths = depths(types == 0);
ages = ages(types == 0);
% load the downcore discrete-depth median age from the pickfile
load(pickfile, 'Adiscagemed', 'discdepth');
% ---- (6b): Now for the plotting
figure(1)
clf
% plot a cloud of single specimen values
plot(ages/1000,depths,'b.','markersize',1,'color',[0.8 0.8 0.8]) % divide age by 1000 to make ka (easier for the plot labels)
% plot the discrete-depth median ages on top
hold on
plot(Adiscagemed/1000,discdepth,'k-','linewidth',2);
% make the figure nicer
set(gca,'ydir','reverse')
xlabel('ka')
ylabel('Depth (cm)')
grid on



% (7) Figure of carrier signal vs simulated core depth
% ---- (7a): Get some stuff to plot --------------------
% Load single specimen data to plot:
load(runfile, 'depths','carrierA','types')
carrierA = carrierA(:,1); % get the first carrierA signal for Species A (in this example we only have one anyway: d18O)
% get only Species A single specimens (types value of 0)
depths = depths(types == 0);
carrierA = carrierA(types == 0);
% load the downcore discrete-depth carrier signal
load(pickfile,'Adisccarmean','discdepth')
Adisccarmean = Adisccarmean(:,1); % get the first carrierA signal
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


% (8) A histogram of single specimens from a particular depth interval
% ---- (8a): Get some stuff to plot --------------------
% tell matlab where the sediment simulation file is
% which depth interval to do? (interval >=depth1 and <depth2)
depth1 = 100;
depth2 = 101;
% Load single specimen data from the simulation
load(runfile,'depths','ages','carrierA')
% get only the depth interval and species type we are interested in
ind = find(types == 0 & depths >= depth1 & depths < depth2);
ages = ages(ind);
carrierA = carrierA(ind,1); % the first carrierA signal
% ---- (8b): Now for the plotting. Let's make two histograms, one for age and one for d18O
figure(3)
clf

subplot(1,2,1) % subplot 1
hist(ages/1000) % divide by 1000 to get ka
xlabel('ka')
ylabel('n specimens')

subplot(1,2,2) % subplot 2
hist(carrierA)
xlabel('\delta^1^8O')
ylabel('n specimens')
set(gca,'xdir','reverse')



%% (9) A calibration plot (e.g. similar to Fig. 3 in the manuscript)
% The figure is too complex for the Octave plotting capabilities.
% Perhaps plot in external plotting software?
% See the Matlab tutorial file for an idea of how it is made.