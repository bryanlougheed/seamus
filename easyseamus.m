function [depthout, proxyout] = easyseamus(proxyin, timein, agedepth, bd, abu)
%[depthout, proxyout] = easyseamus(proxyin, timein, agedepth, bd, abu)
%
% EasySEAMUS is like SEAMUS, but carried out probabilistically on the mean signal, 
% i.e. without explicitly simulating single forams. 
% Much less resource intensive than SEAMUS, but also less fun?
%
% Output has been validated against both TURBO2 output and mean signal output of
% a SEAMUS single foram simulation.
%
% Garbage in (time domain)
% ------------------------
% proxyin  = vector of proxy data
% timein   = age of proxy data in ka, same dims as proxyin
% agedepth = age-depth relationship of the core archive 
%            Two columns. Col 1 = age in ka, Col 2 = depth in cm
% bd  = Bioturbation depths in time domain. Two columns. Col 1 = age in ka, Col 2 = BD in cm
% abu = Abundance in time domain. Two columns. Col 1 = age in ka, Col 2 = abundance between 0 and 1.
%
% Garbage out (core depth domain)
% -------------------------------
% depthout = series of central depth values of discrete 1 cm depths
% proxyout = the mean value of the proxy data associated with depthout
%
% B.C. Lougheed, September 2020
% bryan.lougheed@geo.uu.se

% make in yr instead of ka
timein = timein*1000;
agedepth(:,1) = agedepth(:,1) * 1000;
bd(:,1) = bd(:,1) * 1000;
abu(:,1) = abu(:,1) * 1000;

% sort everything young to old
abu = sortrows(abu);
bd = sortrows(bd);
agedepth = sortrows(agedepth);
[timein, ind] = sort(timein);
proxyin = proxyin(ind);

% make inputs annual resolution
proxya = interp1(timein,proxyin,min(timein):1:max(timein));
timea = min(timein):1:max(timein);
agedeptha(:,1) = timea;
agedeptha(:,2) = interp1(agedepth(:,1),agedepth(:,2),timea,'linear','extrap');
sar = diff(agedeptha(:,2)) ./ diff(agedeptha(:,1));
sar = [sar(1); sar];
bda(:,1) = timea;
bda(:,2) = interp1(bd(:,1),bd(:,2),timea);
bda(bda(:,1)<=bd(1,1),2) = bd(1,2); % extend with first value
bda(bda(:,1)>=bd(end,1),2) = bd(end,2); % extend with final value
bda(bda(:,2)<0.000001,2) = 0.000001; % represent zero BD as super small
abua(:,1) = timea;
abua(:,2) = interp1(abu(:,1),abu(:,2),timea);
abua(abua(:,1)<=abu(1,1),2) = abu(1,2);
abua(abua(:,1)>=abu(end,1),2) = abu(end,2);

% do bioturbation simulation
depthints = round(min(agedepth(:,2))):1:round(max(agedepth(:,2)));
bddepth = interp1(agedeptha(:,2),bda(:,2),depthints);
sardepth = interp1(agedeptha(:,2),sar,depthints);
timedepth = (depthints-bddepth) ./ sardepth; % -bddepth because adjusted to focal point for bioturbation
% extend proxy and abundance data for the deepest depth (where bioturbation can go beyond the end of the record)
ppri = exp(-((0:round(bddepth(end)*6/sardepth(end))).*sardepth(end))./bddepth(end));
proxya(end:end+numel(ppri)) = proxya(end); 
abua(end:end+numel(ppri),1) = NaN; 
abua(end:end+numel(ppri),2) = abua(end); 
% the discrete depths
depthout = NaN(numel(depthints)-1,1);
proxyout = NaN(size(depthout));
for i = 1:numel(depthints)-1
	depthout(i) = (depthints(i) + depthints(i+1)) / 2;
	[~, ind] = min(abs(timedepth(i)-timea));
	ind = ind(1); % om det finns två lika
	% prior age distribution of this discrete depth following Berger and Heath (1968). 5 bioturbation depths should be enough.
	ppri = exp(-((0:round(bddepth(i)*5/sardepth(i))).*sardepth(i))./bddepth(i));
	ppri = ppri/sum(ppri); % normalise
	ppri = ppri .* abua(ind:ind+numel(ppri)-1,2)'; % also incorporate abundance
	ppri = ppri/sum(ppri); % normalise
	proxyout(i) = sum(ppri .* proxya(ind:ind+numel(ppri)-1),2);  % expected proxy signal value for this discrete depth
end
