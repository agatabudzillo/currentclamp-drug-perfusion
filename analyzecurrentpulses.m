% Agata Budzillo

% August 12, 2014
% Setup: a folder with at least one file containing IV square pulses. All the files in the folder 
% should be from a single cell.  Files ideally contain both a current and voltage trace. 
% Parameters in abfload function may need to be modified depending on the user's Clampex channel naming 
% conventions. Note that if there is no current trace the user will need to hard code the start and end 
% times for the current pulse, "startpulseindex" & "endpulseindex". If you do this you may also wish to 
% eliminate the safety offset I introduce to make sure I don't capture stimulation artifact in my 
% measurements. You should just be able to set "allowanceforpulse" to 0. Note that input resistance 
% calculation assumes current is in pA and voltage is in mV.  If you want
% to suppress the spike finding warnings, see comment below

% Associated functions:
% abfload(v2),agatasfindpeaks,findpulsetime,getcellconditions,getspikeintervals,plotspacedsweeps,pvpmod

% Modified DJP and CAW 3 Sept 2014
% Manually give pulse start and end times and amplitudes for Andy's data


clear; close all;
warning ('off','resoftwareonserver:agatasfindpeaks:noPeaks');
warning ('off','resoftwareonserver:agatasfindpeaks:largeMinPeakHeight');
% warnings will have to be updated depending on your current directory.  After you run the code and a warning
% w = warning('query','last') will give you the warning info and w.identifier will print the identifier.

%% CONSTANTS

% Data format
currentChannelRecorded = 1; % If current was recorded, this should be 1 (Agata's case)

% FILE RELATED
if currentChannelRecorded == 1
    voltagechannelname = 'IN 0 CC 1';
    currentchannelname = 'IN 1 CC 1';
else
    voltagechannelname = 'IN 0';
    %currentchannelname = 'IN 1 CC 1';
end
abffilenamelength=8; % (filename size in Clampex)

% MEASUREMENT RELATED
numfilesthreshold = 3; % the number of files to use when displaying sweep used to determine spike threshold
minpeakdistancems = 2.5; % minimum time in ms between spikes
testcurrentsweep=1; % the current sweep# used for finding pulse times
allowanceforpulsems = 5; % time after peak of current pulse in ms at which to begin analysis
sagwindow = 100; % in ms
steadystatewindow = 500; % in ms
maxcurrentforIVfit = 1; % maximum current step pulse used for IV slope fit
polynomialdegree = 1; % for linear IV slope fit

% GRAPH RELATED
maxplottableFR = 20; % sweeps with a FR above this will be ignored for image clarity when plotting example overlayed sweeps
voltagemin = -120; % normalized yaxis for plotting overlayed sweeps
voltagemax = 30;
%% ACCESS FILES TO BE ANALYZED, GET EXPERIMENTAL CELL CONDITIONS, SET UP DATA STRUCT, GET Fs, FIND SPIKE SORTING THRESHOLD

filedir = uigetdir;
cd(filedir); 

files = dir('*.abf');
numfiles = size(files,1);
[conditions,numconditions] = getcellconditions(numfiles);
colormap_conditions = hsv(numconditions); % check out other pretty colormaps!
conditionindices=conditions(1:numconditions,2);
      
for f = 1:numconditions
    mkdir(conditions{f,1});
    conditionstarts(f) = conditionindices{f};
end
mkdir('plots');

datastruct = struct('filename',[],'condition',[],'sampling_frequency',[],...
    'filestartsec',[],'time',[],'voltage',[],'restingpotential',[],...
    'inputresistance',[],'spikes',{},'spikeintervals',{},'firingrate',[],...
    'firingrateNANfree',[], 'avgFiringRate', [], 'earlyresponseMIN',{},'earlyresponseAVG',{},...
    'lateresponseAVG',{},'IVslope',[],'IVintercept',[],'currentpulse',[],...
    'currentstart',[],'currentend',[],'currentsize',[]);

[samplesweep,sampling_interval] = plotspacedsweeps(numfiles,numfilesthreshold);
[samplederiv,sampling_interval] = plotspacedsweepderivs(numfiles,numfilesthreshold);

useVThresh = input('Enter 1 to use voltage threshold or anything else for derivative');
if (useVThresh == 1)
    minpeakheight = input('Enter Spike Peak Threshold: '); % spike sorting threshold
else
    derivThresh = input('Enter Spike Slope Threshold: '); % spike deriv threshold
end

if (currentChannelRecorded == 0)
    startpulseindex = input('Enter start index for current pulse: ');
    endpulseindex = input('Enter end index for current pulse: ');
    firstCurrentAmplitude = input('Enter amplitude of first current pulse (pA): ');
    currentAmplitudeIncrement = input('Enter increment of current amplitude (pA): ');
end

close all;

%% CONVERSIONS AND SET UP OF PARAMETERS

nA_pAconversionfactor = 1000;
us_sconversionfactor = 1E6;
s_msconversionfactor = 1000;
ms_sconversionfactor = 1/s_msconversionfactor;
gigaohm_megaohmconversionfactor = 1E3; % for calculating input resistance from IV slope

sampling_interval = sampling_interval/(us_sconversionfactor); % now in s
sampling_frequency = 1/sampling_interval;
time = (0:sampling_interval:sampling_interval*(length(samplesweep)-1))*s_msconversionfactor;

minpeakdistance = minpeakdistancems * ms_sconversionfactor * sampling_frequency;
allowanceforpulse = allowanceforpulsems * ms_sconversionfactor * sampling_frequency;

sagwindow_samples = sagwindow * ms_sconversionfactor * sampling_frequency; 
steadystatewindow_samples = steadystatewindow * ms_sconversionfactor * sampling_frequency;

colorcondition = 1;
currentcondition = [];

% INITIALIZE FIGURES WHICH WILL BE ACCESSED AFTER ANALYSIS OF EACH FILE
vrest = figure('Visible','off'); rinput = figure('Visible','off'); earlyiv = figure('Visible','off');
lateiv = figure('Visible','off'); fi = figure('Visible','off'); sag = figure('Visible','off');
rinput = figure('Visible','off');
saturation_filenum= linspace(0,0.9,numfiles);

%% READ AND ANALYZE EACH FILE IN THE FOLDER
for n = 1:numfiles
    [data,si,header]=abfload(strcat(filedir,'/',files(n).name),'channels', {voltagechannelname}, 'sweeps','a');
    
    % If Agata's file format read current from trace
    if (currentChannelRecorded == 1)
        current = abfload(strcat(filedir,'/',files(n).name),'channels', {currentchannelname}, 'sweeps','a');
    end
    numsweeps = size(data,3);
    filestartsec = header.lFileStartTime;
    filenumber = files(n).name(1:abffilenamelength);
    endfilenumber = str2num(filenumber((end-2):end));
    
    % REFER BACK TO USER DEFINED EXPERIMENTAL CONDITIONS TO DEFINE CURRENT FILE
    % find highest condition start time that this file is equal or greater than
    colorcondition = max(find(endfilenumber >= conditionstarts));
    currentcondition = conditions{colorcondition,1};
    
    % RESET ARRAYS FOR STORING SWEEP SPECIFIC INFORMATION
    currentsize = zeros(numsweeps,1);
    startpulsems = [];
    endpulsems = [];
    restingpotential = zeros(numsweeps,1);
    spikes = cell(numsweeps,[]);
    spikeintervals = cell(numsweeps,[]);
    firingrate = zeros(numsweeps,1);
    avgFiringRate = zeros(numsweeps, 1);
    earlyresponseMIN = zeros(numsweeps,1);
    earlyresponseAVG = zeros(numsweeps,1);
    lateresponseAVG = zeros(numsweeps,1);
    
    % FIND CURRENT PULSE START AND END TIMES from the 1st sweep's CURRENT trace in the .abf file 
    if (currentChannelRecorded == 1)
        [startpulseindex, endpulseindex] = findpulsetime(current(:,:,testcurrentsweep));  
    end
    
    startpulsems = (startpulseindex*sampling_interval)*s_msconversionfactor;
    endpulsems = (endpulseindex*sampling_interval)*s_msconversionfactor;
    pulselength = endpulsems-startpulsems;
        
    % INITIALIZE FIGURE (WILL BE ACCESSED AFTER ANALYSIS OF EACH SWEEP)
    overlaysweeps=figure('visible','off');
    saturation_sweepnum= gray(numsweeps); 
 
    
    %% READ AND ANALYZE EACH SWEEP IN THE CURRENT FILE 
    for p = 1:numsweeps     
        % DETERMINE THE PARAMETERS OF THE CURRENT INJECTION, AND ISOLATE THE DESIRED PARTS OF THE RESPONSE
        if (currentChannelRecorded == 1)
            baselinecurrent = mean(current(1:(startpulseindex-allowanceforpulse),:,p));
            pulsecurrent = mean(current((startpulseindex+allowanceforpulse):(endpulseindex-allowanceforpulse),:,p));
            currentsize(p) = nA_pAconversionfactor*(pulsecurrent-baselinecurrent);
        else
            currentsize(p) = firstCurrentAmplitude + (p-1)*currentAmplitudeIncrement;
        end
        
        % INTRODUCE A SHORT TIME OFFSET TO MAKE SURE THERE IS NO STIMULATION ARTIFACT IN MEASUREMENTS
        prestartpulse = startpulseindex - allowanceforpulse;
        poststartpulse = startpulseindex + allowanceforpulse;
        preendpulse = endpulseindex - allowanceforpulse;

        restingpotential(p) = mean(data(1:prestartpulse,:,p));
        earlyresponse = data(poststartpulse:(poststartpulse + sagwindow_samples),p);
        lateresponse = data((preendpulse - steadystatewindow_samples):preendpulse,:,p);
        wholeresponse = data(poststartpulse:preendpulse,:,p);
        
        % FIND SPIKES, CALCULATE INTERSPIKE INTERVALS
        if (useVThresh == 1)
            [pks,locs]=agatasfindpeaks(wholeresponse,'minpeakheight',minpeakheight,'minpeakdistance',minpeakdistance);
        else
            deriv = wholeresponse(2:end) - wholeresponse(1:end-1);
            [pks, locs] = findpeaks(deriv, 'minpeakheight', derivThresh, 'minpeakdistance', minpeakdistance);
        end
        %[pks,locs]=findpeaks(wholeresponse,'minpeakheight',minpeakheight,'minpeakdistance',minpeakdistance);
        locsms = (locs*sampling_interval)*s_msconversionfactor + allowanceforpulsems;
        %[spike_intervals] = getspikeintervals(locs);
        spike_intervals = locs(2:end) - locs(1:end-1); %% DJP added
        spikeintervalsms = (spike_intervals*sampling_interval)*s_msconversionfactor;
        
        % FIND SUBTHRESHOLD RESPONSES - If there are not spikes - measure voltage.  "0" = no voltage measurement
        if (useVThresh ~= 1)
            minpeakheight = derivThresh;
        end
        [pks2,locs2] = agatasfindpeaks(earlyresponse,'minpeakheight',minpeakheight,'minpeakdistance',minpeakdistance);
        [pks3,locs3] = agatasfindpeaks(lateresponse,'minpeakheight',minpeakheight,'minpeakdistance',minpeakdistance);
  
            if isempty(locs2)
                if currentsize(p) < 0
                    [minvals,minindices] = min(earlyresponse);
                    minindices=min(minindices);
                    earlymin = mean(data((poststartpulse+minindices-1):(poststartpulse+minindices+1)));
                    earlyavg = mean(earlyresponse);
                else
                    earlymin = 0;
                    earlyavg = mean(earlyresponse);
                end    
            else
                earlymin = 0;
                earlyavg = 0;
            end 
        
            if isempty(locs3)
                lateavg = mean(lateresponse);
            else
                lateavg = 0;
            end

        % WRAP UP SWEEP SPECIFIC MEASUREMENTS
        firingrate(p) = 1/(mean(spikeintervalsms)*ms_sconversionfactor);
        avgFiringRate(p) = length(locs)/(pulselength*ms_sconversionfactor);
        spikes{p} = locsms;
        spikeintervals{p} = spikeintervalsms;
        earlyresponseMIN(p) = earlymin;
        earlyresponseAVG(p) = earlyavg;
        lateresponseAVG(p) = lateavg;
        
        % UPDATE FIGURE OF OVERLAYED SWEEPS
        figure(overlaysweeps); set(overlaysweeps,'visible','off');

        if firingrate(p) >= maxplottableFR
            ;
        else
            hold on
            plot(time,data(:,:,p),'Color',(saturation_sweepnum(p,:)));
            hold on
        end

end
        % WRAP UP FILE SPECIFIC MEASUREMENTS
        restingpotential_file = mean(restingpotential);
        frONEspikeindex = find(isinf(firingrate));
        frNOspikesindex = find(isnan(firingrate));
        NANfreefiringrate = firingrate;
        NANfreefiringrate(frONEspikeindex) = 0;
        NANfreefiringrate(frNOspikesindex) = 0;
        
        % LABEL, SAVE AS PDF & CLOSE FIGURE WITH OVERLAYED SWEEPS
        figure(overlaysweeps); set(overlaysweeps,'visible','off');
        hold off
        xlabel('Time (sec)');
        ylabel('Membrane Potential (mV)');
        ylim([voltagemin voltagemax])
        titlei=[filenumber,' ','sweeps'];
        filenamei=['sweeps','-',filenumber];
        title(titlei,'Interpreter','none')
        saveas(overlaysweeps,strcat(currentcondition,'/',filenamei),'pdf');
        close(overlaysweeps);
        
        %% UPDATE FIGURES WITH CURRENT FILE'S MEASUREMENTS
        % RESTING POTENTIAL
        figure(vrest);set(vrest,'visible','off');
        hold on
        plot(filestartsec,restingpotential_file,'Color',colormap_conditions(colorcondition,:),'Marker','o');
        allhandles_vrest(colorcondition) = plot(filestartsec,restingpotential_file,'Color',colormap_conditions(colorcondition,:),'Marker','o');
        
        % EARLY PART OF VOLTAGE RESPONSE
        nonzeroindicesearlyiv=find(earlyresponseAVG);
        figure(earlyiv);set(earlyiv,'visible','off');
        hold on
        plot(currentsize(nonzeroindicesearlyiv),earlyresponseAVG(nonzeroindicesearlyiv),'Color',colormap_conditions(colorcondition,:),'Marker','o');
        allhandles_earlyiv(colorcondition) = plot(currentsize(nonzeroindicesearlyiv),earlyresponseAVG(nonzeroindicesearlyiv),'Color',colormap_conditions(colorcondition,:),'Marker','o');

        % LATE PART OF VOLTAGE RESPONSE
        nonzeroindiceslateiv=find(lateresponseAVG);
        figure(lateiv);set(lateiv,'visible','off');
        hold on
        xIV = currentsize(nonzeroindiceslateiv);
        yIV = lateresponseAVG(nonzeroindiceslateiv);
        plot(xIV,yIV,'Color',colormap_conditions(colorcondition,:),'Marker','o');
        allhandles_lateiv(colorcondition) =  plot(currentsize(nonzeroindiceslateiv),lateresponseAVG(nonzeroindiceslateiv),'Color',colormap_conditions(colorcondition,:),'Marker','o');
        
        % CALCULATE & PLOT INPUT RESISTANCE 
        fitindices = find(xIV < maxcurrentforIVfit);
        xIVforfit = currentsize(fitindices);
        yIVforfit = lateresponseAVG(fitindices);
        IVfit = polyfit(xIVforfit,yIVforfit,polynomialdegree);
        IVslope = IVfit(1);
        IVintercept = IVfit(2);
        Rinput_megaOhms = IVslope * gigaohm_megaohmconversionfactor;
        
        figure(rinput);set(rinput,'visible','off');
        hold on
        plot(filestartsec,Rinput_megaOhms,'Color',colormap_conditions(colorcondition,:),'Marker','o');
        allhandles_rinput(colorcondition) = plot(filestartsec,Rinput_megaOhms,'Color',colormap_conditions(colorcondition,:),'Marker','o');
        
        % FIRING RATE
        nonzeroindicesNANfreefiringrate=find(NANfreefiringrate);
        if isempty(nonzeroindicesNANfreefiringrate)
            ;
        else
            figure(fi);set(fi,'visible','off');
            hold on
            plot(currentsize(nonzeroindicesNANfreefiringrate),NANfreefiringrate(nonzeroindicesNANfreefiringrate),'Color',colormap_conditions(colorcondition,:),'Marker','o');
            allhandles_fi(colorcondition) = plot(currentsize(nonzeroindicesNANfreefiringrate),NANfreefiringrate(nonzeroindicesNANfreefiringrate),'Color',colormap_conditions(colorcondition,:),'Marker','o'); 
        end

         if isempty(nonzeroindicesNANfreefiringrate)
            ;
        else
            figure(fi);set(fi,'visible','off');
            hold on
            plot(currentsize(nonzeroindicesNANfreefiringrate),avgFiringRate(nonzeroindicesNANfreefiringrate),'Color',colormap_conditions(colorcondition,:),'Marker','o');
            allhandles_fi(colorcondition) = plot(currentsize(nonzeroindicesNANfreefiringrate),avgFiringRate(nonzeroindicesNANfreefiringrate),'Color',colormap_conditions(colorcondition,:),'Marker','o'); 
         end
         
         % HYPERPOLARIZATION INDUCED SAG
        nonzeroindicesearlyresponseMIN=find(earlyresponseMIN);
        if isempty(nonzeroindicesearlyresponseMIN)
            ;
        else
        figure(sag);set(sag,'visible','off');
        hold on
        plot(currentsize(nonzeroindicesearlyresponseMIN),earlyresponseMIN(nonzeroindicesearlyresponseMIN),'Color',colormap_conditions(colorcondition,:),'Marker','o');
        allhandles_sag(colorcondition) = plot(currentsize(nonzeroindicesearlyresponseMIN),earlyresponseMIN(nonzeroindicesearlyresponseMIN),'Color',colormap_conditions(colorcondition,:),'Marker','o'); 
        end
        
        %% FILL IN STRUCT WITH MEASUREMENTS
        datastruct(n).filename = filenumber;
        datastruct(n).condition = currentcondition;
        datastruct(n).sampling_frequency = sampling_frequency;
        datastruct(n).filestartsec = filestartsec;
        datastruct(n).time = time;
        datastruct(n).voltage = data;
        if (currentChannelRecorded == 1)
            datastruct(n).currentpulse = current;
        else
            datastruct(n).currentpulse = 0; % FIXME: should be an array
        end
        datastruct(n).currentsize = currentsize;
        datastruct(n).currentstart = mean(startpulsems);
        datastruct(n).currentend = mean(endpulsems);
        datastruct(n).restingpotential = restingpotential_file;
        datastruct(n).inputresistance = Rinput_megaOhms;
        datastruct(n).spikes = spikes;
        datastruct(n).spikeintervals = spikeintervals;
        datastruct(n).firingrate = firingrate;
        datastruct(n).avgFiringRate = avgFiringRate;
        datastruct(n).firingrateNANfree = NANfreefiringrate;
        datastruct(n).earlyresponseMIN = earlyresponseMIN;
        datastruct(n).earlyresponseAVG = earlyresponseAVG;
        datastruct(n).lateresponseAVG = lateresponseAVG;
        datastruct(n).IVslope = IVslope;
        datastruct(n).IVintercept = IVintercept;
       
        % MOVE THE .ABF FILE TO ITS APPROPRIATE FOLDER
        %movefile(strcat(filenumber,'.abf'),strcat(currentcondition,'/'));
end

save('measurementsALL','datastruct');

% LABEL FIGURES, SAVE AS PDFs & CLOSE 
figure(vrest)
hold off
xlabel('Time (s)');
ylabel('Resting Membrane Potential (mV)');
legend(allhandles_vrest,conditions{:,1});
titlevrest=['Vrest'];
title(titlevrest,'Interpreter','none')
saveas(vrest,strcat('plots/',titlevrest),'pdf');
%close(vrest);

figure(earlyiv)
hold off
xlabel('Current (pA)');
ylabel('Membrane Potential (mV)');
legend(allhandles_earlyiv,conditions{:,1},'location','NorthWest');
titleearlyiv=['earlyIV'];
title(titleearlyiv,'Interpreter','none')
saveas(earlyiv,strcat('plots/',titleearlyiv),'pdf');
close(earlyiv);

figure(lateiv)
hold off
xlabel('Current (pA)');
ylabel('Membrane Potential (mV)');
legend(allhandles_lateiv,conditions{:,1},'location','NorthWest');
titlelateiv=['lateIV'];
title(titlelateiv,'Interpreter','none')
saveas(lateiv,strcat('plots/',titlelateiv),'pdf');
%close(lateiv);
        
figure(rinput)
hold off
xlabel('Time (s)');
ylabel('Input Resistance (MOhm)');
legend(allhandles_rinput,conditions{:,1});
titlerinput=['rinput'];
title(titlerinput,'Interpreter','none')
saveas(rinput,strcat('plots/',titlerinput),'pdf');
%close(lateiv);

figure(fi)
hold off
xlabel('Current (pA)');
ylabel('Avg Firing Rate (Hz)');
legend(allhandles_fi(find(allhandles_fi)),conditions{find(allhandles_fi),1},'location','NorthWest');
titlefi=['FI'];
title(titlefi,'Interpreter','none')
saveas(fi,strcat('plots/',titlefi),'pdf');
%close(fi);

figure(sag)
hold off
xlabel('Current (pA)');
ylabel('Membrane Potential (mV)');
legend(allhandles_sag(find(allhandles_sag)),conditions{find(allhandles_sag),1});
titlesag=['sag'];
title(titlesag,'Interpreter','none')
saveas(sag,strcat('plots/',titlesag),'pdf');
close(sag);