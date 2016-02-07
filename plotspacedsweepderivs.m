function [samplesweep, sampling_interval] = plotspacedsweepderivs(numfiles,numfilesthreshold)
% plotspacedsweeps takes a number of files, as well as a desired number of 
% plots, generates evenly spaced figures.  Also outputs the sampling
% interval;

filedir = pwd;
files = dir('*.abf');

if numfilesthreshold >= numfiles
        numfilesthreshold = numfiles;
        files2plot = [1:1:numfiles];
        else files2plot = ceil(linspace(1,numfiles,numfilesthreshold));
end

for m = 1:numfilesthreshold
        plottingfile = files2plot(m);
        [data,sampling_interval]=abfload(strcat(filedir,'/',files(plottingfile).name));
        sweeps(m)=size(data,3);
        figure;
        samplesweep=(data(:,:,sweeps(m)));
        samplederiv = samplesweep(2:end) - samplesweep(1:end-1); % unscaled derivative
        plot(samplederiv);
end
    

end

