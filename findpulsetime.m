function [ startpulseindex, endpulseindex] = findpulsetime( currenttrace)

%findpulsetime takes a current trace that contains a square pulse and
% outputs the start and end times of the pulse.

        derivpulse = diff(currenttrace);
        [minimum,startpulseindex] = min(derivpulse);
        [maximum,endpulseindex] = max(derivpulse);
        if startpulseindex <= endpulseindex
            ;
        else
            a = startpulseindex;
            b = endpulseindex;
            startpulseindex = b;
            endpulseindex = a;
        end
  
end

