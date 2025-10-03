function [Mxs, peaks_output] = ea_cycles(wave,bins,bcl,drawPlot)
%==========================================================================
% ea_cycles.m
%
% Plot the shoreline and find the start of each swash cycle
% Use the bore collapse location as the start of the swash cycle
% 
% Author: B. Davidson
% Last Updated: 3 October 2025
%==========================================================================

    for i = 1:size(wave,1)
        p(i,:) = polyfit(bins,wave(i,:),4); %polynomial fit of the shoreline
    end

    % Define explicit swash cycles
    Mxs = medfilt1(mean(wave,2),6); %mean shoreline
    frames = 1:length(Mxs); %all frames
    %We don't care about drawdown (since it is messy) so I want to define the
    %start of the swash by the most offshore location that we see the wave:
    peaks = islocalmax(Mxs); %find local max
    %we only want the points at the top, so lets threshold by mean + x * sd
    m = mean(Mxs);
    sd = std(Mxs);
    thresh = m+1*sd;
    peaks(Mxs<thresh) = 0;
    pFrames = frames(peaks); %
    pFramesLocs = Mxs(find(peaks));
    peaks_output = [pFrames' pFramesLocs];

    if drawPlot == 1
        plot(frames,Mxs) %plot mean shoreline as a function of frame number
        xlabel('time [frames]')
        ylabel('X (positive is offshore) [px]')
        set(gca,'FontSize',15)
        hold on
        plot([min(frames) max(frames)],[thresh thresh]) %peak threshold
        plot(pFrames,pFramesLocs,'ro') %maximum locations - peakFrames
        plot(frames,ones(size(frames))*bcl,'b--') %Bore Collapse Location

        legend({'mean shoreline position','peak threshold','Peak Locations','Bore Collapse Location'},'location','south')
    end
end