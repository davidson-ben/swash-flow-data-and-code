function [frame_time,time_cyc] = ensemble_avg_frames
%==========================================================================
% ensemble_avg_frames.m
%
% Use shoreline motion and bore collapse location to define
% swash zone.  Assign each frame a time on the swash cycle time, 0:2
% (seconds).
% 
% Author: B. Davidson
% Last Updated: 3 October 2025
%==========================================================================

    
    fPath = "../Results";


    load("../1_camera_preprocessing/A00_SWL.mat")
    xref = X_swl;


    
    load("../Results/shoreline.mat",'bins','wave') %shoreline location

    dt = 0.0330;
    last_frame = size(wave,1);
    frame_time = 0:dt:(last_frame-1)*dt;

%% cycles
    [Mxs, peaks] = ea_cycles(wave,bins,xref,0);

%% Remap time
    step = 60; %With T = 2 second waves, we will force the cycles into 60 'sample' groups for ensemble averageing
    
    cyc = length(peaks(:,1))-2; %number of individual cycles - space between each peak frame
    time_cyc = nan(size(frame_time));
    
    for c = 1:cyc %loop through each cycle
        draw = 1;
        if draw == 1
            figure(2)
            clf
            hold on
            cycInq = peaks(c,1):peaks(c+2,1); %2 cycle inquery frames to look at
            plot(cycInq,Mxs(cycInq)) %plot mean shoreline over the inquery frames
            plot(peaks(c:c+2,1),peaks(c:c+2,2),'ro') %plot the mean shoreline over this region
            plot(peaks(c:c+2,1),ones(3,1)*xref,'b--') %plot bore collapse location
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find start of cycle
        x_start = Mxs(peaks(c,1):peaks(c,1)+30); %start of cycle shoreline
    
        if max(x_start)<xref %if shoreline starts after bcl
            disp('Start is after bcl')
            if abs(max(x_start)-xref) <=50 %if start is 'close enough' to bcl ~50px == ~5cm
                disp('Start is close enough to bcl')
                start = peaks(c,1); %just make the first frame of the cycle the start of the cycle
            else
                error('Start is too far from bcl')
            end
        else %shoreline starts before bcl (prefered)
            f_start = peaks(c,1):peaks(c,1)+30; %frames corresponding to start of cycle shoreline
            [uniqueX, uniqueIdX] = unique(x_start, 'first'); %remove duplicates - keep first occurance
            uniqueF = f_start(uniqueIdX); %corresponding x values
            start = interp1(uniqueX,uniqueF,xref); %time [in frames] of shoreline crossing bcl - not necessarily integer!
        end
    
        if draw == 1
            plot(f_start,x_start,'k-','linewidth',1) %plot start of cycle region
            plot(start,xref,'m*') %mark bcl crossing location - start
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find end of cycle
        x_end = Mxs(peaks(c+1,1):peaks(c+1,1)+30); %end of cycle shoreline
    
        if max(x_end)<xref %if shoreline starts after bcl
            disp('End is after bcl')
            if abs(max(x_end)-xref)<=50 %if end is 'close enough' to bcl ~50px == ~5cm
                disp('End is close enough to bcl')
                ed = peaks(c+1,1); %just make the last frame of the cycle the start of the next
            else
                error('End if too far from bcl')
            end
        else %shoreline starts before bcl (prefered)
            f_end = peaks(c+1,1):peaks(c+1,1)+30; %frames corresponding to start of cycle shoreline
            [uniqueX, uniqueIdX] = unique(x_end, 'first'); %remove duplicates - keep first occurance
            uniqueF = f_end(uniqueIdX); %corresponding x values
            ed = interp1(uniqueX,uniqueF,xref); %time [in frames] of shoreline crossing bcl - not necessarily integer!
        end
    
        if draw == 1
            plot(f_end,x_end,'k-','linewidth',1) %plot end of cycle region
            plot(ed,xref,'m*') %mark bcl crossing location - end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        cyc_dur(c) = ed-start; %cycle duration [frames] - likely not integer
    
        steps = linspace(start,ed,61); % time steps (frames) spaced out with 61 steps from start of current cycle to end of current cycle.  The end of the current cycle is the start of the next cycle, so we truncate off the last step as that will be the next starting step (frames)
        t_nd = linspace(0,2,61); %dimensionless time for the cycle
        xq = round(steps(1)):round(steps(end));
        %Give Each frame a new time stamp
        time_cyc(xq) = interp1(steps,t_nd,xq);
    end


end

