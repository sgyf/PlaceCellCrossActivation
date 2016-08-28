function BoxTrackSeqTimeIntervalCompare_EvtPeakAngle
% Part of BoxTrackSeqTimeIntervalCompare. This function compute cells' peak firing angle in each
% turn events, time interval between neighbouring peak angles in each turn event.
%
% Created on 8/25/2016, XM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Choose date and session, import SpikeDB

clear; clc;
global ep_time t angles

prompt   = {'Input full folders path (eg, E:\Dan11\Dan11-2013-10-17), delimitered by white space',...
    'Video Frame interval',...
    'Angular bin size',...
    'Smooth parameter-Half kernel length',...
    'Smooth parameter-kernel step',...
    'Smooth parameter-Sigma'};
dlgtitle = 'Parameter to compute cell active place';
nline    = [1,85;1,85;1,85;1,85;1,85;1,85];
default  = {'','1','6','5','1','2'};
answer   = inputdlg(prompt,dlgtitle,nline,default,'on');

if ~isempty(answer)   % If cancel buttom on inputdlg is hit, return
    datefolderlist = strsplit(answer{1});   % string cellarray of dates
else
    disp('-->No Date folder selected, abort execution!');  disp('********************')
    return
end
D = strsplit(datefolderlist{1},'\');

% Map proper spikeDB to current function workspace
apath   = fullfile(D{1},D{2},[D{2},'_Database\']);
[spkdbFile,spkdbPath] = uigetfile([apath,'*.mat'],'Select ONE spikeDB .mat File','MultiSelect','on');
DBobj   = matfile([spkdbPath,spkdbFile]); % map proper spikeDB to current function workspace
list    = who(DBobj);                     % variable list in matfile
spikeDB = DBobj.(list{1});                % use parenses to access the unknown variable in matfile


%% Compute all cells' peak firing angle in all events
for DNUM = 1:numel(datefolderlist)        % DNUM = date number
    fprintf('-->Start processing %s\n',datefolderlist{DNUM});
    date = datefolderlist{DNUM};   date = strsplit(date,'\');   date = date{3};
    dataind = strcmp(spikeDB.pinfo.general.datedir,date);
    cellid_date = spikeDB.data.grouplist.groupindex;
    cellid_date = cellid_date{1}(dataind);
    
    % Read RED and GREEN diode .pos file for all sessions of interest
    % ===============================================================
    posPath = [datefolderlist{DNUM},'\final\pos\'];
    poslist = dir([posPath,'*.pos']);
    str     = {poslist.name};   % .pos file cellarray
    k = 1;
    while 1    % read in times, X,Y for each session, one session at a time
        [s,v] = listdlg('Name','CellActiveAngle','PromptString','Select ONE pair of RED and Green diole .pos File',...
            'ListString',str,'ListSize',[250,350]);
        if v == 0
            break
        end
        [T{k},Xgreen{k},Ygreen{k}] = readPosFile(posPath,str{s(1)});         % green diode coordinate, position time stamp T
        [~,Xred{k},Yred{k}] = readPosFile(posPath,str{s(2)});                % red diode coordinate
        tok{k} = strtok(str{s(2)},'_');   tok{k} = strrep(tok{k},'VT1','');  % 'tok' = session name
        k = k + 1;                                                           % k = session number
    end
    
    % Cycle thru every session of interest
    for SNUM = 1:numel(T)       % SNUM = session number
        
        % Compute angle distribution in degree in current session, t and angle DO NOT change with event name
        %----------------------------------------------------------------------------------------------------
        int   = str2double(answer{2});  % frame interval
        Xgreen{SNUM} = Xgreen{SNUM}(1:int:end); Ygreen{SNUM} = Ygreen{SNUM}(1:int:end); % X,Y coordinate taken every 'int'-th frame
        Xred{SNUM}   = Xred{SNUM}(1:int:end);   Yred{SNUM}   = Yred{SNUM}(1:int:end);
        t     = T{SNUM}(1:int:end);     % timestamps of sampled frames
        angles = atan2d((Yred{SNUM}-Ygreen{SNUM}),(Xred{SNUM}-Xgreen{SNUM}));  % head direction, [-180,180]
        
        % Find event files in current session
        %--------------------------------------
        evtPath = [datefolderlist{DNUM},'\final\events\'];
        evtlist = dir([evtPath,tok{SNUM},'_c*','.evt']);   % structure array of turn event files
        
        % Cycle thru every event file
        %--------------------------------------
        for ENUM = 1:numel(evtlist)      % ENUM = event file number
            
            % Get spike timestamps of cells active in current event file
            %------------------------------------------------------------
            [~,~,ep_time] = readEvtFile(evtPath,evtlist(ENUM).name);    % evt start/end time in ep_time
            C = strsplit(evtlist(ENUM).name,{'_','.'});                 % multiple delimiter,C{2} = 'cwturn' or 'ccwturn'
            
            % Get sampled frames' timestamps and corresponding head direction within all events
            %-----------------------------------------------------------------------------------
            t_ind = arrayfun(@(x,y) find(t>x & t<y),ep_time(:,1),ep_time(:,2),'UniformOutput',0);  % within-event spikes index
            t_ind_vec = cell2mat(t_ind);
            angle_evt = cellfun(@(x) angles(x),t_ind,'UniformOutput',0); % head direction cell array within events
            angle_evt_vec = angles(t_ind_vec);                           % vectorized head direction arrays
            
            spiketime  = spikeDB.data.spike.spiketime(cellid_date);  % spikes of cells active in evt
            
            % Find spikes' corresponding head directions for each cell
            %-----------------------------------------------------------
            disp('-->Finding spikes'' corresponding head directions')
            [Aspike_evt,eventID_evt] = cellfun(@(x) spk2ang(x),spiketime,'UniformOutput',0);
            noevtspkid   = cellfun(@isempty,Aspike_evt);  % cell ID with no spike in any evt
            Aspike_cell  = Aspike_evt(~noevtspkid);       % angle trains of cells having evt spike
            eventID_cell = eventID_evt(~noevtspkid);      % cellarray of event ID vectors with at east 1 spike in any evt
            cellid_date   = cellid_date(~noevtspkid);       % if EvtMeanRate used to select cell, this line meaningless
            
            % Convert head direction angle to turned angle[0,360], starting from 0 degree
            %--------------------------------------------------------------------------------
            % Image in MATLAB is flipped upside-down relative to video image. Set true north in video
            % (true south in matlab coordinates, which is -90 returned by atan2d) as turning onset:
            % if turn CW in MATLAB coordinate:   -90=>-180/180=>90=>0=>-90, turnning angle 0=>360;
            % if turn CCW in MATLAB coordinate:  -90=>0=>90=>180/-180=>-90, turnning angle 0=>-360.
            
            disp('-->Convert head direction angle to turned angle[0,360], starting from 0 degree')
            % For angle_evt_vec-overal spike angles across all events
            id = (angle_evt_vec >= -180) & (angle_evt_vec <= -90);  % convert head directions corresponding to each frame within all events
            if strcmpi(C{2},'cwturn')
                angle_evt_vec(id)  = angle_evt_vec(id).*(-1) - 90;
                angle_evt_vec(~id) = angle_evt_vec(~id).*(-1) + 270;
            elseif strcmpi(C{2},'ccwturn')
                angle_evt_vec(id)  = angle_evt_vec(id).*(-1) - 450;
                angle_evt_vec(~id) = angle_evt_vec(~id).*(-1) - 90;
            end
            
            % For angle_evt
            if strcmpi(C{2},'cwturn')
                for k = 1:numel(angle_evt)
                    id = (angle_evt{k} >= -180) & (angle_evt{k} <= -90);
                    angle_evt{k}(id)  = angle_evt{k}(id).*(-1) - 90;
                    angle_evt{k}(~id) = angle_evt{k}(~id).*(-1) + 270;
                end
            elseif strcmpi(C{2},'ccwturn')
                for k = 1:numel(angle_evt)
                    id = (angle_evt{k} >= -180) & (angle_evt{k} <= -90);
                    angle_evt{k}(id)  = angle_evt{k}(id).*(-1) - 450;
                    angle_evt{k}(~id) = angle_evt{k}(~id).*(-1) - 90;
                end
            end
            
            % For Aspike_cell-head direction at each spike
            Aturn_cell = cell(1,numel(Aspike_cell));  % initiate turn angle cellarray
            if strcmpi(C{2},'cwturn')
                for k = 1:numel(Aspike_cell)              % converting algorithm
                    id = (Aspike_cell{k} >= -180) & (Aspike_cell{k} <= -90);
                    Aturn_cell{k}(id)  = Aspike_cell{k}(id).*(-1) - 90;
                    Aturn_cell{k}(~id) = Aspike_cell{k}(~id).*(-1) + 270;
                end
            elseif strcmpi(C{2},'ccwturn')
                for k = 1:numel(Aspike_cell)
                    id = (Aspike_cell{k} >= -180) & (Aspike_cell{k} <= -90);
                    Aturn_cell{k}(id)  = Aspike_cell{k}(id).*(-1) - 450;
                    Aturn_cell{k}(~id) = Aspike_cell{k}(~id).*(-1) - 90;
                end
            end
            
            % Find peak firing rate angle for each cell in each events
            %-----------------------------------------------------------
            angbin = str2double(answer{3});
            if strcmpi(C{2},'cwturn')
                edge = 0:angbin:360;     % angle bin of animal head direction
                figedge = edge(2:end);   % used for polar, campass, cartesian plot
            elseif strcmpi(C{2},'ccwturn')
                edge = -360:angbin:0;    figedge = edge(1:end-1);
            end
            
            % Smoothing kernel
            kerlen = str2double(answer{4});   kerstep = str2double(answer{5});   kersigma = str2double(answer{6});
            kernel = pdf('Normal',-kerlen:kerstep:kerlen,0,kersigma);      % smooth kernel, may need tweak parameter
            
            ang_pk_evt = cell(1,size(ep_time,1));
            for ep = 1:size(ep_time,1)   % campass plot for smoothed angle for each event, cell k
                
                Aturn_cell_evt = cell(1,numel(Aturn_cell));
                fr_cell = cell(1,numel(Aturn_cell));
                for k = 1:numel(Aturn_cell)
                    s = eventID_cell{k} == ep;
                    Aturn_cell_evt{k} = Aturn_cell{k}(s);
                    if isempty(Aturn_cell_evt{k})
                        continue
                    else
                        N_ep       = histcounts(angle_evt{ep},edge,'Normalization','probability');  % stay prob in each space bin for each evt
                        %N_ep       = conv(N_ep,kernel,'same');
                        t_evt{ep}   = ep_time(ep,2) - ep_time(ep,1); % time spent in turn event ep
                        t_evt{ep}   = N_ep.*t_evt{ep};
                        binspk_ep  = histcounts(Aturn_cell_evt{k},edge);
                        %binspk_ep  = conv(binspk_ep,kernel,'same');
                        fr_cell{k} = binspk_ep./t_evt{ep};           % mean firing rate in each head direction bin
                        
                        % Handle inf/nan in fr_cell{k}:
                        % if no spike in current bin--binspk_ep=0,fr_cell=0;
                        % if no position data in current bin--t_evt_ep=0,fr_cell=inf;
                        % if both are 0--fr_cell=NaN.
                        % After replace inf and NaN to 0, if cell has 0 in all bins (sum(fr_cell{k})==0), need to treat
                        % tr_cell{k} = [], otherwise fr_cell{k}./max(fr_cell{k}) generate all NaN, the returned index
                        % for [pk_cell,ang_pk_evt{ep}] = cellfun(@max,fr_cell,'UniformOutput',0) will be 1,that's WRONG.
                        % After treat tr_cell{k} = [], the returned index for fr_cell{k}./max(fr_cell{k}) is [].
                        % This [] index will replaced by NaN so this cell has NO preferred peak.
                        fr_cell{k}(isinf(fr_cell{k})|isnan(fr_cell{k})) = 0;
                        if sum(fr_cell{k}) > 0
                            fr_cell{k} = conv(fr_cell{k},kernel,'same');
                        else
                            fr_cell{k} = [];
                        end
                    end
                    fr_cell{k} = fr_cell{k}./max(fr_cell{k});  % normalize max prob to 1
                end
                
                % Handle cells with empty spike train and order cells within each event
                [pk_cell,ang_pk_evt{ep}] = cellfun(@max,fr_cell,'UniformOutput',0);
                emptyid = cellfun(@isempty,pk_cell);
                [ang_pk_evt{ep}{emptyid}] = deal(nan);
                ang_pk_evt{ep} = cell2mat(ang_pk_evt{ep});
                ang_pk_evt{ep}(~emptyid) = figedge(ang_pk_evt{ep}(~emptyid));
                
            end
            
            % Export variables to workspace
            %----------------------------------------------------------
            fprintf('-->Saved %s %s %s ordered cells figures to workspace\n',datefolderlist{DNUM},C{1},C{2});
            PeakAngle = struct('day',date,'direction',C{2},'cellID_act',cellid_date,'Angles',figedge,...
                'TimeInBinsEvt',{t_evt},'PeakAngleEvt',{ang_pk_evt});
            D   = strsplit(datefolderlist{DNUM},'\');    % get animalname+date in D{3}
            D   = strsplit(D{3},'-');
            var = [D{1},'_',D{2},'_',D{3},'_',D{4},C{1},C{2},'_PeakAngle'];
            assignin('base',var,PeakAngle)
            
        end % end of all event files
    end % end of all sessions
    fprintf('-->Finish processing %s\n',datefolderlist{DNUM});
end % end of all dates
disp('********************')

end % EOF


%% Subfunctions
function [Aspike_evt,eventID_evt] = spk2ang(spiketime)
global ep_time t angles

spike_evt   = [];   % initiate spike vector for all cell
eventID_evt = [];   % initiate each spike's event ID for cell k

% Find spikes' timestamps within each event,k
for k = 1:size(ep_time,1)
    s = spiketime(spiketime>ep_time(k,1) & spiketime<ep_time(k,2));  % within-event spikes
    eventID_evt = [eventID_evt; k*ones(length(s),1)];                % get Y-axis=number of events
    spike_evt   = [spike_evt; s];
end

% Exclude spikes without pre or post location bound, ensure size(eventID_evt)==size(spike_evt)
pre  = arrayfun(@(x) find(t<x,1,'last'),spike_evt,'UniformOutput',0);  % array of bounds of spikes in spike_evt
post = arrayfun(@(x) find(t>x,1,'first'),spike_evt,'UniformOutput',0);

bd = cellfun(@(x,y) ~isempty(x) && ~isempty(y),pre,post,'UniformOutput',0);
bd = cell2mat(bd);    % logical index vector of spikes with pre and post bound
eventID_evt(~bd) = [];    spike_evt(~bd)   = [];

% Interpolate head direction angle corresponding to each spike
T_bound     = [t([pre{bd}]), t([post{bd}])];         % time bound for spikes in s
angle_bound = [angles([pre{bd}]), angles([pre{bd}])];  % angle bound for spikes in s
Aspike_evt  = [];                                    % initiate head direction corresponding to each spike for cell k
for n = 1:numel(spike_evt)
    A          = interp1(T_bound(n,:), angle_bound(n,:), spike_evt(n));
    Aspike_evt = [Aspike_evt; A];                    % head direction angle
end

end % end of spk2ang


