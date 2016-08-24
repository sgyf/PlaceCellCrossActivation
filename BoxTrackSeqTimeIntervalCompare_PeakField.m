function BoxTrackSeqTimeIntervalCompare_PeakField
% Part of BoxTrackSeqTimeIntervalCompare. This function compute all cells' session averaged peak
% firing location on track, time interval between neighbouring peaks.
%
% Created on 8/24/2016, XM.haha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Choose date and session, import SpikeDB

clear; clc;
global ep_time ep_unwanted t pos

prompt   = {'Input full folders path (eg, E:\Dan11\Dan11-2013-10-17), delimitered by white space',...
    'Position bin size(pix)',...
    'Position smooth parameter-Half kernel length',...
    'Position smooth parameter-kernel step',...
    'Position smooth parameter-Sigma'};
dlgtitle = 'Parameter to compute cell active place';
nline    = [1,85;1,85;1,85;1,85;1,85];
default  = {'','6','5','1','2'};
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


%% Compute all cells' peak firing location on track
for DNUM = 1:numel(datefolderlist)        % DNUM = date number
    fprintf('-->Start processing %s\n',datefolderlist{DNUM});
    date = datefolderlist{DNUM};   date = strsplit(date,'\');   date = date{3};
    dataind = strcmp(spikeDB.pinfo.general.datedir,date);
    cellid_date = spikeDB.data.grouplist.groupindex;
    cellid_date = cellid_date{1}(dataind);
    
    % Read ONE reliable .pos file for all sessions of interest
    %----------------------------------------------------------
    posPath = [datefolderlist{DNUM},'\final\pos\'];
    poslist = dir([posPath,'*.pos']);
    str1    = {poslist.name};   % .pos file cellarray
    k = 1;
    while 1    % read in times, X,Y for each session, one session at a time
        [s,v] = listdlg('Name','CellActivePlace','PromptString','Select ONE reliable .pos File',...
            'ListString',str1,'ListSize',[250,350]);
        if v == 0
            break
        end
        [T{k},X{k},Y{k}] = readPosFile(posPath,str1{s});    % diode coordinate, position time stamp T
        posfile{k} = str1{s};   % .pos file name
        tok{k} = strtok(str1{s},'_');   tok{k} = strrep(tok{k},'VT1','');  % 'tok' = session name
        k = k + 1;                                                         % k = session number
    end
    
    % Cycle thru every session of interest
    %-----------------------------------------------------------
    for SNUM = 1:numel(T)       % SNUM = session number
        fprintf('-->Start processing .pos file: %s\n',posfile{SNUM});
        
        t = T{SNUM};   pos = X{SNUM};   % for subfunction spk2pos use
        
        % Read start/end timestamps for all unwanted events
        evtPath = [datefolderlist{DNUM},'\final\events\'];
        evtlist = dir([evtPath,'*.evt']);
        str2    = {evtlist.name};     % .evt file cellarray
        [s,v]   = listdlg('PromptString','Select Unwanted .evt File(ripple,stop...) to exclude spikes within',...
            'ListString',str2,'ListSize',[250,350]);
        if v ~= 0
            for k = 1:numel(s)
                [~,~,ep_unwanted{k}] = readEvtFile(evtPath,str2{s(k)});  % event start/end timestamp
            end
        end
        
        % Cycle thru every event file
        %-------------------------------------------------------------
        while 1
            
            % Get spike timestamps of cells active in current event file
            [evtFile,evtPath] = uigetfile([evtPath,'*.evt'],'Select Event File');
            if ~ischar(evtFile)
                break
            end
            [~,~,ep_time] = readEvtFile(evtPath,evtFile);     % evt start/end time in ep_time
            C = strsplit(evtFile,{'_','.'});                  % multiple delimiter,C{2} = 'leftright' or 'rightleft'
            
            % Get sampled frames' timestamps and corresponding head direction within all events
            t_ind = arrayfun(@(x,y) find(t>x & t<y),ep_time(:,1),ep_time(:,2),'UniformOutput',0);  % within-event spikes index
            t_ind_vec = cell2mat(t_ind);
            %             xi_evt = cellfun(@(x) pos(x),t_ind,'UniformOutput',0); % head direction cell array within events
            xi_evt_vec = pos(t_ind_vec);                           % vectorized head direction arrays
            
            spiketime  = spikeDB.data.spike.spiketime(cellid_date);  % spikes of cells active in evt
            
            % Compute smoothed firing rate
            %---------------------------------------------------------
            % Find spikes' corresponding x-position for each cell
            disp('-->Find spikes'' corresponding x-positions');
            [Xspike_evt,~] = cellfun(@(x) spk2pos(x),spiketime,'UniformOutput',0);
            noevtspkid   = cellfun(@isempty,Xspike_evt);   % cell ID with no event spike
            Xspike_cell  = Xspike_evt(~noevtspkid);
            %             eventID_cell = eventID_evt(~noevtspkid);
            cellid_date   = cellid_date(~noevtspkid);
            
            % Smooth spike trains
            posbin  = str2double(answer{2});
            edge    = min(xi_evt_vec)-10:posbin:max(xi_evt_vec)+10;      % give 10 pix flexibility @ two ends
            figedge = edge(2:end);    % x-axis range for later plot to match FR values
            kerlen  = str2double(answer{3});   kerstep = str2double(answer{4});   kersigma = str2double(answer{5});
            kernel  = pdf('Normal',-kerlen:kerstep:kerlen,0,kersigma);   % smooth kernel, may need tweak parameter
            
            % Use binned spike count to compute smooth firing curve
            %         clname = spikeDB.pinfo.general.clname(cellid_act);  % TT name of cells of interest, no empty event-spikes
            fr = cell(1,numel(Xspike_cell));                    % initiate FR cellarray
            for k = 1:numel(Xspike_cell)
                N      = histcounts(xi_evt_vec,edge,'Normalization','probability');  % overall stay prob in each space bin
                % N      = conv(N,kernel,'same');
                t_evt  = sum(ep_time(:,2) - ep_time(:,1));  % total time spent within run events
                t_evt  = N.*t_evt;                          % time spent in each position bin
                binspk = histcounts(Xspike_cell{k},edge);   % spike number in each bin
                % binspk  = conv(binspk,kernel,'same');         % smooth spikes distribution in bins
                fr{k}  = binspk./t_evt;                     % mean firing rate in each bin
                fr{k}(isinf(fr{k})|isnan(fr{k})) = 0;       % handle nan due to zeros in t_evt, or 'conv' generate all NaNs
                fr{k}  = conv(fr{k},kernel,'same');
            end
            
            % Order cells based on place field position
            %----------------------------------------------------------
            disp('-->Order cells based on place field positins');
            fr = cellfun(@(x) x/max(x),fr,'UniformOutput',0);  % normalize max firing prob to 1 for each cell
            [~,xi_peak_id] = cellfun(@max,fr);
            xi_peak = figedge(xi_peak_id);
            
            % Export variables to workspace
            %----------------------------------------------------------
            fprintf('-->Export to workspace variables: %s %s %s\n',datefolderlist{DNUM},C{1},C{2})
            PeakField = struct('day',{date},'direction',C{2},'cellID_act',cellid_date,'cellorder',xi_peak_id,...
                'TimeInBins',t_evt,'PositionPix',figedge,'PeakField_avg',xi_peak);
            D   = strsplit(datefolderlist{DNUM},'\');    % get animalname+date in D{3}
            D   = strsplit(D{3},'-');
            var = [D{1},'_',D{2},'_',D{3},'_',D{4},C{1},C{2},'_PeakField'];
            assignin('base',var,PeakField)
        end
    end
    fprintf('-->Finish processing %s\n',datefolderlist{DNUM});
end % end of all dates
disp('********************')

end % EOF


%% Subfunctions
function [Xspike_evt,eventID_evt] = spk2pos(spiketime)
global ep_time ep_unwanted t pos

% Find spikes' timestamps within each event,k
Tspike_evt  = [];  % initiate each spike's timestamp for cell k
eventID_evt = [];  % initiate each spike's event ID for cell k
for k = 1:size(ep_time,1)
    s = spiketime(spiketime>ep_time(k,1) & spiketime<ep_time(k,2));  % within-event spikes
    eventID_evt = [eventID_evt; k*ones(length(s),1)];                % get Y-axis=number of events
    Tspike_evt  = [Tspike_evt; s];
end

% Find spikes' timestamps within unwanted events
T_unwanted = [];  % initiate Ripple evt spike timestamp for cell k
for k = 1:numel(ep_unwanted)
    for m = 1:size(ep_unwanted{k},1)
        s = spiketime(spiketime>ep_unwanted{k}(m,1) & spiketime<ep_unwanted{k}(m,2));  % unwanted spikes
        T_unwanted = [T_unwanted; s];
    end
end

% Exclude both spike timestamps in unwanted evt, and their corresponding eventID
id = [];
for tripp = T_unwanted'   % vector after '=' must be row vector
    id = [id;find(Tspike_evt == tripp)];
end
eventID_evt(id) = [];
Tspike_evt      = setdiff(Tspike_evt,T_unwanted);

% Exclude spikes without pre or post location bound,ensure size(eventID_evt)==size(Tspike_evt)
pre  = arrayfun(@(x) find(t<x,1,'last'),Tspike_evt,'UniformOutput',0);    % cellarray of bounds of spikes in s
post = arrayfun(@(x) find(t>x,1,'first'),Tspike_evt,'UniformOutput',0);

bd = cellfun(@(x,y) ~isempty(x) && ~isempty(y),pre,post,'UniformOutput',0);
bd = cell2mat(bd);    % logical index of spikes with pre and post bound
eventID_evt(~bd) = [];    Tspike_evt(~bd)  = [];

% Interpolate x-position corresponding to each spike
T_bound    = [t([pre{bd}]), t([post{bd}])];     % time bound for spikes in s
X_bound    = [pos([pre{bd}]), pos([pre{bd}])];  % angle bound for spikes in s
Xspike_evt = [];  % initiate head direction corresponding to each spike for cell k

for n = 1:numel(Tspike_evt)
    x_inbound  = interp1(T_bound(n,:), X_bound(n,:), Tspike_evt(n));
    Xspike_evt = [Xspike_evt; x_inbound];
end

end % end of spk2pos



%% Subfunction spk2ang
function [Aspike_evt,eventID_evt] = spk2ang(spiketime)
global ep_time t angle

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
angle_bound = [angle([pre{bd}]), angle([pre{bd}])];  % angle bound for spikes in s
Aspike_evt  = [];                                    % initiate head direction corresponding to each spike for cell k
for n = 1:numel(spike_evt)
    A          = interp1(T_bound(n,:), angle_bound(n,:), spike_evt(n));
    Aspike_evt = [Aspike_evt; A];                    % head direction angle
end

end % end of spk2ang


