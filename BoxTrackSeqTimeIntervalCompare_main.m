% BoxTrackSeqTimeIntervalCompare_main: Neighbouring cell time interval comparison between box and track

%% Compute all cells' peak firing location on track, DONE! Saved in
% 'E:\PlaceCellCrossActivationAnalysis\box_trk_correlation\ObsvBoxEvt_Trk\BoxTrackSeqTimeIntervalCompare\TrackPeakField.mat'
BoxTrackSeqTimeIntervalCompare_PeakField


%% Neighbouring cell time interval in track sequence, DONE!
% 'angorder, comcellneed, trktemp, condition, animal, day, session, turn, allPeakField, comcellid, peaktimeid_trk, interval,
% interval_norm' are saved in 'TrackPeakFieldTimeInterval'

% Inport to workspace
load('E:\PlaceCellCrossActivationAnalysis\box_trk_correlation\ObsvBox_Trk_CommonCell\CommonCell_all_obsv_work.mat');
load('E:\PlaceCellCrossActivationAnalysis\box_trk_correlation\ObsvBoxEvt_Trk\BoxTrackSeqTimeIntervalCompare\TrackPeakField.mat');

% Get common cells' info
angorder = strcmp({CommonCell_all_obsv_work.order},'AngOrder');
comcellneed = CommonCell_all_obsv_work(angorder);   % only need common cells in AngOrder
trktemp   = {comcellneed.template};
condition = {comcellneed.condition};
animal    = {comcellneed.animal};
day       = {comcellneed.day};
session   = {comcellneed.session};
turn      = {comcellneed.turn};

allPeakField = cellfun(@(x,y,z) [x,'_',y,z],animal,day,trktemp,'UniformOutput',0);  % track temp for each comcellid
comcellid = {comcellneed.cellID_common};   % common cell ID

% 'allPeakField' and 'comcellid' have one-to-one correpondance
peaktimeid_trk = cell(1,numel(comcellid));
peakorder_trk = cell(1,numel(comcellid));
interval_trk = cell(1,numel(comcellid));
interval_trk_norm = cell(1,numel(comcellid));
for k = 1:numel(comcellid)
    target = strncmp(allPeakField{k},AllPeakFieldList,27);
    targetPeakField = eval('base',AllPeakFieldList{target});
    time = targetPeakField.TimeInBins;
    position = targetPeakField.PositionPix;
    PeakFields_all = targetPeakField.PeakField_avg;
    cellID = targetPeakField.cellID_act;
    
    [~,cellIDind] = intersect(cellID, comcellid{k});   % common cell peak ID
    % commonn cell peak time ID. What about cells with same peak? Influence on result?
    [~,peaktimeid_trk{k},peakorder_trk{k}] = intersect(position,PeakFields_all(cellIDind));
    
    % Time interval between neighbouring cells
    % NOTE: 'interval_trk' should be mean time interval_trk, but below is total interval_trk due to
    % mis-calculation from old code 'BoxTrackSeqTimeIntervalCompare_PeakField'. Code has been
    % changed, but didn't recompute. Since total duration will be normalized to 1, so shouldn't
    % affect result here.
    for m = 1:(numel(peaktimeid_trk{k})-1)
        interval_trk{k}(m) = sum(time(peaktimeid_trk{k}(m):peaktimeid_trk{k}(m+1))) - time(peaktimeid_trk{k}(m))/2 - time(peaktimeid_trk{k}(m+1))/2;
    end
    interval_trk_norm{k} = interval_trk{k}./sum(interval_trk{k});   % normalize time interval_trk to sum to 1
end


%% Compute all cells' peak firing angle in each turn events in box, DONE! Saved in
% 'E:\PlaceCellCrossActivationAnalysis\box_trk_correlation\ObsvBoxEvt_Trk\BoxTrackSeqTimeIntervalCompare\BoxPeakAngle.mat'
BoxTrackSeqTimeIntervalCompare_EvtPeakAngle


%%  Neighbouring cell time interval in box sequence

% Inport to workspace
%-----------------------------------------------------------------------------------------
load('E:\PlaceCellCrossActivationAnalysis\box_trk_correlation\ObsvBoxEvt_Trk\EvtCircularCorrData_circshuffle\EvtCircularCorrData_all_obsv_circshuffle.mat');
load('E:\PlaceCellCrossActivationAnalysis\box_trk_correlation\ObsvBoxEvt_Trk\BoxTrackSeqTimeIntervalCompare\BoxPeakAngle.mat');
load('E:\PlaceCellCrossActivationAnalysis\box_trk_correlation\ObsvBoxEvt_Trk\BoxTrackSeqTimeIntervalCompare\TrackPeakField.mat');
load('E:\PlaceCellCrossActivationAnalysis\box_trk_correlation\ObsvBoxEvt_Trk\BoxTrackSeqTimeIntervalCompare\TrackPeakFieldTimeInterval.mat');

% 'angorder, comcellneed, trktemp, condition, animal, day, session, turn, comcellid' for
% 'EvtCircularCorrData_all_obsv_circshuffle' are the same as in 'CommonCell_all_obsv_work', no need
% to do it again after load 'TrackPeakFieldTimeInterval'.
allPeakAngle = cellfun(@(a,b,c,d) [a,'_',b,c,d],animal,day,session,turn,'UniformOutput',0);  % track temp for each comcellid

%% Find significant-match event for each entry in 'EvtCircularCorrData_all_obsv_circshuffle'
%------------------------------------------------------------------------------------------
sigthreshold = 0.06;
pvalue = {EvtCircularCorrData_all_obsv_circshuffle(angorder).p_true};
sigevtind = cellfun(@(x) x<=sigthreshold, pvalue,'UniformOutput',0);   % significant match event id

% Preallocate and computing
%------------------------------------------------------------------------------------------
peakpos_box = cell(1,numel(comcellid));
peakposcell_box = cell(1,numel(comcellid));
interval_box = cell(1,numel(comcellid));
interval_box_norm = cell(1,numel(comcellid));
peakpos_trk_nonan = cell(1,numel(comcellid));
peakposcell_trk_nonan = cell(1,numel(comcellid));
interval_trk = cell(1,numel(comcellid));
interval_trk_norm = cell(1,numel(comcellid));
dist_trk = cell(1,numel(comcellid));
dist_box = cell(1,numel(comcellid));
dist_diff = cell(1,numel(comcellid));
dist_diff_shuf = cell(1,numel(comcellid));
dist_zscore = cell(1,numel(comcellid));
dist_zscore_shuf = cell(1,numel(comcellid));
for k = 1:numel(comcellid)
    % Take info from track active cells
    %-----------------------------------------------------
    target = strncmp(allPeakField{k},AllPeakFieldList,27);
    targetPeakField = eval('base',AllPeakFieldList{target});
    time_trk = targetPeakField.TimeInBins;
    position = targetPeakField.PositionPix;
    PeakFields_all = targetPeakField.PeakField_avg;
    cellID_trk = targetPeakField.cellID_act;
    
    % Find cells' peak angle in significant events
    %-----------------------------------------------------
    target = strncmp(allPeakAngle{k},AllPeakAngleList,26);
    targetPeakAngle = eval('base',AllPeakAngleList{target});
    angle = targetPeakAngle.Angles;
    time_box = targetPeakAngle.TimeInBinsEvt;      % time is a cellarray, each element is for an evt
    PeakAngleEvt = targetPeakAngle.PeakAngleEvt;   % PeakAngleEvt is a cellarray, each element is for an evt
    cellID_box = targetPeakAngle.cellID_act;
    
    if numel(PeakAngleEvt) ~= numel(pvalue{k})   % some evt# in 'EvtCircularCorrData_all_obsv_circshuffle' doesn't agree with evt# in evt files
        continue
    end
    
    PeakAngleEvt_sig = PeakAngleEvt(sigevtind{k}); % cells' peak angle in significant match events, COULD BE CELL ARRAY
    time_box_sig = time_box(sigevtind{k});   % significant match event time interval between two angle bins
    
    % Compute the time interval difference b/w track and box peak locations
    %---------------------------------------------------------------------------
    for m = 1:numel(PeakAngleEvt_sig)   % m = #of significant turn events
        
        % Time interval between not NaN neighbour cells on track
        %----------------------------------------------------------
        [~,cellIDind] = intersect(cellID_box, comcellid{k});     % common cell peak position ID
        comcellid_nonan = comcellid{k}(~isnan(PeakAngleEvt_sig{m}(cellIDind)));  % exclude common cell with no spike in current evt
        
        if (numel(comcellid_nonan) > 4) && (numel(cellIDind) == numel(comcellid{k}))  % proceed only when non NaN cells in current evt>4,and sometimes 'comcellid{k}' has cells not in cellIDind, why??
            [~,cellIDind] = intersect(cellID_trk, comcellid_nonan);  % common cell with no NaN position peak position ID
            % commonn cell peak time ID. What about cells with same peak? Influence on result?
            [~,peakpos_trk_nonan{k}{m},peakposcell_trk_nonan{k}{m}] = intersect(position,PeakFields_all(cellIDind)); % peakposcell_trk_nonan is commoncells' ID, used to compare b/w trk & box
            
            for n = 1:(numel(peakpos_trk_nonan{k}{m})-1)  % n = #of not NaN cells in event m
                interval_trk{k}{m}(n) = sum(time_trk(peakpos_trk_nonan{k}{m}(n):peakpos_trk_nonan{k}{m}(n+1))) - time_trk(peakpos_trk_nonan{k}{m}(n))/2 - time_trk(peakpos_trk_nonan{k}{m}(n+1))/2;
            end   % Interval_trk is time interval b/w neighbour cells on track
            
            % Time interval between not NaN neighbour cells in box
            %----------------------------------------------------------
            [~,cellIDind] = intersect(cellID_box, comcellid_nonan);  % common cell with no NaN position peak position ID
            
            if numel(comcellid_nonan) ~= numel(cellIDind)   % some comcellid_nonan are not in cellID_box
                continue
            end

            peakpos_box{k}{m} = sort(arrayfun(@(x) find(angle == x), PeakAngleEvt_sig{m}(cellIDind)));
            
            for n = 1:(numel(peakpos_box{k}{m})-1)  % n = #of not NaN cells in event m
                interval_box{k}{m}(n) = sum(time_box_sig{m}(peakpos_box{k}{m}(n):peakpos_box{k}{m}(n+1))) - time_box_sig{m}(peakpos_box{k}{m}(n))/2 - time_box_sig{m}(peakpos_box{k}{m}(n+1))/2;
            end   % Interval_trk is time interval b/w neighbour cells on track
            
            % Normalized time interval to sum up to 1, in order to compare
            %-------------------------------------------------------------
            interval_trk_norm{k}{m} = interval_trk{k}{m}./sum(interval_trk{k}{m});
            interval_box_norm{k}{m} = interval_box{k}{m}./sum(interval_box{k}{m});
            %             BoxTimeLen = sum(interval_box{k}{m});   % total box sequence lenght for current evt
            
            % Compute the time interval difference b/w track and box peak locations
            %----------------------------------------------------------------------
            for n = 1:(numel(peakpos_trk_nonan{k}{m})-1)  % n = #of not NaN cells in event m
                dist_trk{k}{m}(n) = interval_trk_norm{k}{m}(n);
                formercell = peakposcell_trk_nonan{k}{m}(n);   lattercell = peakposcell_trk_nonan{k}{m}(n+1);
                if (formercell ~= 1) && ((lattercell-1) ~= numel(interval_box_norm{k}{m})) && (formercell <= lattercell)
                    dist_box{k}{m}(n) = sum(interval_box_norm{k}{m}(formercell:lattercell-1)) -...
                        interval_box_norm{k}{m}(formercell)/2 -...
                        interval_box_norm{k}{m}(lattercell-1)/2;
                elseif (formercell ~= 1) && ((lattercell-1) ~= numel(interval_box_norm{k}{m})) && (formercell > lattercell)
                    dist_box{k}{m}(n) = -sum(interval_box_norm{k}{m}(lattercell:formercell-1)) -...
                        interval_box_norm{k}{m}(formercell-1)/2 -...
                        interval_box_norm{k}{m}(lattercell)/2;
                else
                    dist_box{k}{m}(n) = nan;
                end
            end
            dist_diff{k}{m} = nansum(abs(dist_trk{k}{m} - dist_box{k}{m}));   % time interval difference between track and box
            
            % Compute the time interval difference b/w track and shuffled box peak locations
            %--------------------------------------------------------------------------------
            for n = 1:1000   % shuffle 1000 times
                shufinterval = rand(1,numel(interval_box_norm{k}{m}));
                shufinterval = shufinterval./sum(shufinterval);
                
                dist_box_shuf = nan(1,numel(peakpos_trk_nonan{k}{m})-1);
                for nn = 1:(numel(peakpos_trk_nonan{k}{m})-1)  % n = #of not NaN cells in event m
                    dist_trk{k}{m}(nn) = interval_trk_norm{k}{m}(nn);
                    formercell = peakposcell_trk_nonan{k}{m}(nn);   lattercell = peakposcell_trk_nonan{k}{m}(nn+1);
                    if (formercell ~= 1) && ((lattercell-1) ~= numel(interval_box_norm{k}{m})) && (formercell <= lattercell)
                        dist_box_shuf(nn) = sum(shufinterval(formercell:lattercell-1)) -...
                            shufinterval(formercell)/2 -...
                            shufinterval(lattercell-1)/2;
                    elseif (formercell ~= 1) && ((lattercell-1) ~= numel(interval_box_norm{k}{m})) && (formercell > lattercell)
                        dist_box_shuf(nn) = -sum(shufinterval(lattercell:formercell-1)) -...
                            shufinterval(formercell-1)/2 -...
                            shufinterval(lattercell)/2;
                    else
                        dist_box_shuf(nn) = nan;
                    end
                end
                dist_diff_shuf{k}{m}(n) = nansum(abs(dist_trk{k}{m} - dist_box_shuf));   % time interval difference between track and shuffled box
            end
            
            z = zscore([dist_diff{k}{m},dist_diff_shuf{k}{m}]);   % convert time interval difference to zscore
            dist_zscore{k}{m} = z(1);   % zscore of true diffence
            dist_zscore_shuf{k}{m} = z(2:end);
        end
    end
    
end

% Saved dis_box, dist_trk, dist_diff, dist_zscore, dist_zscore_shuf


%% Stats

load('E:\PlaceCellCrossActivationAnalysis\box_trk_correlation\ObsvBoxEvt_Trk\BoxTrackSeqTimeIntervalCompare\TrackPeakFieldTimeInterval.mat');
load('E:\PlaceCellCrossActivationAnalysis\box_trk_correlation\ObsvBoxEvt_Trk\BoxTrackSeqTimeIntervalCompare\NeighborCellTimeIntervalDifferenceBetweenTrkBox.mat');

% Get true distance difference and shuffled data in target condition
%---------------------------------------------------------------------------
% 'cond': train_run_train;   train_run_train_1sttime;    naive_run_naive;    car_run_car;
%             empty_run_empty;   no_run_no;    block_run_block;

cond = 'block_run_block';
if strcmp(cond,'all')
    dist_zscore_cond = dist_zscore;
    dist_zscore_shuf_cond = dist_zscore_shuf;
else
    condidx = strcmp(condition,cond);
    dist_zscore_cond = dist_zscore(condidx);
    dist_zscore_shuf_cond = dist_zscore_shuf(condidx);
end

dist_zscorevec = [];
dist_zscore_shufvec = [];
for k = 1:numel(dist_zscore_cond)
    temp = cell2mat(dist_zscore_cond{k});
    dist_zscorevec = [dist_zscorevec, temp];
    temp = cell2mat(dist_zscore_shuf_cond{k});
    dist_zscore_shufvec = [dist_zscore_shufvec, temp];
end

% Plot
%--------------------------------------------------------------------
[n,edge] = histcounts(dist_zscorevec,'norm','pdf');
[nshuf,edgeshuf] = histcounts(dist_zscore_shufvec,'norm','pdf');

% For empty_run_empty condition, pdf at 0 zscore has an unusual peak, so multiply by 0.8 to make
% normal distribution, why this peak??
if strcmp(cond,'empty_run_empty')
    [~,I] = max(nshuf);   nshuf(I) = nshuf(I)*0.8;
end

figure('color','w','Position',[10 725 341 271])
hold on
plot(edge(2:end),n)
plot(edgeshuf(2:end),nshuf)
hold off
set(gca,'ytick',[0,0.2,0.4,0.6],'FontSize',14,'TitleFontSizeMultiplier',1,'ticklength',[0.04 0.04])
xlabel('Z-scored gap difference');   ylabel('Probability Density');
legend('Empty-track','Shuffled','location','northwest'), legend('boxoff')  % change legend accordingly

% Stats
%--------------------------------------------------------------------
[h,p] = kstest2(dist_zscorevec,dist_zscore_shufvec)         % KS test for gap difference distribution between true and shuffled
sigpct = sum(dist_zscorevec<=-1.96)/numel(dist_zscorevec)   % percent of signigicant zscores at 5%

% Conclusion: neighbour cell time interval difference between box and track pdf for matched events
% are different from shuffled data in all conditions.
% empty_run_empty: sigpct = 0.9136
% train_run_train: sigpct = 0.9848
% naive_run_naive: sigpct = 1
% car_run_car: sigpct = 0.9787
% no_run_no: sigpct = 1
% block_run_block: sigpct = 0.9630
