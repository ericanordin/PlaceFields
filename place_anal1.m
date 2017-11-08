sessions = [[8068 53]
    [8068 55]];
numSessions = size(sessions, 1);

%Pre-allocate avg_map
avg_map = cell(numSessions, 1);
all_maps_sum = cell(numSessions,1);
maxRotations = 359; %Highest degree of rotation to be checked, in increments of 1 degree

for sess_i = 1:numSessions
    
    if 1
        rat = sessions(sess_i, 1);
        cur_sess = sessions(sess_i, 2);
        sess_info = session_info(rat, cur_sess, 0); % the last 1 extracts depths
        saved_centers_file = ['R' num2str(rat) 'Sess' num2str(cur_sess) 'PlaceFieldCenters.mat'];
        
        
        
        % Load the events file to obtain the epochs info
        disp(['loading Events file: ' sess_info.events_file]);
        [event_data, epoch_times] = read_events4(sess_info.events_file);
        
        % load video tracker data
        disp(['loading PVD files: ' sess_info.pos_file1 ' & ' sess_info.pos_file2]);
        %pv_data: Video data, position of rat
        %C1: Time at various positions
        %C2: x-position
        %C3: y-position
        pv_data = load(sess_info.pos_file1);
        pv_data = pv_data(logical([diff(pv_data(:,1))>0; 1]),:);
        pv_data = pv_data(:,1:3); % discard some data to save memory
        disp('done');
        
        tfile_list = FindFiles('*.t', 'StartingDirectory', sess_info.tfile_path);
        ts_data = LoadSpikes(tfile_list);
        keep_celli = true(1,length(ts_data));
        %keep_celli: 1 marks tetrodes to keep; 0 marks tetrodes to ignore
        
        %Pre-allocation of arrays
        cell_info(1:length(ts_data)) = struct('tet_id', 0, 'cell_num', 0, 'filename', '', 'counts', [0 0 0]);
        spike_data = cell(1,length(ts_data));
        spike_pos_data = cell(1,length(ts_data));
        
        for i = 1:length(ts_data)
            
            % extract tetrode and cell number from filename
            [path, name, ext] = fileparts(tfile_list{i});
            tet_name_i = strfind(name, 'TT')+2;
            [tt_id, cell_id] = strread(name, 'TT%dnsc_%d');
            if tt_id>6  % block tetrodes 7-12, which are in frontal cortex
                keep_celli(i) = 0;
            end
            cell_info(i).tet_id = tt_id;
            cell_info(i).cell_num = cell_id;
            cell_info(i).filename = tfile_list{i};
            
            times = Data(ts_data{i}); % extract data from the silly ts object
            times = times*100;  % convert to microseconds
            
            spike_data{i} = times;
            
            % count spikes in each epoch
            %epoch 1 = sleep, 2 = run, 3 = sleep.
            %epoch_times: row = epoch, column = start/stop time
            epoch_counts = NaN*zeros(1,3);
            for j = 1:size(epoch_times,1)
                epoch_counts(j) = sum(times>epoch_times(j,1) & times<epoch_times(j,2));
                %times>epoch_times(j,1) produces an array declaring whether
                %times is greater for EACH value of the times array
                %Removes any times not in the epoch. Only concerns run time epoch,
                %not sleep.
            end
            cell_info(i).counts = epoch_counts;
            if cell_info(i).counts(2)<50 %Cells with less than 50 spikes are cut (non-place cells)
                keep_celli(i) = 0; %Whether cell should be processed
                disp(['Cell TT' num2str(cell_info(i).tet_id) '_' num2str(cell_info(i).cell_num) ' rejected: ' num2str(cell_info(i).counts(2)) ' spikes']);
            end
            
            disp(['Loaded: ' tfile_list{i}]);
            
            %Spike data is collected more often than position data
            %Interpolates the position of the rat for spikes that do not
            %correspond to recorded positions.
            %pv_data(:,1) = timestamp {x}
            %[pv_data(:,2) pv_data(:,3)] = x-position y-position {f(x)}
            %times {x'}: new values to be interpolated into f(x) data
            %Output = [x y] interpolated data array
            posit = interp1(pv_data(:,1), [pv_data(:,2) pv_data(:,3)], times);
            spike_pos_data{i} = [times posit]; %3 columns: times, x, y
            
        end
        spike_pos_data = spike_pos_data(keep_celli);
        cell_info = cell_info(keep_celli);
        spike_data = spike_data(keep_celli);
        ncells = length(spike_data);
        
    end
    
    
    nbins = 50;
    smoothfac = .75; %Smoothing factor, determines the degree to which the neighboring points impact the central point (kernel values)
    locx.min = min(pv_data(:,2));
    locx.max = max(pv_data(:,2));
    locy.min = min(pv_data(:,3));
    locy.max = max(pv_data(:,3));
    %NOT Matlab's built in histogram2 function
    %Occupancy map: how much time spent in a given location
    occ_map = histogram2(pv_data(:,2)', pv_data(:,3)', [locx.min locx.max nbins; locy.min locy.max nbins]);
    occ_map(occ_map==0) = NaN;
    
    if exist(saved_centers_file, 'file')
        cur_ncells = ncells;
        load(saved_centers_file) %loading overwrites ncells with value from file
        if cur_ncells ~= ncells
            disp(['Cells in stored file ' saved_centers_file ' dont match current data.  Recompute centers.']);
            return;
        end
    else
        centerx = ones(1, length(spike_data));
        centery = ones(1, length(spike_data));
        for i = 1:ncells
            %Firing map: how many spikes fired in a given location
            %Maps out the entire recording area
            firing_map = histogram2(spike_pos_data{i}(:,2)', spike_pos_data{i}(:,3)', [locx.min locx.max nbins; locy.min locy.max nbins]);
            %Conversion from number of spikes into spikes per time
            spkrate_map = firing_map./occ_map; %element-wise division
            spkrate_map(isnan(occ_map)) = 0; %All values that are NaN in occ_map are made 0 in spkrate_map
            subplot(1,2,1) %smoothed
            %Averaging over local area with gaussian curve.
            %9 = Width of the filter in pixels (where to cut off the filter - 4 on each side of the central pixel)
            imagesc(smooth(spkrate_map,smoothfac,9,smoothfac, 9));
            subplot(1,2,2) %unsmoothed
            imagesc(spkrate_map);
            subplot(1,2,1);
            title('click on center of best place field');
            [x, y] = ginput(1);
            if isempty(x) %User presses 'Enter' because there is no centre
                centerx(i) = NaN;
                centery(i) = NaN;
            else
                centerx(i) = x;
                centery(i) = y;
            end
        end
        save(saved_centers_file, 'centerx', 'centery', 'ncells');
        
    end
    
    %Adds buffer space for shifting place fields back and forth while adding
    %them for maximum overlap
    %9x larger than actual recording area
    full_map = NaN*zeros(3*nbins, 3*nbins);
    nb_h = fix(nbins/2);
    cut_map = zeros(2*nb_h + 1, 2*nb_h + 1, ncells);
    
    for i = 1:ncells
        firing_map = histogram2(spike_pos_data{i}(:,2)', spike_pos_data{i}(:,3)', [locx.min locx.max nbins; locy.min locy.max nbins]);
        spkrate_map = firing_map./occ_map;
        %spkrate_map(isnan(occ_map)) = 0;
        %Inserts place cell map into map with buffer
        full_map(nbins+1:2*nbins,nbins+1:2*nbins) = spkrate_map;
        
        %Shifts marked centre for full_map pixel marks
        if ~isnan(centerx(i))
            new_centerx = round(centerx(i) + nbins);
            new_centery = round(centery(i) + nbins);
            
            %cut_map centres over the place field and removes extraneous
            %data far beyond the place cell. Can contain valid (inside
            %recording area) and invalid (outside recording area) data
            %(ignored during averaging).
            cut_map(:,:,i) = full_map(new_centery-nb_h:new_centery+nb_h, new_centerx-nb_h:new_centerx+nb_h);
            
        else
            cut_map(:,:,i) = -1;
            %All values in a real cut_map are positive. If a check reveals
            %a negative value in cut_map, that map will not be used in
            %map_sum.
        end
        %cut_map = array of trimmed down maps. If cut_map(i) = NaN, cell i
        %is not a place cell
        
    end
    
    %Assumes that cut_map is square.
    sideLength = size(cut_map, 1);
    
    %Size of the buffered map is double the size of cut_map
    %so regardless of the angle of rotation, the rotated buffer_map can be
    %cut down to its original size and overlaid with another buffered_map
    %without losing any data or including 0s added by rotation
    sidePlusBuffer = sideLength*2;
    
    if mod(sidePlusBuffer-sideLength,2) == 1
        %Difference is odd and cut_map will be unable to centre properly in
        %buffered_map. Adds 1 to enable centring.
        sidePlusBuffer = sidePlusBuffer + 1;
    end
    
    buffered_map = NaN*zeros(sidePlusBuffer, sidePlusBuffer, ncells);
    map_sum = zeros(sidePlusBuffer, sidePlusBuffer);
    map_count = zeros(sidePlusBuffer, sidePlusBuffer);
    
    bufferWidth = (sidePlusBuffer-sideLength)/2;
    centredLocation = (bufferWidth+1:sideLength+bufferWidth);


    
    %Template to find corners in predictable locations following rotation.
    %Corners are marked with 1s to allow them to be found once rotation occurs.
    unrotatedCorners = zeros(sidePlusBuffer, sidePlusBuffer);
    unrotatedCorners(bufferWidth+1, bufferWidth+1) = 1; %Upper left corner
    unrotatedCorners(bufferWidth+1, sideLength+bufferWidth) = 1; %Upper right corner
    unrotatedCorners(sideLength+bufferWidth, bufferWidth+1) = 1; %Lower left corner
    unrotatedCorners(sideLength+bufferWidth, sideLength+bufferWidth) = 1; %Lower right corner
    
    %Makes matrix which mimics the rotation of the buffered array for 
    %1:maxRotation degrees.
    rotatedCorners = cell(maxRotations,1);
    for i = 1:maxRotations
        rotatedCorners{i} = imrotate(unrotatedCorners, i);
    end

   
    for i = 1:ncells
        %Insert cut_maps into centre of buffered_maps
        if cut_map(1,1,i) < 0
            %NaN centre. Ignore map.
            buffered_map(:, :, i) = -1;
        else
            buffered_map(centredLocation, centredLocation, i) = cut_map(:, :, i);
        end
    end
    
    hasBaseMap = 0; %Indicates whether map_sum has been set to a
    %non-rotated map to enable rotational comparisons.
    
    for i = 1:ncells
        %Build nan_i after image rotation and cut down
        if buffered_map(1,1,i) ~= -1
            if hasBaseMap == 0
                %Skip rotation
                max_overlap = buffered_map(:,:,i);
                nan_i = isnan(max_overlap(:,:)); 
                max_overlap(nan_i) = 0;  % turn all NaN's to zero
                hasBaseMap = 1;
            else
                max_overlap = buffered_map(:,:,i);
                nan_i = isnan(max_overlap(:,:)); 
                max_overlap(nan_i) = 0;  % turn all NaN's to zero
                
                max_corr = diag(corrcoef(max_overlap,map_sum),1);
                %Correlation between max_overlap and map_sum
                angle_vs_corr = zeros(maxRotations+1, 2);
                angle_vs_corr(1,1) = 0;
                angle_vs_corr(1,2) = max_corr;

                for angle = 1:maxRotations
                    
                    [temp_map, nan_i_temp] = rotateAndPrep(buffered_map(:,:,i), ...
                        angle, sidePlusBuffer, rotatedCorners{angle});
                    
                    temp_corr = diag(corrcoef(temp_map,map_sum),1);
                    %Correlation between temp_map and map_sum
                    
                    angle_vs_corr(angle+1, 1) = angle;
                    angle_vs_corr(angle+1, 2) = temp_corr;
                    %Check whether correlation from temp_map is better than 
                    %max_overlap. If yes, set max_overlap to rotated matrix.
                    
                    if temp_corr > max_corr
                        max_overlap = temp_map;
                        nan_i = nan_i_temp;
                        max_corr = temp_corr;
                    end
                end
                
                
                f = figure;
                %figure('Position', [1400 500 1000 1000]); %CCBN
                set(f,'Position', [100 100 1000 1000]); %Home
                subplot(2,1,1);%, 'Xlim', [0 359]);                
                plot(angle_vs_corr(:,1), angle_vs_corr(:,2));
                xlim([0 359]);
                title('Angle vs Correlation');
                subplot(2,3,4);
                imagesc(map_sum);
                title('Running total');
                subplot(2,3,5);
                imagesc(buffered_map(:,:,i));
                title('Unrotated cell');
                subplot(2,3,6);
                imagesc(max_overlap);
                title('Rotated for max overlap');
                pause;
                
            end

            map_sum = map_sum + max_overlap;
            map_count = map_count + ~nan_i;
            
        end
        
    end
    
    disp(['Rat ' num2str(rat) ' Sess ' num2str(cur_sess)]);
    disp(['Total Cells: ' num2str(ncells) ' Cells With Counted Place Fields: ' num2str(sum(~isnan(centerx)))]);
    
    figure;
    avg_map{sess_i} = map_sum./map_count;
    nan_i_avg = isnan(avg_map{sess_i}(:,:)); 
    avg_map{sess_i}(nan_i_avg) = 0;
    subplot(1,2,1);
    imagesc(avg_map{sess_i});
    axis equal;
    title(['Rat ' num2str(rat) ' Sess ' num2str(cur_sess) ' Avg Field']);
    
    subplot(1,2,2);
    imagesc(smooth(avg_map{sess_i}, smoothfac, 9, smoothfac, 9));
    axis equal;
    
end

figure;
nmaps = length(avg_map);
total_sum = zeros(size(avg_map{1}));
for i = 1:length(avg_map)
    total_sum = total_sum + avg_map{i};
end

grand_avg_map = total_sum./nmaps;

imagesc(grand_avg_map);
axis equal;
title(['Grand Avg Place Field, ' num2str(nmaps) ' Sessions']);

