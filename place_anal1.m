sessions = [[8068 53]
            [8068 55]];

for sess_i = 1:2;
   
   if 1
   rat = sessions(sess_i, 1);
	cur_sess = sessions(sess_i, 2);
   sess_info = session_info(rat, cur_sess, 0); % the last 1 extracts depths
   saved_centers_file = ['R' num2str(rat) 'Sess' num2str(cur_sess) 'PlaceFieldCenters.mat'];
   
   % Load the events file to obtain the epochs info
   disp(['loading Events file: ' sess_info.events_file]);
	[event_data epoch_times] = read_events4(sess_info.events_file);

   % load video tracker data
   disp(['loading PVD files: ' sess_info.pos_file1 ' & ' sess_info.pos_file2]);
   pv_data = load(sess_info.pos_file1);
   pv_data = pv_data(logical([diff(pv_data(:,1))>0; 1]),:);
   pv_data = pv_data(:,1:3); % discard some data to save memory
   disp('done');	

   tfile_list = FindFiles('*.t', 'StartingDirectory', sess_info.tfile_path);
   ts_data = LoadSpikes(tfile_list);
   keep_celli = logical(ones(1,length(ts_data)));

   for i = 1:length(ts_data)

      % extract tetrode and cell number from filename
      [path name ext] = fileparts(tfile_list{i});
      tet_name_i = findstr(name, 'TT')+2;
      [tt_id cell_id] = strread(name, 'TT%dnsc_%d');   
      if tt_id>6  % block tetrodes 7-12, which are in frontal cortex
         keep_celli(i) = 0;
      end;
      cell_info(i).tet_id = tt_id;
      cell_info(i).cell_num = cell_id;
      cell_info(i).filename = tfile_list{i};

      times = Data(ts_data{i}); % extract data from the silly ts object
      times = times*100;  % convert to microseconds

      spike_data{i} = times;

      % count spikes in each epoch
      epoch_counts = NaN*zeros(1,3);
      for j = 1:size(epoch_times,1)
         epoch_counts(j) = sum(times>epoch_times(j,1) & times<epoch_times(j,2));
      end;
      cell_info(i).counts = epoch_counts;
      if cell_info(i).counts(2)<50
         keep_celli(i) = 0;
         disp(['Cell TT' num2str(cell_info(i).tet_id) '_' num2str(cell_info(i).cell_num) ' rejected: ' num2str(cell_info(i).counts(2)) ' spikes']);
      end;

      %posit = interp1(pv_data(:,1), [pv_data(:,2) pv_data(:,3)], times);
      %spike_data{i} = [times posit];
      %spike_dist{i} = interp1(pv_data(:,1), cumdist, times);
      disp(['Loaded: ' tfile_list{i}]);

      posit = interp1(pv_data(:,1), [pv_data(:,2) pv_data(:,3)], times);
      spike_pos_data{i} = [times posit];

   end;  % if exist(sess_info.tfile_path, 'dir')
   spike_pos_data = spike_pos_data(keep_celli);
   cell_info = cell_info(keep_celli);
   spike_data = spike_data(keep_celli);
   ncells = length(spike_data);
   %disp(['Cells after screening: ' num2str(ncells)]);
   
   end;  % if 0
   
   
   nbins = 50;
   smoothfac = .75;
   locx.min = min(pv_data(:,2));
   locx.max = max(pv_data(:,2));
   locy.min = min(pv_data(:,3));
   locy.max = max(pv_data(:,3));
   occ_map = histogram2(pv_data(:,2)', pv_data(:,3)', [locx.min locx.max nbins; locy.min locy.max nbins]);
   occ_map(occ_map==0) = NaN;
   
   if exist(saved_centers_file)
      cur_ncells = ncells;
      load(saved_centers_file)
      if cur_ncells ~= ncells
         disp(['Cells in stored file ' saved_centers_file ' dont match current data.  Recompute centers.']); 
         return;
      end;
   else
      centerx = ones(1, length(spike_data));
      centery = ones(1, length(spike_data));
      for i = 1:ncells
         firing_map = histogram2(spike_pos_data{i}(:,2)', spike_pos_data{i}(:,3)', [locx.min locx.max nbins; locy.min locy.max nbins]);
         spkrate_map = firing_map./occ_map;
         spkrate_map(isnan(occ_map)) = 0;
         subplot(1,2,1)
         imagesc(smooth(spkrate_map,smoothfac,9,smoothfac, 9));
         subplot(1,2,2)
         imagesc(spkrate_map);
         subplot(1,2,1);
         title('click on center of best place field');
         [x, y] = ginput(1);
         if isempty(x)
            centerx(i) = NaN;
            centery(i) = NaN;
         else
            centerx(i) = x;
            centery(i) = y;
         end;
      end;
      save(saved_centers_file, 'centerx', 'centery', 'ncells');
      
   end; % if exist(saved_centers_file)
   
   full_map = NaN*zeros(3*nbins, 3*nbins);
   nb_h = int8(fix(nbins/2));
   map_sum = zeros(2*nb_h + 1, 2*nb_h + 1);
   map_count = zeros(2*nb_h + 1, 2*nb_h + 1);
   for i = 1:ncells
      firing_map = histogram2(spike_pos_data{i}(:,2)', spike_pos_data{i}(:,3)', [locx.min locx.max nbins; locy.min locy.max nbins]);
      spkrate_map = firing_map./occ_map;
      %spkrate_map(isnan(occ_map)) = 0;
      full_map(nbins+1:2*nbins,nbins+1:2*nbins) = spkrate_map; 

      %imagesc(smooth(full_map,smoothfac,9,smoothfac, 9));
      %hold on;
      
      if ~isnan(centerx(i))
         new_centerx = int8(centerx(i) + nbins);
         new_centery = int8(centery(i) + nbins);

%          subplot(1,2,1);
%          imagesc(full_map);
%          hold on;
%          plot(new_centerx, new_centery,'xw', 'markersize', 10); 
%          hold off;
         
%          subplot(1,2,2);
          cut_map = full_map(new_centery-nb_h:new_centery+nb_h, new_centerx-nb_h:new_centerx+nb_h);
%          imagesc(cut_map);
         %pause;
      
         nan_i = isnan(cut_map);
         temp_map = cut_map;
         temp_map(nan_i) = 0;  % turn all NaN's to zero
         map_sum = map_sum + temp_map;
         map_count = map_count + ~nan_i;
         
      end;
   
   end;  % for i = 1:ncells
   
   disp(['Rat ' num2str(rat) ' Sess ' num2str(cur_sess)]);
   disp(['Total Cells: ' num2str(ncells) ' Cells With Counted Place Fields: ' num2str(sum(~isnan(centerx)))]);
   
   figure;
   avg_map{sess_i} = map_sum./map_count;
   subplot(1,2,1);
   imagesc(avg_map{sess_i});
   axis equal;
   title(['Rat ' num2str(rat) ' Sess ' num2str(cur_sess) ' Avg Field']);
   
   subplot(1,2,2);
   imagesc(smooth(avg_map{sess_i}, smoothfac, 9, smoothfac, 9));
   axis equal;

end;

figure;
nmaps = length(avg_map);
map_sum = zeros(size(avg_map{1}));
for i = 1:length(avg_map)
   map_sum = map_sum + avg_map{i};
end;

grand_avg_map = map_sum./nmaps;

imagesc(grand_avg_map);
axis equal;
title(['Grand Avg Place Field, ' num2str(nmaps) ' Sessions']);

    