folder = './data_files/CapessoExample/pathloss/';

for i_ = 1:size(cell_pathloss_data,3)
    
    % Generate .LOS File
    filename = ['BTS' num2str(i_) '#001FL1.LOS'];
    los_file_ = fopen([folder filename], 'w');
    a_ = cell_pathloss_data(:,:,i_).*16;
    b_ = a_.';
    b_ = b_(:);
    fwrite(los_file_, b_, 'int16');
    fclose(los_file_);
    
    clear 'a_'; 
    
    if i_~=1
        copyfile([folder 'BTS1#001FL1.PAR'], [folder 'BTS' num2str(i_) '#001FL1.PAR']);
    end
end