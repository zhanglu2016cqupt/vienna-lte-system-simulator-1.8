% Take pathlossmaps as reference size of zero matrix
% Rows and cols exchanged because Matrix is transposed after readout !
myDefaultDTM = zeros(473, 412);
myfile_ = fopen('./data_files/CapessoExample/dtm/exampleDTM.BIL', 'w');
fwrite(myfile_, myDefaultDTM, 'int16');
fclose(myfile_);