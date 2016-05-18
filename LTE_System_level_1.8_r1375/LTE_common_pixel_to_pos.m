function [ pos ]                      = LTE_common_pixel_to_pos( pos_pixel, roi_min, data_res)
% Converts a pixel position to position to absolute values
% (c) Josep Colom Ikuno, INTHFT, 2008
% input:    pos_pixel   ... [x,y] Pixel position (1-indexed)
%           roi_min     ... [x,y] Lower-leftmost corner of the ROI
%           data_res    ... meters/pixel, resolution of the map
% output:   pos         ... [x,y] Position to convert

pos(:,1) = (pos_pixel(:,1)-1)*data_res+roi_min(1);
pos(:,2) = (pos_pixel(:,2)-1)*data_res+roi_min(2);
