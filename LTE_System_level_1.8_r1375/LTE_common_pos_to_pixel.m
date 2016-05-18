function [ pos_pixel pos_pixel_exact] = LTE_common_pos_to_pixel( pos,       roi_min, data_res)
% Converts a position in absolute values to a pixel position
% (c) Josep Colom Ikuno, INTHFT, 2008
% input:    pos             ... [x,y] Position to convert
%           roi_min         ... [x,y] Lower-leftmost corner of the ROI
%           data_res        ... meters/pixel, resolution of the map
% output:   pos_pixel       ... [x,y] Pixel position (1-indexed)
%           pos_pixel_exact ... [x,y] Pixel position, exact value. Can be
%                               used to perform interpolation. As with
%                               pos_pixel, it is one-indexed.
%
% NOTE: this code is replicated in the get_pathloss functions for both the
% macroscopic pathloss and shadow fading. Since those functions are called
% thousands of times, it was faster to include this couple of lines there
% instead of calling the function.

pos_pixel(:,1) = floor((pos(:,1)-roi_min(1))/data_res)+1;
pos_pixel(:,2) = floor((pos(:,2)-roi_min(2))/data_res)+1;

pos_pixel_exact(:,1) = ((pos(:,1)-roi_min(1))/data_res)+1;
pos_pixel_exact(:,2) = ((pos(:,2)-roi_min(2))/data_res)+1;

