function pl_NLOS = pathloss_urbanmicro(obj,distance)
% function to evaluate the microcell LOS and NLOS pathloss based on the COST231
% Walfish-Ikegami model, see TR25.996 and COST 231 book
% (c) Martin Wrulich, INTHFT
%           distance ... actual distance in m
% output:   pl_NLOS  ... NLOS pathloss in dB

distance = distance/1000;          % Calculations are done in Km
frequency = obj.frequency/1000000; % Calculations are done in freq in MHz

% NLOS component
L_0 = 32.4 + 20*log10(distance) + 20*log10(frequency); %free space loss

delta_hmobile = obj.h_roof - obj.h_mobile;
if ((0 <= obj.phi) && (obj.phi < 35))
    L_ori = -10 + 0.354*obj.phi;
elseif ((35 <= obj.phi) && (obj.phi < 55))
    L_ori = 2.5 + 0.075*(obj.phi - 35);
elseif ((55 <= obj.phi) && (obj.phi < 90))
    L_ori = 4.0 - 0.114*(obj.phi - 55);
else
    error('Wrong wave propagation angle.');
end
L_rts = -16.9 - 10*log10(obj.w) + 10*log10(frequency) + ...
    20*log10(delta_hmobile) + L_ori; %roof-top-to-street diffracton and scatter loss

delta_hbase = obj.h_base - obj.h_roof;
if obj.h_base > obj.h_roof
    L_bsh = -18*log10(1 + delta_hbase);
else
    L_bsh = 0;
end
if obj.h_base > obj.h_roof
    k_a = 54;
elseif ((distance >= 0.5) && (obj.h_base <= obj.h_roof))
    k_a = 54 - 0.8*delta_hbase;
elseif ((distance < 0.5) && (obj.h_base <= obj.h_roof))
    k_a = 54 - 0.8*delta_hbase*distance/0.5;
end
if obj.h_base > obj.h_roof
    k_d = 18;
else
    k_d = 18 - 15*delta_hbase/obj.h_roof;
end
k_f = -4 + 1.5*(frequency/925 - 1); %for metropolitan areas
%for the other two, the COST231-WI model does not apply ;)
L_msd = L_bsh + k_a + k_d*log10(distance) + k_f*log10(frequency) - 9*log10(obj.b);
%multiple screen diffraction loss

pl_NLOS = L_0 + max(L_rts + L_msd,0); %NLOS pathloss
end