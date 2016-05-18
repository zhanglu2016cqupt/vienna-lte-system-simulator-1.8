function CQI_params = LTE_common_get_CQI_params(LTE_config,CQI)
% Given a CQI value, this function returns the associated CQI parameteres,
% as specified in TS36.213, table 7.2.3-1 (4-bit CQI table)
% [CQI_params] = LTE_common_get_CQI_params(CQI)
% Author: Josep Colom Ikuno, josep.colom@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input:    CQI         ... [1 x N]double - CQI value
% 
% output:   CQI_params  ... [1 x N]struct settings associated to
%                           that CQI value.
%
% date of creation: 2008/08/11
% last changes:
%    13/02/2009: improved function call. Now CQI values are stored in
%    LTE_params and the function can be called with a CQI vector as an
%    argument.

% CQI table
%
% CQI	Mod     rate (x 1024)	rate (1/n)
% 1     QPSK    78              13.13
% 2     QPSK	120             8.53
% 3     QPSK	193             5.31
% 4     QPSK	308             3.32
% 5     QPSK	449             2.28
% 6     QPSK	602             1.70
% 7     16QAM	378             2.71
% 8     16QAM	490             2.09
% 9     16QAM	616             1.66
% 10	64QAM	466             2.20
% 11	64QAM	567             1.81
% 12	64QAM	666             1.54
% 13	64QAM	772             1.33
% 14	64QAM	873             1.17
% 15	64QAM	948             1.08

% Some extra code would be needed if CQIs<1 would be introduced, which is
% not the case now.

if strcmp(CQI,'range')
    CQI_params = [1 15];
else
    CQI_params = LTE_config.CQI_params(CQI);
end
