classdef RRH < handle
    % Defines a Remote Radio Head (RRH)
    % (c) Josep Colom Ikuno, INTHFT, 2013

    properties
        parent_eNodeB      % eNodeB to which this RRH belongs
        id                 % id shared with the eNodeBs (identifies the
                           % pathloss map used also)
        site_id            % identifies the site (i.e., the shadow fading map)
        pos                % Position in [m]
        azimuth            % Sector antenna's azimuth
        antenna_type       % The type of antenna
        antenna            % Sector antenna
        nTX                % Number of antennas
        site_type = 'RRH'; 

        macroscopic_pathloss_model % Macroscopic pathloss model to be used. Empty if none is used (e.g. imported data)
    end
end
