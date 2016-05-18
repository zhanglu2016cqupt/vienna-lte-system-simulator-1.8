% WINNER Phase II channel model
% Version 1.3, Aug. 15, 2008
%
%% Channel model functions
%   wim                          - WINNER channel model (D1.1.2)
%     wim_core                   - Computes the big channel matrix generation loop i.e. the formula in [2, Eq. 4.14, 4.17, 4.19]. 
%   scenpartables                - Set WIM parameters for WINNER scenarios
%   pathloss                     - Pathloss models for 2GHz and 5GHz 
%
%% Helper functions for model Initialization
%   wimparset                    - Model parameter configuration for WIM
%   layoutparset                 - Layout parameter configuration for WIM
%
%% Determination of link-level parameters
%   layout2link                  - Computes and converts layout to link parameters. The following helper functions are used:
%     StationDistXY              - Distance between stations in XY plane
%     StationDirectionXY         - Link angular direction in Global-Coordinate-System, w.r.t Y-axis 
%     StationVelocityXY          - MS velocity 
%
%% Utility functions
%   StationNumElements           - Retrieve number of array elements
%   NTlayout                     - Visualisation of network layout
%
%% Generation of structural parameters
%   generate_bulk_par            - Generation of RANDOM WIM bulk parameters
%     LScorrelation.m            - Correlation of Large-Scale-Parameters (LSP)
%     LOSprobability.m           - Probability of LoS condition
%     offset_matrix_generation.m - Generation of fixed angle offsets
%     ScenarioMapping.m          - Retrive WINNER scnario label
%     struct_generation.m        - Manipulation with bulk_parameters
%   fixedAoas.m                  - FIXED AoAs for CDL model
%   fixedAods.m                  - FIXED AoDs for CDL model
%   fixedPdp.m                   - FIXED PDP for CDL model
%
%% 3D Antenna Array Model functions
%   AntennaArray.m               - 3D-AA model construction
%     ArrayPreprocess.m          - construction w. element field-pattern rotatation in ACS
%   AntennaResponse.m            - Computation of array response
%
%   Aperture_Calc.p              - EADF representation of radiation pattern
%     BP2Aperture.m              - Interface toward Aperture_Calc.p for 2D field pattern
%     BP2Aperture1D.m            - Interface toward Aperture_Calc.p for 1D field pattern
%     G_Calc3D_simple.p          - Calculation for 2D field pattern (3D-AA)
%     G_Calc1D_simple.p          - Calculation for 1D field pattern
%   antenna_pol_vect.m           - Computation of polarization vectors
%
%   mycart2sph.m                 - Transfmation of carthesian into spherical coordinates
%   mysph2cart.m                 - Transfmation of spherical into carthesian coordinates
%   rotate_vector.m              - 3D rotation
%   unrotate_vector.m            - 3D inverse rotation
%
%% Miscellaneous functions
%   cas                          - Circular angle spread (3GPP TR 25.996)
%   ds                           - RMS delay spread 
%   dipole                       - Field pattern of half wavelength dipole
%
%% Examples
%   example_EADF_approx.m        - Representation of radiation pattern
%   example_syntetic_arrays.m    - Array struct generation
%   example_channel_matrix.m     - Generation of channel matrix