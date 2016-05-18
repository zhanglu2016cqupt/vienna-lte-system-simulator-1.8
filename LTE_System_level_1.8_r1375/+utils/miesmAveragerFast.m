classdef miesmAveragerFast < utils.sinrAverager
    % Implements a MIESM averager. Needs a precalculated table containing
    % the BICM capacity for each of the constellations used, which are here
    % classified according to their modulation orders. Alternative
    % implementation to the original MIESM SINR averager that works faster
    % This MIESM averager DOES support calibration
    % (c) Josep Colom Ikuno, INTHFT, 2010

   properties
       length
       BICM_matrix_x
       BICM_matrix_y
       BICM_matrix_inv_x
       BICM_matrix_inv_y
       resolution = 0.01;
       resolution_inv
       multipliers
       SINR_range
       BICM_range
       nCQIs
       betas
       betas_dB
   end

   methods
       function obj = miesmAveragerFast(LTE_config,mat_file_to_load,betas,varargin)
           % Create the inverse I curve (the same as with the previous implementation)
           load(mat_file_to_load);
           
           % Upsample the BICM curves
           SNR_vector = BICM_capacity_tables(1).SNR(1):obj.resolution:BICM_capacity_tables(1).SNR(end);
           
           CQI_range  = LTE_common_get_CQI_params(LTE_config,'range');
           CQI_tables = LTE_common_get_CQI_params(LTE_config,CQI_range(1):CQI_range(2));
           n_CQIs     = length(CQI_tables);
           mod_orders = [CQI_tables.modulation_order];
           
           for i_=1:length(BICM_capacity_tables)
               % Overwrite loaded values with upsampled ones
               BICM_capacity_tables(i_).I   = interp1(BICM_capacity_tables(i_).SNR,BICM_capacity_tables(i_).I,SNR_vector);
               BICM_capacity_tables(i_).SNR = SNR_vector;
               
               % Generate the inverse mapping
               BICM_capacity_tables(i_).I(1) = 0;
               BICM_capacity_tables(i_).I(end) = ceil(BICM_capacity_tables(i_).I(end));
               [b,m,n]=unique(BICM_capacity_tables(i_).I,'first');
               
               % To have only unique values
               BICM_capacity_tables(i_).I_inv = b;
               BICM_capacity_tables(i_).SNR_inv = BICM_capacity_tables(i_).SNR(m);
           end
           
           % Assume that the SNR is the same one for all data. It should be
           % automaticallythe case if the BICM capacity script of this
           % simulator was used.
           min_C = 0;
           max_C = max([BICM_capacity_tables.I]);
           BICM_matrix_x     = SNR_vector;
           BICM_matrix_y     = zeros(length(BICM_matrix_x),n_CQIs);
           BICM_matrix_inv_x = linspace(min_C,max_C,length(BICM_matrix_x));
           BICM_matrix_inv_y = zeros(length(BICM_matrix_inv_x),n_CQIs);

           obj.nCQIs = length(CQI_range(1):CQI_range(2));
           if length(betas)~=obj.nCQIs
               error('length of beta calibration parameters must be %d',obj.nCQIs);
           end
           
           betas        = betas(:);        % Store in column format
           obj.betas    = betas;
           obj.betas_dB = 10*log10(betas); % Precalculate also the value in dB
           
           for cqi_ = CQI_range(1):CQI_range(2)
               m_j          = CQI_tables(cqi_).modulation_order;
               struct_index = find([BICM_capacity_tables.m_j]==m_j,1,'first');
               BICM_matrix_y(:,cqi_)     = BICM_capacity_tables(struct_index).I;
               BICM_matrix_inv_y(:,cqi_) = interp1(BICM_capacity_tables(struct_index).I_inv,BICM_capacity_tables(struct_index).SNR_inv,BICM_matrix_inv_x);
               
               % Small fix not to have a NaN at the maximum capacity point
               max_capacity = max(BICM_capacity_tables(struct_index).I);
               [C,I] = min(abs(BICM_matrix_inv_x-max_capacity));
               BICM_matrix_inv_y(I,cqi_) = interp1(BICM_matrix_inv_x(1:(I-1)),BICM_matrix_inv_y(1:(I-1),cqi_),BICM_matrix_inv_x(I),'linear','extrap');
           end
           
           obj.length            = length(BICM_matrix_x);
           obj.BICM_matrix_x     = BICM_matrix_x;
           obj.BICM_matrix_y     = BICM_matrix_y;
           obj.BICM_matrix_inv_x = BICM_matrix_inv_x;
           obj.BICM_matrix_inv_y = BICM_matrix_inv_y;
           obj.multipliers       = reshape(0:length(BICM_matrix_x):length(BICM_matrix_x)*(n_CQIs-1),[],1);
           
           % Assume equally-spaced SNRs
           obj.resolution     = (BICM_matrix_x(2)-BICM_matrix_x(1));
           obj.resolution_inv = (BICM_matrix_inv_x(2)-BICM_matrix_inv_x(1));
           
           obj.SINR_range = [min(BICM_matrix_x) max(BICM_matrix_x)];
           obj.BICM_range = [min(BICM_matrix_inv_x) max(BICM_matrix_inv_x)];
           
           % Some (optional) plotting
           if isempty(varargin)
               plot_capacity = false;
           else
               plot_capacity = varargin{1};
           end
           if plot_capacity
               for m_j_idx=1:length(BICM_capacity_tables)
                   BICM_capacity_tables_all(m_j_idx,:) = BICM_capacity_tables(m_j_idx).I;
                   displaynames{m_j_idx} = sprintf('BICM capacity, %d-QAM',2^BICM_capacity_tables(m_j_idx).m_j);
               end
               figure;
               % Assume that all BICM capacities are calcualted over the same
               % SNR range (so I can plot it like this and Matlab automatically
               % puts the colors there)
               plot(BICM_capacity_tables(m_j_idx).SNR,BICM_capacity_tables_all);
               legend(displaynames,'Location','Best');
               xlabel('SNR [dB]');
               ylabel('BICM Capacity I_{m_j}(\gamma)');
               title('BICM capacity');
               grid on;
           end
       end
       
       % Average the give SINR vector. varargin contains the following:
       %   - varargin{1} = MCS -> values in the range 0:15
       % ALL INPUT VALUES IN LINEAR unless an extra parameter is set to "true"
       % Allows for calculations with multiple beta values
       function [effective_SINR_dB, Is_mean, Is_std] = average(obj,SINR_vector,MCSs,varargin)   
           if isempty(varargin)
               input_in_dB  = false;
           else
               input_in_dB = varargin{1};
           end
           
           % Put SINR in a row vector if it is a column vector
           if size(SINR_vector,1)~=1 && size(SINR_vector,2)==1
               SINR_vector = reshape(SINR_vector,1,[]);
           end
           
           % Data needed beforehand
           SNR_length      = length(obj.BICM_matrix_x);
           multipliers     = obj.multipliers(MCSs);
           multipliers_mat = multipliers(:,ones(length(SINR_vector),1));
           
           if input_in_dB
               SINR_vector_dB = SINR_vector;
           else
               SINR_vector_dB = 10*log10(SINR_vector);
           end

           MC_betas_vect_dB    = obj.betas_dB(MCSs);
           SINR_vector_mat_dB  = SINR_vector_dB(ones(length(MCSs),1),:);
           MCS_betas_dB        = MC_betas_vect_dB(:,ones(length(SINR_vector_dB),1));
           SINR_vector_mat_log = SINR_vector_mat_dB - MCS_betas_dB;         
           
           SNR_idxs_mat        = round((SINR_vector_mat_log-obj.SINR_range(1))/obj.resolution + 1);
           SNR_idxs_mat(SNR_idxs_mat<1) = 1;
           SNR_idxs_mat(SNR_idxs_mat>SNR_length) = SNR_length;

           % Convert to BICM capacity
           Is_mat  = obj.BICM_matrix_y(SNR_idxs_mat+multipliers_mat);
           Is_mean = mean(Is_mat,2);
           Is_std  = std(Is_mat,0,2);
           
           % Inverse mapping
           Is_idxs = round((Is_mean-obj.BICM_range(1))/obj.resolution_inv + 1);
           Is_idxs(Is_idxs<1) = 1; % Safeguard agains negative indices
           obj_length = obj.length;
           Is_idxs(Is_idxs>obj_length) = obj_length;
           effective_SINR_dB = obj.BICM_matrix_inv_y(Is_idxs+multipliers)+MC_betas_vect_dB; % Version with scaling
           if  ~sum(SINR_vector_dB < obj.SINR_range(end)) % all SINRs above BICM range --> use linear averaging instead
               effective_SINR_dB(:) = 10*log10(mean(10.^(SINR_vector_dB/10)));
           end
       end
       
       % Average the give MI vector which contains N ROWS from different layers to one codeword. varargin contains the following:
       %   - varargin{1} = MCS -> value (SINGLE ONE!) in the range 0:15
       function [effective_SINR_dB] = average_codeword(obj,Is_input,MCSs)
           % Data needed beforehand
           rank_number      = size(Is_input,2);
           multipliers      = obj.multipliers(MCSs);
           MC_betas_vect_dB = obj.betas_dB(MCSs);
           Is_mean          = mean(Is_input,1);
           
           % Inverse mapping
           Is_idxs = round((Is_mean-obj.BICM_range(1))/obj.resolution_inv + 1);
           Is_idxs(Is_idxs<1) = 1; % Safeguard agains negative indices
           obj_length = obj.length;
           Is_idxs(Is_idxs>obj_length) = obj_length;
           Is_idxs_finite = isfinite(Is_idxs); % It may happen that in FFR some values are NaN
           effective_SINR_dB = -Inf(size(Is_idxs));
           effective_SINR_dB(Is_idxs_finite) = obj.BICM_matrix_inv_y(Is_idxs(Is_idxs_finite)+multipliers)+MC_betas_vect_dB; % Version with scaling
       end
       
       % Average the give MI vector for the RI calculation of SM. varargin contains the following:
       %   - varargin{1} = MCS -> values in the range 0:15
       function [effective_SINR_dB] = average_for_RI(obj,Is_mean,MCSs)
           if length(MCSs)~=size(Is_mean,1)
               error('Dimensions do not agree');
           end
           
           % Data needed beforehand
           rank_number      = size(Is_mean,2);
           multipliers      = obj.multipliers(MCSs);
           MC_betas_vect_dB = obj.betas_dB(MCSs);
           
           % Inverse mapping
           Is_idxs = round((Is_mean-obj.BICM_range(1))/obj.resolution_inv + 1);
           Is_idxs(Is_idxs<1) = 1; % Safeguard agains negative indices
           obj_length = obj.length;
           Is_idxs(Is_idxs>obj_length) = obj_length;
           effective_SINR_dB = obj.BICM_matrix_inv_y(Is_idxs+multipliers(:,ones(1,rank_number)))+MC_betas_vect_dB(:,ones(1,rank_number)); % Version with scaling
       end
       
       function [Is_sum Is_sorted_idxs] = order_SINRs_per_MI(obj,Is_mat)
           disp(); % TO DO
       end
       
       % Converts a matrix of SINRs (in dB!) to MI (BICM Capacity)
       % MCSs is a vector of MCSs to be used (i.e. different beta values)
       function I_MCSs = SINR_to_I(obj,SINR_dB,MCSs)
           SINR_dims             = size(SINR_dB);
           SINR_dims_numel       = length(SINR_dims);
           SINR_dB_vect          = reshape(SINR_dB(:),1,[]);
           SINR_dB_vect_idxs     = isfinite(SINR_dB_vect);
           SINR_dB_vect_filtered = SINR_dB_vect(SINR_dB_vect_idxs);
           I_MCSs                = zeros([SINR_dims length(MCSs)]);
           
           % Data needed beforehand
           SNR_length      = length(obj.BICM_matrix_x);
           multipliers     = obj.multipliers(MCSs);
           multipliers_mat = multipliers(:,ones(length(SINR_dB_vect_filtered),1));
           
           % Convert to BICM capacity
           MC_betas_vect_dB    = obj.betas_dB(MCSs);
           SINR_vector_mat_dB  = SINR_dB_vect_filtered(ones(length(MCSs),1),:);
           MCS_betas_dB        = MC_betas_vect_dB(:,ones(length(SINR_dB_vect_filtered),1));
           SINR_vector_mat_log = SINR_vector_mat_dB - MCS_betas_dB;
           
           SNR_idxs_mat        = round((SINR_vector_mat_log-obj.SINR_range(1))/obj.resolution + 1);
           SNR_idxs_mat(SNR_idxs_mat<1) = 1;
           SNR_idxs_mat(SNR_idxs_mat>SNR_length) = SNR_length;
           
           % BICM capacity table lookup
           Is_mat  = obj.BICM_matrix_y(SNR_idxs_mat+multipliers_mat);
           
           I_out_vect = NaN(size(Is_mat,1),length(SINR_dB_vect));
           I_out_vect(:,SINR_dB_vect_idxs) = Is_mat;
           
           % Fill in the output
           for MCS_idx = 1:length(MCSs)
               current_MCS_Is = reshape(I_out_vect(MCS_idx,:),SINR_dims);
               switch SINR_dims_numel
                   case 2
                       I_MCSs(:,:,MCS_idx) = current_MCS_Is;
                   case 3
                       I_MCSs(:,:,:,MCS_idx) = current_MCS_Is;
                   otherwise
                       error('I still didnt manage to learn how to handle arbitrary dimension numbers... sorry!');
               end
           end
       end
       
       % Average a vector of MI in the given dimension and convert back to
       % the SINR domain. The last dimension is assumed to be of length equal to MCSs
       function effective_SINR_dB_output = average_I(obj,Is,dim,MCSs)
           Is_mean          = mean(Is,dim);
           Is_mean_dims     = size(Is_mean);
           if Is_mean_dims(end)~=length(MCSs)
               error('dimensions do not match! length(MCSs) has to be equal to the size of the last dimension');
           end
           Is_mean_vect     = Is_mean(:);
           multipliers      = kron(obj.multipliers(MCSs),ones(prod(Is_mean_dims(1:(end-1))),1));
           MC_betas_vect_dB = kron(obj.betas_dB(MCSs)   ,ones(prod(Is_mean_dims(1:(end-1))),1));
           
           % Inverse mapping
           Is_idxs = round((Is_mean_vect-obj.BICM_range(1))/obj.resolution_inv + 1);
           Is_idxs(Is_idxs<1) = 1; % Safeguard agains negative indices
           obj_length = obj.length;
           Is_idxs(Is_idxs>obj_length) = obj_length;
           effective_SINR_dB = obj.BICM_matrix_inv_y(Is_idxs+multipliers)+MC_betas_vect_dB; % Version with scaling
           
           effective_SINR_dB_output = reshape(effective_SINR_dB,Is_mean_dims);
       end
   end
end 
