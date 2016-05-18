function process_UE_feedback_at_eNodeB( eNodeBs, config, noise_power)
% Feedback calculation at the eNodeB side
% (c) Josep Colom Ikuno, jcolom@nt.tuwien.ac.at ITC, 2013
% www.nt.tuwien.ac.at

all_feedback = [eNodeBs.last_received_feedback];

noise_W_half_RB = noise_power/2;

tx_mode = -1;
% Assumes that all of the UEs employ the same transmit mode
for i=1:length(all_feedback)
    if ~isempty(all_feedback(i).tx_mode)
        tx_mode = all_feedback(i).tx_mode(1);
        break;
    end
end
if tx_mode==-1
   error('No feedback received at any eNodeB'); 
end

if tx_mode==0
    % Feedback still did not arrive: do nothing. It will be dealt with at
    % the scheduler with default values.
elseif tx_mode<100
    % In these "standard" TX modes, the needed feedback is already there
else
    switch tx_mode
        case 400
            for b_=1:length(all_feedback)
                % CLSM
                CSI_feedbacks = all_feedback(b_).CSI;
                nUEs          = length(CSI_feedbacks);
                nRBs          = eNodeBs(b_).RB_grid.n_RB;
                
                precoding_codebook = [eNodeBs(b_).scheduler.Rel8_codebook{:,eNodeBs(b_).total_nTX,4}];
                eNodeBs(b_).last_received_feedback.CQI = zeros(nRBs,2,nUEs);
                eNodeBs(b_).last_received_feedback.RI  = zeros(1,nRBs); % Number of layers
                
                % Fill in CQI and RI values (needed by the scheduler to
                % know what MCSs to apply and how many layers are there
                for u_=1:length(CSI_feedbacks)
                    [PMI_feedback,TB_SINR_feedback_dB] = feedback_calculation.CLSM_PMI_feedback(...
                        precoding_codebook,...
                        CSI_feedbacks(u_).H_0, CSI_feedbacks(u_).TX_W_RB_0, noise_W_half_RB,...
                        CSI_feedbacks(u_).H_i, CSI_feedbacks(u_).TX_W_RB_i, CSI_feedbacks(u_).W_i, CSI_feedbacks(u_).ICI_power );

                    feedback = feedback_calculation.CQI_and_rank(...
                        config,...
                        4,...
                        TB_SINR_feedback_dB,...
                        nRBs,... % The PMI feedback won't be empty, so no need for this variable
                        PMI_feedback,...
                        [],... % DL_signaling only needed for the adaptive RI case
                        eNodeBs(b_).scheduler.SINR_averager,...
                        eNodeBs(b_).scheduler.CQI_mapper,...
                        eNodeBs(b_).unquantized_CQI_feedback,...
                        false); % No adaptive RI
                    % Fill in CQI and precoder, and RI (number of layers) in the feedback structure
                    eNodeBs(b_).last_received_feedback.CQI(:,1:size(feedback.CQI,1),u_) = feedback.CQI';
                    eNodeBs(b_).last_received_feedback.precoder(u_).W = ...
                        precoding_codebook(feedback.RI).W(:,:,feedback.PMI);
                    eNodeBs(b_).last_received_feedback.RI(u_) = feedback.RI;
                end
            end
        otherwise
            error('eNodeB feedback processing for TX mode %d not implemented',tx_mode);
    end
    
end

end

