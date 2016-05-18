classdef simulationResults < handle
    % Needed to import LL simulation results
    % Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
    % (c) 2009 by INTHFT
    % www.nt.tuwien.ac.at

    properties
        cell_specific  % Cell specific traces
        UE_specific    % UE-specific traces
        nUE
        maxStreams
        SNR_vector     % SNR vector used for this simulation
    end

    methods
        % Class constructor. Preallocate all necessary space
        function obj = simulationResults(eNodeB_count,UE_count,N_subframes,SNR_vector,maxStreams,nRx,nTx,trace_subcarrier_SNR,Ntot)

            obj.nUE = UE_count;
            obj.maxStreams = maxStreams;

            SNR_vector_length = length(SNR_vector);

            % Preallocate cell specific traces
            obj.cell_specific = results.cellSpecificTraces(1,1,1,1,1,Ntot,false);
            for b_ = 1:eNodeB_count
                obj.cell_specific(b_) = results.cellSpecificTraces(N_subframes,SNR_vector_length,maxStreams,nRx,nTx,Ntot,trace_subcarrier_SNR);
            end

            % Preallocate UE-specific traces
            obj.UE_specific = results.ueSpecificTraces(1,1,1);
            for u_ = 1:UE_count
                obj.UE_specific(u_) = results.ueSpecificTraces(N_subframes,SNR_vector_length,maxStreams);
            end
        end

        % Process results from this TTI. We will assume that only 1 eNodeB is in the simulation resutls file
        function process_TTI_results(obj,BS_output,UE_output,subframe_i,SNR_i)
            % Loop over all UEs and streams
            for uu = 1:obj.nUE
                for stream_i = 1:BS_output.UE_signaling(uu).MCS_and_scheduling.nCodewords

                    % Update ACK (UE traces) and biterrors (cell traces)
                    obj.UE_specific(uu).ACK(subframe_i,SNR_i,stream_i)             = UE_output(uu).ACK(stream_i);
                    obj.UE_specific(uu).rv_idx(subframe_i,SNR_i,stream_i)          = UE_output(uu).rv_idx(stream_i);
                    
                    % Update biterrors (UE traces)
                    obj.UE_specific(uu).biterrors_coded(subframe_i,SNR_i,stream_i)   = sum(abs(UE_output(uu).rx_data_bits{stream_i}  - BS_output.genie(uu).data_bits{stream_i}));
                    obj.UE_specific(uu).biterrors_uncoded(subframe_i,SNR_i,stream_i) = sum(abs(UE_output(uu).rx_coded_bits{stream_i} - BS_output.genie(uu).sent_bits{stream_i}));

                    % Update blocksize (UE traces)
                    obj.UE_specific(uu).blocksize_coded(subframe_i,SNR_i,stream_i)   = length(UE_output(uu).rx_data_bits{stream_i});
                    obj.UE_specific(uu).blocksize_uncoded(subframe_i,SNR_i,stream_i) = length(UE_output(uu).rx_coded_bits{stream_i});

                    % Update FER and throughput (UE traces)
                    
                    % Coded
                    if UE_output(uu).ACK(stream_i)
                        obj.UE_specific(uu).throughput_coded(subframe_i,SNR_i,stream_i) = length(UE_output(uu).rx_data_bits{stream_i});
                        obj.UE_specific(uu).FER_coded(subframe_i,SNR_i,stream_i) = 0;
                    else
                        obj.UE_specific(uu).FER_coded(subframe_i,SNR_i,stream_i) = 1;
                    end
                    
                    % Uncoded
                    if (obj.UE_specific(uu).biterrors_uncoded(subframe_i,SNR_i,stream_i)==0)
                        obj.UE_specific(uu).throughput_uncoded(subframe_i,SNR_i,stream_i) = length(UE_output(uu).rx_coded_bits{stream_i});
                        obj.UE_specific(uu).FER_uncoded(subframe_i,SNR_i,stream_i)        = 0;
                    else
                        obj.UE_specific(uu).FER_uncoded(subframe_i,SNR_i,stream_i) = 1;
                    end

                    % Update what codewords were used (valid positions in the traces)
                    obj.UE_specific(uu).used_codewords(subframe_i,SNR_i,stream_i) = 1;
                end

                % Update cell coded and uncoded bit errors
                obj.cell_specific.biterrors_coded(subframe_i,SNR_i,:)   = obj.cell_specific.biterrors_coded(subframe_i,SNR_i,:)   + obj.UE_specific(uu).biterrors_coded(subframe_i,SNR_i,:);
                obj.cell_specific.biterrors_uncoded(subframe_i,SNR_i,:) = obj.cell_specific.biterrors_uncoded(subframe_i,SNR_i,:) + obj.UE_specific(uu).biterrors_uncoded(subframe_i,SNR_i,:);

                % Update blocksize (cell traces)
                obj.cell_specific.blocksize_coded(subframe_i,SNR_i,:)   = obj.cell_specific.blocksize_coded(subframe_i,SNR_i,stream_i)   + obj.UE_specific(uu).blocksize_coded(subframe_i,SNR_i,:);
                obj.cell_specific.blocksize_uncoded(subframe_i,SNR_i,:) = obj.cell_specific.blocksize_uncoded(subframe_i,SNR_i,stream_i) + obj.UE_specific(uu).blocksize_uncoded(subframe_i,SNR_i,:);

                % Update FER and throughput (cell traces)
                obj.cell_specific.throughput_coded(subframe_i,SNR_i,:)  = obj.cell_specific.throughput_coded(subframe_i,SNR_i,stream_i) + obj.UE_specific(uu).throughput_coded(subframe_i,SNR_i,:);
                obj.cell_specific.throughput_uncoded(subframe_i,SNR_i,:)= obj.cell_specific.throughput_uncoded(subframe_i,SNR_i,stream_i) + obj.UE_specific(uu).throughput_uncoded(subframe_i,SNR_i,:);
                obj.cell_specific.FER_coded(subframe_i,SNR_i,:)         = obj.cell_specific.FER_coded(subframe_i,SNR_i,:)   + uint16(obj.UE_specific(uu).FER_coded(subframe_i,SNR_i,:));
                obj.cell_specific.FER_uncoded(subframe_i,SNR_i,:)       = obj.cell_specific.FER_uncoded(subframe_i,SNR_i,:) + uint16(obj.UE_specific(uu).FER_uncoded(subframe_i,SNR_i,:));
                
                % The number of codewords is the maximum (bitwise OR) of the used codewords matrix
                obj.cell_specific.used_codewords(subframe_i,SNR_i,:) = obj.cell_specific.used_codewords(subframe_i,SNR_i,:) + uint16(obj.UE_specific(uu).used_codewords(subframe_i,SNR_i,:));
                
                channel_error_help = zeros(size(obj.cell_specific.channel_error(subframe_i,SNR_i,:,:)));
                channel_error_help(1,1,:,:) = UE_output(uu).channel_estimation_error;
                obj.cell_specific.channel_error(subframe_i,SNR_i,:,:) = obj.cell_specific.channel_error(subframe_i,SNR_i,:,:) + channel_error_help;
            end
            obj.cell_specific.channel_error(subframe_i,SNR_i,:,:) = obj.cell_specific.channel_error(subframe_i,SNR_i,:,:)/obj.nUE;
            
        end

        % Calculate simulations aggregates.
        function calculate_sim_aggregates(obj,elemets_to_remove)
            %remove results, which havent used estimated autocorrelations matrix
            obj.cell_specific.biterrors_coded(1:elemets_to_remove,:,:) = [];
            obj.cell_specific.biterrors_uncoded(1:elemets_to_remove,:,:) = [];
            obj.cell_specific.FER_coded(1:elemets_to_remove,:,:) = [];
            obj.cell_specific.channel_error(1:elemets_to_remove,:,:,:) = [];
            obj.cell_specific.blocksize_coded(1:elemets_to_remove,:,:) = [];
            obj.cell_specific.blocksize_uncoded(1:elemets_to_remove,:,:) = [];
            obj.cell_specific.used_codewords(1:elemets_to_remove,:,:) = [];
            obj.cell_specific.throughput_coded(1:elemets_to_remove,:,:) = [];
            obj.cell_specific.throughput_uncoded(1:elemets_to_remove,:,:) = [];
            obj.cell_specific.FER_uncoded(1:elemets_to_remove,:,:) = [];
            
            % Cell specific results
            obj.cell_specific.BER_coded   = squeeze(sum(obj.cell_specific.biterrors_coded,1)   ./ sum(obj.cell_specific.blocksize_coded,1));
            obj.cell_specific.BER_uncoded = squeeze(sum(obj.cell_specific.biterrors_uncoded,1) ./ sum(obj.cell_specific.blocksize_uncoded,1));
            obj.cell_specific.BLER         = squeeze(sum(obj.cell_specific.FER_coded,1) ./ sum(obj.cell_specific.used_codewords,1));

            obj.cell_specific.BER_coded_overall   = squeeze(sum(sum(obj.cell_specific.biterrors_coded,1),3)   ./ sum(sum(obj.cell_specific.blocksize_coded,1),3));
            obj.cell_specific.BER_uncoded_overall = squeeze(sum(sum(obj.cell_specific.biterrors_uncoded,1),3) ./ sum(sum(obj.cell_specific.blocksize_uncoded,1),3));
            obj.cell_specific.BLER_overall = squeeze(sum(sum(obj.cell_specific.FER_coded,1),3) ./ sum(sum(obj.cell_specific.used_codewords,1),3));
            obj.cell_specific.MSE_overall = mean(mean(mean(obj.cell_specific.channel_error,4),3),1);
            
            % UE-specific results
            for uu = 1:obj.nUE
                obj.UE_specific(uu).biterrors_coded(1:elemets_to_remove,:,:) = [];
                obj.UE_specific(uu).biterrors_uncoded(1:elemets_to_remove,:,:) = [];
                obj.UE_specific(uu).FER_coded(1:elemets_to_remove,:,:) = [];
                obj.UE_specific(uu).blocksize_coded(1:elemets_to_remove,:,:) = [];
                obj.UE_specific(uu).blocksize_uncoded(1:elemets_to_remove,:,:) = [];
                obj.UE_specific(uu).used_codewords(1:elemets_to_remove,:,:) = [];
                
                obj.UE_specific(uu).BER_coded   = sum(obj.UE_specific(uu).biterrors_coded,1)   ./ sum(obj.UE_specific(uu).blocksize_coded,1);
                obj.UE_specific(uu).BER_uncoded = sum(obj.UE_specific(uu).biterrors_uncoded,1) ./ sum(obj.UE_specific(uu).blocksize_uncoded,1);
                obj.UE_specific(uu).BLER         = squeeze(sum(obj.UE_specific(uu).FER_coded,1) ./ sum(obj.UE_specific(uu).used_codewords,1));
                
                obj.UE_specific(uu).BER_coded_overall   = squeeze(sum(sum(obj.UE_specific(uu).biterrors_coded,1),3)   ./ sum(sum(obj.UE_specific(uu).blocksize_coded,1),3));
                obj.UE_specific(uu).BER_uncoded_overall = squeeze(sum(sum(obj.UE_specific(uu).biterrors_uncoded,1),3) ./ sum(sum(obj.UE_specific(uu).blocksize_uncoded,1),3));
                obj.UE_specific(uu).BLER_overall        = squeeze(sum(sum(obj.UE_specific(uu).FER_coded,1),3) ./ sum(sum(obj.UE_specific(uu).used_codewords,1),3));
            end
        end
    end
end
