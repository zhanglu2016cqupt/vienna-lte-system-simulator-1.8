classdef miscUtils
    % Implements miscellaneous functions related to the trace generation.
    % Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at.
    % (c) 2011 by INTHFT
    % www.nt.tuwien.ac.at
    
    properties
    end
    
    methods(Static)
        function LTE_params = LTE_params_function
            % Re-create needed load_parameters data from Link level for the generation of the precoding matrices.
            % (c) Josep Colom Ikuno, INTHFT, 2008/2011
            % www.nt.tuwien.ac.at
            
            
            %% Create the Codebook for Precoding
            
            % Transmit diversity
            LTE_params.Z{1} =  [1, 0, 1i,  0;
                0,-1,  0, 1i;
                0, 1,  0, 1i;
                1, 0,-1i,  0];
            LTE_params.Z{2} =  [1, 0, 0, 0, 1i,  0,  0, 0;
                0, 0, 0, 0,  0,  0,  0, 0;
                0,-1, 0, 0,  0, 1i,  0, 0;
                0, 0, 0, 0,  0,  0,  0, 0;
                0, 1, 0, 0,  0, 1i,  0, 0;
                0, 0, 0, 0,  0,  0,  0, 0;
                1, 0, 0, 0,-1i,  0,  0, 0;
                0, 0, 0, 0,  0,  0,  0, 0;
                0, 0, 0, 0,  0,  0,  0, 0;
                0, 0, 1, 0,  0,  0, 1i, 0;
                0, 0, 0, 0,  0,  0,  0, 0;
                0, 0, 0,-1,  0,  0,  0,1i;
                0, 0, 0, 0,  0,  0,  0, 0;
                0, 0, 0, 1,  0,  0,  0,1i;
                0, 0, 0, 0,  0,  0,  0, 0;
                0, 0, 1, 0,  0,  0,-1i, 0];
            
            % Spatial multiplexing
            U_temp = [  1,-1,-1,-1;     % Matrix corresponding to vectors u0 ... u15 in Table 6.3.4.2.3-2
                1,-1i,1,1i;
                1,1,-1,1;
                1,1i,1,-1i;
                1,(-1-1i)/sqrt(2), -1i,(1-1i)/sqrt(2);
                1,(1-1i)/sqrt(2), 1i,(-1-1i)/sqrt(2);
                1,(1+1i)/sqrt(2), -1i,(-1+1i)/sqrt(2);
                1,(-1+1i)/sqrt(2), 1i,(1+1i)/sqrt(2);
                1,-1,1,1;
                1,-1i,-1,-1i;
                1,1,1,-1;
                1,1i,-1,1i;
                1,-1,-1,1;
                1,-1,1,-1;
                1,1,-1,-1;
                1,1,1,1;].';
            Wn = zeros(4,4,16);
            for ii = 1:16
                LTE_params.Wn(:,:,ii)=diag(ones(1,4))-2*U_temp(:,ii)*U_temp(:,ii)'/(U_temp(:,ii)'*U_temp(:,ii));
            end
            
            % W Matrix according to Table 6.3.4.2.3-1
            %  LTE_params.W{1} = cat(3,[1;0],[0;1],[1/sqrt(2);1/sqrt(2)],[1/sqrt(2);-1/sqrt(2)],...
            %         [1/sqrt(2);1i/sqrt(2)],[1/sqrt(2);-1i/sqrt(2)]);
            LTE_params.W{1} = cat(3,[1/sqrt(2);1/sqrt(2)],[1/sqrt(2);-1/sqrt(2)],...
                [1/sqrt(2);1i/sqrt(2)],[1/sqrt(2);-1i/sqrt(2)]);
            LTE_params.W{2} = cat(3,1/sqrt(2)*[1,0;0,1],1/(2)*[1,1;1,-1],1/(2)*[1,1;1i,-1i]);
            
            % Large delay CDD
            LTE_params.U_l{1} = 1;
            LTE_params.U_l{2} = 1/sqrt(2)*[1,1;1,exp(-1i*pi)];
            LTE_params.U_l{3} = 1/sqrt(3)*[1,1,1;1,exp(-1i*2*pi/3),exp(-1i*4*pi/3);1,exp(-1i*4*pi/3),exp(-1i*8*pi/3)];
            LTE_params.U_l{4} = 1/2*[1,1,1,1;1,exp(-1i*2*pi/4),exp(-1i*4*pi/4),exp(-1i*6*pi/4);...
                1,exp(-1i*4*pi/4),exp(-1i*8*pi/4),exp(-1i*12*pi/4);...
                1,exp(-1i*6*pi/4),exp(-1i*12*pi/4),exp(-1i*18*pi/4)];
            LTE_params.D_l{1} = 1;
            LTE_params.D_l{2} = [1,0;0,exp(-1i*pi)];
            LTE_params.D_l{3} = [1,0,0;0,exp(-1i*2*pi/3),0;0,0,exp(-1i*4*pi/3)];
            LTE_params.D_l{4} = [1,0,0,0;0,exp(-1i*2*pi/4),0,0;0,0,exp(-1i*4*pi/4),0;0,0,0,exp(-1i*6*pi/4)];
            
            % Note that as of v.8.3.0, small delay CDD is removed from the standard
            % (28/05/08	RAN_40	RP-080432	0043	-	Removal of small-delay CDD
            
            % Precoding matrix W columns to take for each layer mapping
            LTE_params.mapping{1} = ones(16,1);
            LTE_params.mapping{2}=[1 4;1 2;1 2;1 2;1 4;1 4;1 3;1 3;1 2;1 4;1 3;1 3;1 2;1 3;1 3;1 2];
            LTE_params.mapping{3}=[1 2 4;1 2 3;1 2 3;1 2 3;1 2 4;1 2 4;1 3 4;1 3 4;1 2 4;1 3 4;1 2 3;1 3 4;1 2 3;1 2 3;1 2 3;1 2 3];
            LTE_params.mapping{4}=[1 2 3 4;1 2 3 4;3 2 1 4;3 2 1 4;1 2 3 4;1 2 3 4;1 3 2 4;
                1 3 2 4;1 2 3 4;1 2 3 4;1 3 2 4;1 3 2 4;1 2 3 4;1 3 2 4;3 2 1 4;1 2 3 4];
        end
        
        function LTE_params_rel10 = LTE_params_function_rel10
            
           %Specifie Codebook entries for TxMode 9 according TS 36.211 V10.0.0 
           %Table 6.3.4.2.3-3 to 6.3.4.2.3-10. Due to the user-specific reference
           %signals other Codebooks might be possible
           
           %cf. TS 36.211 V10.0.0: i1 = m, i2 = n
           
           
           %1 Layer:
           index = 0;
           for k = 0:3
               for n = 0:3
                   for m = 0:15 
                       
                        phi = exp(1i*pi*n/2);
                        v   = [1;exp(1i*2*pi*(2*m+k)/32);exp(1i*4*pi*(2*m+k)/32);exp(1i*6*pi*(2*m+k)/32)];
                        LTE_params_rel10{1}.W(:,:,index+1) = 1/sqrt(8) * [v;phi.*v];
                        LTE_params_rel10{1}.codebook_index(index+1) = index;
                        index = index + 1;
                                           
                   end               
               end
           end
            
            
            %2 Layers:
            index_1 = [0 1 2 3 0 1 0 1];
            index_2 = [0 1 2 3 1 2 3 3];
            index = 0;
            for k = 1:8
               for n = 0:1
                  for m = 0:15     
                                           
                        phi = exp(1i*pi*n/2);
                        v   = [1;exp(1i*2*pi*(2*m+index_1(k))/32);exp(1i*4*pi*(2*m+index_1(k))/32);exp(1i*6*pi*(2*m+index_1(k))/32)];
                        v_2 = [1;exp(1i*2*pi*(2*m+index_2(k))/32);exp(1i*4*pi*(2*m+index_2(k))/32);exp(1i*6*pi*(2*m+index_2(k))/32)];
                        LTE_params_rel10{2}.W(:,:,index+1) = 1/4 * [v v_2;phi.*v -phi*v_2];
                        LTE_params_rel10{2}.codebook_index(index+1) = index;
                        index = index + 1;
                                           
                   end               
               end
            end
            
            
            %3 Layers:
            index_1 = [0 8 2  10  4 12  6 14];
            index_2 = [8 0 10  2 12  4 14  6];
            index_3 = [8 0 10  2 12  4 14  6];
            index = 0;
            for k = 1:8
                  for m = 0:3     
                                           
                        v   = [1;exp(1i*2*pi*(8*m+index_1(k))/32);exp(1i*4*pi*(8*m+index_1(k))/32);exp(1i*6*pi*(8*m+index_1(k))/32)];
                        v_2 = [1;exp(1i*2*pi*(8*m+index_2(k))/32);exp(1i*4*pi*(8*m+index_2(k))/32);exp(1i*6*pi*(8*m+index_2(k))/32)];
                        v_3 = [1;exp(1i*2*pi*(8*m+index_3(k))/32);exp(1i*4*pi*(8*m+index_3(k))/32);exp(1i*6*pi*(8*m+index_3(k))/32)];
                        LTE_params_rel10{3}.W(:,:,index+1) = 1/sqrt(24) * [v v_2 v_3;v -v_2 -v_3];
                        LTE_params_rel10{3}.codebook_index(index+1) = index;
                        index = index + 1;
                                           
                   end               
            end
            index_1 = [0 8 2  10  4 12  6 14];
            index_2 = [0 0 2  2   4 4   6 6 ];
            index_3 = [8 8 10 10 12 12 14 14];
            for k = 1:8
                 for m = 0:3     
                                           
                        v   = [1;exp(1i*2*pi*(8*m+index_1(k))/32);exp(1i*4*pi*(8*m+index_1(k))/32);exp(1i*6*pi*(8*m+index_1(k))/32)];
                        v_2 = [1;exp(1i*2*pi*(8*m+index_2(k))/32);exp(1i*4*pi*(8*m+index_2(k))/32);exp(1i*6*pi*(8*m+index_2(k))/32)];
                        v_3 = [1;exp(1i*2*pi*(8*m+index_3(k))/32);exp(1i*4*pi*(8*m+index_3(k))/32);exp(1i*6*pi*(8*m+index_3(k))/32)];
                        LTE_params_rel10{3}.W(:,:,index+1) = 1/sqrt(24) * [v v_2 v_3;v v_2 -v_3];
                        LTE_params_rel10{3}.codebook_index(index+1) = index;
                        index = index + 1;
                                           
                   end               
            end
            
            
            %4 Layers:
            index_1 = [0  2  4   6];
            index_2 = [8 10 12  14];
            index = 0;
            for k = 1:4
                for n = 0:1
                   for m = 0:3     
                           
                        phi = exp(1i*pi*n/2);
                        v   = [1;exp(1i*2*pi*(8*m+index_1(k))/32);exp(1i*4*pi*(8*m+index_1(k))/32);exp(1i*6*pi*(8*m+index_1(k))/32)];
                        v_2 = [1;exp(1i*2*pi*(8*m+index_2(k))/32);exp(1i*4*pi*(8*m+index_2(k))/32);exp(1i*6*pi*(8*m+index_2(k))/32)];
                        LTE_params_rel10{4}.W(:,:,index+1) = 1/sqrt(32) * [v v_2 v v_2;phi*v phi*v_2 -phi*v -phi*v_2];
                        LTE_params_rel10{4}.codebook_index(index+1) = index;
                        index = index + 1;
                                           
                   end  
                end
            end
            
            
            %5 Layers:
            for m = 0:3     
                           
                v   = [1;exp(1i*2*pi*(2*m)/32);exp(1i*4*pi*(2*m)/32);exp(1i*6*pi*(2*m)/32)];
                v_2 = [1;exp(1i*2*pi*(2*m+8)/32);exp(1i*4*pi*(2*m+8)/32);exp(1i*6*pi*(2*m+8)/32)];
                v_3 = [1;exp(1i*2*pi*(2*m+16)/32);exp(1i*4*pi*(2*m+16)/32);exp(1i*6*pi*(2*m+16)/32)];
                LTE_params_rel10{5}.W(:,:,m+1) = 1/sqrt(40) * [v v v_2 v_2 v_3;v -v v_2 -v_2 v_3];
                LTE_params_rel10{5}.codebook_index(m+1) = m;
            end
            
            
            %6 Layers:
            for m = 0:3     
                           
                v   = [1;exp(1i*2*pi*(2*m)/32);exp(1i*4*pi*(2*m)/32);exp(1i*6*pi*(2*m)/32)];
                v_2 = [1;exp(1i*2*pi*(2*m+8)/32);exp(1i*4*pi*(2*m+8)/32);exp(1i*6*pi*(2*m+8)/32)];
                v_3 = [1;exp(1i*2*pi*(2*m+16)/32);exp(1i*4*pi*(2*m+16)/32);exp(1i*6*pi*(2*m+16)/32)];
                LTE_params_rel10{6}.W(:,:,m+1) = 1/sqrt(48) * [v v v_2 v_2 v_3 v_3;v -v v_2 -v_2 v_3 -v_3];
                LTE_params_rel10{6}.codebook_index(m+1) = m;
            end
            
           
           %7 Layers:
           for m = 0:3     
                           
                v   = [1;exp(1i*2*pi*(2*m)/32);exp(1i*4*pi*(2*m)/32);exp(1i*6*pi*(2*m)/32)];
                v_2 = [1;exp(1i*2*pi*(2*m+8)/32);exp(1i*4*pi*(2*m+8)/32);exp(1i*6*pi*(2*m+8)/32)];
                v_3 = [1;exp(1i*2*pi*(2*m+16)/32);exp(1i*4*pi*(2*m+16)/32);exp(1i*6*pi*(2*m+16)/32)];
                v_4 = [1;exp(1i*2*pi*(2*m+24)/32);exp(1i*4*pi*(2*m+24)/32);exp(1i*6*pi*(2*m+24)/32)];
                LTE_params_rel10{7}.W(:,:,m+1) = 1/sqrt(56) * [v v v_2 v_2 v_3 v_3 v_4;v -v v_2 -v_2 v_3 -v_3 v_4];
                LTE_params_rel10{7}.codebook_index(m+1) = m;
           end
           
           
           %8 Layers:
           v   = [1;exp(1i*2*pi*(2*m)/32);exp(1i*4*pi*(2*m)/32);exp(1i*6*pi*(2*m)/32)];
           v_2 = [1;exp(1i*2*pi*(2*m+8)/32);exp(1i*4*pi*(2*m+8)/32);exp(1i*6*pi*(2*m+8)/32)];
           v_3 = [1;exp(1i*2*pi*(2*m+16)/32);exp(1i*4*pi*(2*m+16)/32);exp(1i*6*pi*(2*m+16)/32)];
           v_4 = [1;exp(1i*2*pi*(2*m+24)/32);exp(1i*4*pi*(2*m+24)/32);exp(1i*6*pi*(2*m+24)/32)];
           LTE_params_rel10{8}.W(:,:,1) = 1/sqrt(56) * [v v v_2 v_2 v_3 v_3 v_4 v_4;v -v v_2 -v_2 v_3 -v_3 v_4 -v_4];
           LTE_params_rel10{8}.codebook_index(1) = 0;
                                            
            
        end
        
        function precoding_config = get_all_precoding_combinations
            % This small helper function returns all possible precoding options for LTE.
            % (c) Josep Colom Ikuno, INTHFT, 2008
            % www.nt.tuwien.ac.at
            
            precoding_config = cell(8,8,9); % Up to 4 layers, up to 4 TX antennas, 6 TX modes
            LTE_params       = phy_modeling.miscUtils.LTE_params_function;
            LTE_params_rel10 = phy_modeling.miscUtils.LTE_params_function_rel10;
            
            for tx_mode = 1:9
                switch tx_mode
                    case 1
                        % SISO
                        precoding_config{1,1,tx_mode}.tx_mode = tx_mode;
                        precoding_config{1,1,tx_mode}.nAtPort = 1;
                        precoding_config{1,1,tx_mode}.nLayers = 1;
                    case 2
                        % TxD
                        for nAtPort = [2 4]
                            precoding_config{nAtPort,nAtPort,tx_mode}.tx_mode = tx_mode;
                            precoding_config{nAtPort,nAtPort,tx_mode}.nAtPort = nAtPort;
                            precoding_config{nAtPort,nAtPort,tx_mode}.nLayers = nAtPort;
                            
                            %% Codebook setting
                            % We call the precoding matrix of TxD Z: Matrix corresponding to 36.211 section 6.3.4.3
                            precoding_config{nAtPort,nAtPort,tx_mode}.Z       = LTE_params.Z{nAtPort/2};
                            precoding_config{nAtPort,nAtPort,tx_mode}.name    = 'TxD';
                        end
                    case {3 4}
                        % OLSM, CLSM
                        for nAtPort = [2 4]
                            for nLayers = 1:nAtPort
                                precoding_config{nLayers,nAtPort,tx_mode}.tx_mode = tx_mode;
                                precoding_config{nLayers,nAtPort,tx_mode}.nAtPort = nAtPort;
                                precoding_config{nLayers,nAtPort,tx_mode}.nLayers = nLayers;
                                
                                %% Codebook setting
                                switch tx_mode
                                    case 4
                                        % CLSM
                                        switch nAtPort
                                            case 2
                                                switch nLayers
                                                    case 1
                                                        codebook_indexs = 0:3;
                                                    case 2
                                                        codebook_indexs = 1:2;
                                                end
                                            case 4
                                                codebook_indexs = 0:15;
                                        end
                                        
                                        % Closed loop spatial multiplexing, section 6.3.4.2.1
                                        W = zeros(nAtPort,nLayers,length(codebook_indexs));
                                        if (nAtPort == 2)
                                            if (min(codebook_indexs)<0 || max(codebook_indexs)>3) && nLayers ==1
                                                error('Only codebooks 0-3 are defined for %d layers (see TS.36.211, Table 6.3.4.2.3-1)',nLayers);
                                            elseif (min(codebook_indexs)<0 || max(codebook_indexs)>2) && nLayers ==2
                                                error('Only codebooks 0-2 are defined for %d layers (see TS.36.211, Table 6.3.4.2.3-1)',nLayers);
                                            end
                                            for cb_ = 1:length(codebook_indexs)
                                                codebook_index = codebook_indexs(cb_);
                                                W(:,:,cb_) = LTE_params.W{nLayers}(:,:,codebook_index+1);
                                            end
                                        else
                                            if min(codebook_indexs)<0 || max(codebook_indexs)>15
                                                error('Only codebooks 0-15 are defined (see TS.36.211, Table 6.3.4.2.3-2)');
                                            end
                                            for cb_ = 1:length(codebook_indexs)
                                                codebook_index = codebook_indexs(cb_);
                                                W_temp = 1/sqrt(nLayers)*LTE_params.Wn(:,:,codebook_index+1);
                                                W(:,:,cb_) = W_temp(:,LTE_params.mapping{nLayers}(codebook_index+1,:),1);
                                            end
                                        end
                                        
                                        % To avoid crashes with automatically generated code due to unexpected real values when complex ones are expected
                                        if isreal(W)
                                            W = complex(W);
                                        end
                                        
                                        precoding_config{nLayers,nAtPort,tx_mode}.W              = W;
                                        precoding_config{nLayers,nAtPort,tx_mode}.name           = 'CLSM';
                                        precoding_config{nLayers,nAtPort,tx_mode}.codebook_index = codebook_indexs;
                                        
                                    case 3
                                        % OLSM uses codebooks 12-15 in a cyclic way. Thus we set codebook_index to [12 13 14 15]
                                        switch nAtPort
                                            case 4
                                                codebook_index = [12 13 14 15]-1;
                                            case 2
                                                codebook_index = 0;
                                        end
                                        
                                        if (nAtPort == 2)
                                            W = LTE_params.W{nLayers}(:,:,codebook_index+1);
                                        else
                                            W_temp = 1/sqrt(nLayers)*LTE_params.Wn(:,:,codebook_index+1);
                                            % nLayers long Cyclic precoding matrix
                                            W = zeros(4,nLayers,4);
                                            for ii = 13:16
                                                W(:,:,ii-12) = W_temp(:,LTE_params.mapping{nLayers}(ii,:),ii-12);
                                            end
                                        end
                                        
                                        % To avoid crashes with automatically generated code due to unexpected real values when complex ones are expected
                                        if isreal(W)
                                            W = complex(W);
                                        end
                                        
                                        % Open loop spatial multiplexing, section 6.3.4.2.2 (Large CDD)
                                        precoding_config{nLayers,nAtPort,tx_mode}.W              = W;
                                        precoding_config{nLayers,nAtPort,tx_mode}.U              = LTE_params.U_l{nLayers};
                                        precoding_config{nLayers,nAtPort,tx_mode}.D              = LTE_params.D_l{nLayers};
                                        precoding_config{nLayers,nAtPort,tx_mode}.name           = 'OLSM';
                                        precoding_config{nLayers,nAtPort,tx_mode}.codebook_index = codebook_index;
                                end
                            end
                        end
                    case {5 6}
                        % MU-MIMO, Rank-1-SU-MIMO
                        for nAtPort = [2 4]
                            for nLayers = 1:1
                                precoding_config{nLayers,nAtPort,tx_mode}.tx_mode = tx_mode;
                                precoding_config{nLayers,nAtPort,tx_mode}.nAtPort = nAtPort;
                                precoding_config{nLayers,nAtPort,tx_mode}.nLayers = nLayers;
                                
                                %% Codebook setting
                                switch tx_mode
                                    case 6
                                        % Rank-1 CLSM
                                        switch nAtPort
                                            case 2
                                                switch nLayers
                                                    case 1
                                                        codebook_indexs = 0:3;
                                                    case 2
                                                        codebook_indexs = 1:2;
                                                end
                                            case 4
                                                codebook_indexs = 0:15;
                                        end
                                        
                                        % Closed loop spatial multiplexing, section 6.3.4.2.1
                                        W = zeros(nAtPort,nLayers,length(codebook_indexs));
                                        if (nAtPort == 2)
                                            if (min(codebook_indexs)<0 || max(codebook_indexs)>3) && nLayers ==1
                                                error('Only codebooks 0-3 are defined for %d layers (see TS.36.211, Table 6.3.4.2.3-1)',nLayers);
                                            elseif (min(codebook_indexs)<0 || max(codebook_indexs)>2) && nLayers ==2
                                                error('Only codebooks 0-2 are defined for %d layers (see TS.36.211, Table 6.3.4.2.3-1)',nLayers);
                                            end
                                            for cb_ = 1:length(codebook_indexs)
                                                codebook_index = codebook_indexs(cb_);
                                                W(:,:,cb_) = LTE_params.W{nLayers}(:,:,codebook_index+1);
                                            end
                                        else
                                            if min(codebook_indexs)<0 || max(codebook_indexs)>15
                                                error('Only codebooks 0-15 are defined (see TS.36.211, Table 6.3.4.2.3-2)');
                                            end
                                            for cb_ = 1:length(codebook_indexs)
                                                codebook_index = codebook_indexs(cb_);
                                                W_temp = 1/sqrt(nLayers)*LTE_params.Wn(:,:,codebook_index+1);
                                                W(:,:,cb_) = W_temp(:,LTE_params.mapping{nLayers}(codebook_index+1,:),1);
                                            end
                                        end
                                        
                                        % To avoid crashes with automatically generated code due to unexpected real values when complex ones are expected
                                        if isreal(W)
                                            W = complex(W);
                                        end
                                        
                                        precoding_config{nLayers,nAtPort,tx_mode}.W              = W;
                                        precoding_config{nLayers,nAtPort,tx_mode}.name           = 'Rank-1 CLSM';
                                        precoding_config{nLayers,nAtPort,tx_mode}.codebook_index = codebook_indexs;
                                        
                                    case 5
                                        % MU-MIMO
                                        switch nAtPort
                                            case 2
                                                switch nLayers
                                                    case 1
                                                        codebook_indexs = 0:3;
                                                    case 2
                                                        codebook_indexs = 1:2;
                                                end
                                            case 4
                                                codebook_indexs = 0:15;
                                        end
                                        
                                        % Closed loop spatial multiplexing, section 6.3.4.2.1
                                        W = zeros(nAtPort,nLayers,length(codebook_indexs));
                                        if (nAtPort == 2)
                                            if (min(codebook_indexs)<0 || max(codebook_indexs)>3) && nLayers ==1
                                                error('Only codebooks 0-3 are defined for %d layers (see TS.36.211, Table 6.3.4.2.3-1)',nLayers);
                                            elseif (min(codebook_indexs)<0 || max(codebook_indexs)>2) && nLayers ==2
                                                error('Only codebooks 0-2 are defined for %d layers (see TS.36.211, Table 6.3.4.2.3-1)',nLayers);
                                            end
                                            for cb_ = 1:length(codebook_indexs)
                                                codebook_index = codebook_indexs(cb_);
                                                W(:,:,cb_) = LTE_params.W{nLayers}(:,:,codebook_index+1);
                                            end
                                        else
                                            if min(codebook_indexs)<0 || max(codebook_indexs)>15
                                                error('Only codebooks 0-15 are defined (see TS.36.211, Table 6.3.4.2.3-2)');
                                            end
                                            for cb_ = 1:length(codebook_indexs)
                                                codebook_index = codebook_indexs(cb_);
                                                W_temp = 1/sqrt(nLayers)*LTE_params.Wn(:,:,codebook_index+1);
                                                W(:,:,cb_) = W_temp(:,LTE_params.mapping{nLayers}(codebook_index+1,:),1);
                                            end
                                        end
                                        
                                        % To avoid crashes with automatically generated code due to unexpected real values when complex ones are expected
                                        if isreal(W)
                                            W = complex(W);
                                        end
                                        
                                        precoding_config{nLayers,nAtPort,tx_mode}.W              = W;
                                        precoding_config{nLayers,nAtPort,tx_mode}.name           = 'MU-MIMO';
                                        precoding_config{nLayers,nAtPort,tx_mode}.codebook_index = codebook_indexs;
                                        
                                end
                            end
                        end
                        
                    case {7 8}
                        
                    case 9
                                        
                        
                        nAtPort = 8;
                        for nLayers = 1:nAtPort
                                precoding_config{nLayers,nAtPort,tx_mode}.tx_mode = tx_mode;
                                precoding_config{nLayers,nAtPort,tx_mode}.nAtPort = nAtPort;
                                precoding_config{nLayers,nAtPort,tx_mode}.nLayers = nLayers;
                                
                                precoding_config{nLayers,nAtPort,tx_mode}.W              = LTE_params_rel10{nLayers}.W;
                                precoding_config{nLayers,nAtPort,tx_mode}.name           = 'Up-to-8Layer-SU-MIMO';
                                precoding_config{nLayers,nAtPort,tx_mode}.codebook_index = LTE_params_rel10{nLayers}.codebook_index;
                        end
                end
            end
        end
    end
    
end