% MEXes and geenrates code for all of the C-source needed and/or optional for the simalator.
% (c) Josep Colom Ikuno, INTHFT, 2013

clear mex

fprintf('Compiling COST231_urban_micro_pathloss.c\n');
mex ./C-source/COST231_urban_micro_pathloss.c -output ./LTE_aux_COST231_urban_micro_pathloss
fprintf('\n');

% Note: this part will ONLY work if you have the MATLAB coder installed
try
    mex_config = coder.MexCodeConfig;
    mex_config.ExtrinsicCalls            = false;
    mex_config.SaturateOnIntegerOverflow = false;
    mex_config.IntegrityChecks           = false;
    mex_config.ResponsivenessChecks      = false;

    fprintf('C code generation of channel_gain_wrappers.shadowFadingMapClaussen.spatiallyCorrelateMap.m\n');
    cd ./+channel_gain_wrappers/+shadowFadingMapClaussen/
    codegen spatiallyCorrelateMap.m ...
        -config mex_config ...
        -args {double(1), coder.typeof(double(ones(1,1)),[inf 2], [1 0]), coder.typeof(double(ones(1,1)),[inf inf]), coder.typeof(double(ones(1,1)), [inf inf inf]), coder.typeof(double(ones(1,1)), [1 inf], [0 1])}...
        -o spatiallyCorrelateMap_mex ...
        -report
    cd ../..
    fprintf('\n');
    
    % Correctly define the struct containing the precoders
    interfering_precoders = coder.newtype('struct', struct('W',coder.typeof(complex(0), [inf inf])), [inf inf inf]);
    
    fprintf('C code generation of network_elements.UE.calculate_per_layer_interference_power.m\n');
    cd ./+network_elements/+UE/
    codegen calculate_per_layer_interference_power.m ...
        -config mex_config ...
        -args {coder.typeof(complex(double(ones(1,1,1,1))),[inf inf inf inf]), coder.typeof(complex(double(ones(1,1,1,1,1))),[inf inf inf inf inf]), interfering_precoders}...
        -o calculate_per_layer_interference_power_mex ...
        -report
    cd ../..
    fprintf('\n');
    
    fprintf('C code generation of network_elements.UE.calculate_effective_channel_and_receiver.m\n');
    cd ./+network_elements/+UE/
    codegen calculate_effective_channel_and_receiver.m ...
        -config mex_config ...
        -args {coder.typeof(complex(double(ones(1,1,1))),[inf inf inf]), coder.typeof(complex(double(ones(1,1,1))),[inf inf inf])}...
        -o calculate_effective_channel_and_receiver_mex ...
        -report
    cd ../..
    fprintf('\n');
    
    fprintf('C code generation of feedback_calculation.codebook_capacity.non_sparse.m\n');
    cd ./+feedback_calculation/+codebook_capacity/
    codegen non_sparse.m ...
        -config mex_config ...
        -args {coder.typeof(complex(double(ones(1,1,1))),[inf inf inf]), coder.typeof(complex(double(ones(1,1,1))),[inf inf inf]), coder.typeof('a string',[1 inf])}...
        -o non_sparse_mex ...
        -report
    cd ../..
    fprintf('\n');
end