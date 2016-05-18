function BICM_capacity_tables = LTE_common_calculate_BICM_capacity_tables
% Calculate the lookup tables from which the BICM capacity for each
% constellation used in LTE will be searched. This data will be used by the
% MIESM to average the several subcarrier SINRs into an effective SINR.
% (c) Josep Colom Ikuno, INTHFT, 2010
% www.nt.tuwien.ac.at

m_j = [2 4 6];
SNR = -20:1:30;
N_realizations = 500;

BICM_capacity_tables = struct('m_j',[],'I',[],'SNR',[]);
BICM_capacity_tables = repmat(BICM_capacity_tables,[1 length(m_j)]);

%matlabpool open;
for m_j_idx=1:length(m_j)
    BICM_capacity_tables(m_j_idx) = BICM_capacity(m_j(m_j_idx),SNR,N_realizations);
end
%matlabpool close;

figure;
BICM_capacity_tables_all = zeros(length(m_j),length(BICM_capacity_tables(m_j(1)).I));

for m_j_idx=1:length(m_j)
    BICM_capacity_tables_all(m_j_idx,:) = BICM_capacity_tables(m_j_idx).I;
    displaynames{m_j_idx} = sprintf('BICM capacity, %d-QAM',2^m_j(m_j_idx));
end

plot(BICM_capacity_tables(m_j(1)).SNR,BICM_capacity_tables_all);
legend(displaynames,'Location','Best');
xlabel('SNR [dB]');
ylabel('BICM Capacity I_{m_j}(\gamma)');
title('BICM capacity');
grid on;

output_file = sprintf('BICM_capacity_tables_%d_realizations',N_realizations);
savefile = fullfile('./data_files/',output_file);
save(savefile,'BICM_capacity_tables');

function BICM_capacity_table = BICM_capacity(m_j,SNR,N_realizations)
% Calculate the BICM capacity (I_{m_j}(SNR)) for the LTE constellation
% with m_j bits per symbol and averaging over N_realizations

[symbol_alphabets bittables] = get_constellations;

if length(symbol_alphabets)<m_j || isempty(symbol_alphabets{m_j})
    error('Constellation with m_j of %d not found',m_j);
else
    symbol_alphabet = symbol_alphabets{m_j};
    bittable        = bittables{m_j};
end

I = zeros(1,length(SNR));
SNR_lin_vect = 10.^(SNR/10);

for SNR_idx = 1:length(SNR)
    SNR_lin = SNR_lin_vect(SNR_idx);
    E_y = zeros(1,N_realizations);
    Y   = (randn(1,N_realizations)+1i*randn(1,N_realizations))/sqrt(2);
    for i_=1:m_j
        for b_=[0 1]
            gamma_bi_idx   = (bittable(i_,:)==b_);
            gamma          = symbol_alphabet;
            gamma_bi       = symbol_alphabet(gamma_bi_idx);
            Y_mat_gamma    = repmat(Y,[length(gamma) 1]);
            Y_mat_gamma_bi = repmat(Y,[length(gamma_bi) 1]);
            gamma_mat      = repmat(gamma,[1 N_realizations]);
            gamma_bi_mat   = repmat(gamma_bi,[1 N_realizations]);
            for z=gamma_bi.'
                numerator   = sum(exp(-abs(Y_mat_gamma-sqrt(SNR_lin)*(gamma_mat-z)).^2),1);
                denominator = sum(exp(-abs(Y_mat_gamma_bi-sqrt(SNR_lin)*(gamma_bi_mat-z)).^2),1);
                E_y = E_y + log2(numerator./denominator);
            end
        end
    end
    I(SNR_idx) = m_j - mean(E_y)/length(symbol_alphabet);
end

BICM_capacity_table.m_j = m_j;
BICM_capacity_table.I   = I;
BICM_capacity_table.SNR = SNR;

function [SymbolAlphabet bittable] = get_constellations(varargin)
% Symbol alphabet and bittables code taken from the configuration file of the
% LTE Link Level simulator (see the 3GPP 36 series, TS 36.211)
%  www.nt.tuwien.ac.at/ltesimulator

if isempty(varargin)
    plot_constellations = false;
else
    plot_constellations = varargin{1};
end

SymbolAlphabet{1} = [ 1+1j, -1-1j].'/sqrt(2);
SymbolAlphabet{2} = [ 1+1j, 1-1j, -1+1j, -1-1j].'/sqrt(2);
SymbolAlphabet{4} = [
    complex( 1,  1)
    complex( 1,  3)
    complex( 3,  1)
    complex( 3,  3)
    complex( 1, -1)
    complex( 1, -3)
    complex( 3, -1)
    complex( 3, -3)
    complex(-1,  1)
    complex(-1,  3)
    complex(-3,  1)
    complex(-3,  3)
    complex(-1, -1)
    complex(-1, -3)
    complex(-3, -1)
    complex(-3, -3)
    ] / sqrt(10);
SymbolAlphabet{6} = [
    complex( 3,  3)
    complex( 3,  1)
    complex( 1,  3)
    complex( 1,  1)
    complex( 3,  5)
    complex( 3,  7)
    complex( 1,  5)
    complex( 1,  7)
    complex( 5,  3)
    complex( 5,  1)
    complex( 7,  3)
    complex( 7,  1)
    complex( 5,  5)
    complex( 5,  7)
    complex( 7,  5)
    complex( 7,  7) % symbol 0-15
    complex( 3, -3)
    complex( 3, -1)
    complex( 1, -3)
    complex( 1, -1)
    complex( 3, -5)
    complex( 3, -7)
    complex( 1, -5)
    complex( 1, -7)
    complex( 5, -3)
    complex( 5, -1)
    complex( 7, -3)
    complex( 7, -1)
    complex( 5, -5)
    complex( 5, -7)
    complex( 7, -5)
    complex( 7, -7) % symbol 16-31
    complex(-3,  3)
    complex(-3,  1)
    complex(-1,  3)
    complex(-1,  1)
    complex(-3,  5)
    complex(-3,  7)
    complex(-1,  5)
    complex(-1,  7)
    complex(-5,  3)
    complex(-5,  1)
    complex(-7,  3)
    complex(-7,  1)
    complex(-5,  5)
    complex(-5,  7)
    complex(-7,  5)
    complex(-7,  7) % symbol 32-47
    complex(-3, -3)
    complex(-3, -1)
    complex(-1, -3)
    complex(-1, -1)
    complex(-3, -5)
    complex(-3, -7)
    complex(-1, -5)
    complex(-1, -7)
    complex(-5, -3)
    complex(-5, -1)
    complex(-7, -3)
    complex(-7, -1)
    complex(-5, -5)
    complex(-5, -7)
    complex(-7, -5)
    complex(-7, -7) ] / sqrt(42); % symbol 48-63

bittable{1} = logical([0,1]);
bittable{2} = logical([...
    0,1,0,1;
    0,0,1,1 ...
    ]);
bittable{4} = logical([...
    0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1;
    0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
    0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
    0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1 ...
    ]);
bittable{6} = logical([...
    0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1;
    0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;
    0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1;
    0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 ...
    ]);

if plot_constellations
    for i_=1:length(SymbolAlphabet)
        if ~isempty(SymbolAlphabet{i_})
            figure;
            plot(SymbolAlphabet{i_},'.k');
            hold on;
            title(sprintf('%d-QAM constellation',2^i_));
            for sym_=1:length(SymbolAlphabet{i_})
                text(real(SymbolAlphabet{i_}(sym_)),imag(SymbolAlphabet{i_}(sym_)),reshape(num2str(bittable{i_}(:,sym_),'%3.0f'),1,[]));
            end
        end
    end
end