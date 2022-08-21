function SI_UNIWARD(precover, coverPath, stegoPath, payload, Q)

if (Q ~= 75) && (Q ~= 85) && (Q ~=95), error('Input parameter "Q" can be only 75, 85 or 95.'); end;

PC_SPATIAL = double(precover);
temp = load(strcat('default_gray_jpeg_obj_', num2str(Q), '.mat'));
default_gray_jpeg_obj = temp.default_gray_jpeg_obj;
C_STRUCT = default_gray_jpeg_obj;
C_QUANT = C_STRUCT.quant_tables{1};

fun=@dct2;
xi= blkproc(double(PC_SPATIAL)-128,[8 8],fun);
% Quantization
fun = @(x) x./C_QUANT;
DCT_real = blkproc(xi,[8 8],fun);
DCT_rounded = round(DCT_real);

C_STRUCT.coef_arrays{1} = DCT_rounded;

e = DCT_rounded - DCT_real;             % Compute rounding error
sgn_e = sign(e);
sgn_e(e==0) = round(rand(sum(e(:)==0), 1)) * 2 - 1;
change = - sgn_e;

%%-------------------------------
%%      BEGIN costs
%%-------------------------------

wetConst = 10^13;
sgm = 2^(-6);

%% Get 2D wavelet filters - Daubechies 8
% 1D high pass decomposition filter
hpdf = [-0.0544158422, 0.3128715909, -0.6756307363, 0.5853546837, 0.0158291053, -0.2840155430, -0.0004724846, 0.1287474266, 0.0173693010, -0.0440882539, ...
        -0.0139810279, 0.0087460940, 0.0048703530, -0.0003917404, -0.0006754494, -0.0001174768];
% 1D low pass decomposition filter
lpdf = (-1).^(0:numel(hpdf)-1).*fliplr(hpdf);

F{1} = lpdf'*hpdf;
F{2} = hpdf'*lpdf;
F{3} = hpdf'*hpdf;

%% Pre-compute impact in spatial domain when a jpeg coefficient is changed
spatialImpact = cell(8, 8);
for bcoord_i=1:8
    for bcoord_j=1:8
        testCoeffs = zeros(8, 8);
        testCoeffs(bcoord_i, bcoord_j) = 1;
        spatialImpact{bcoord_i, bcoord_j} = idct2(testCoeffs)*C_QUANT(bcoord_i, bcoord_j);
    end
end

%% Pre compute impact on wavelet coefficients when a jpeg coefficient is changed
waveletImpact = cell(numel(F), 8, 8);
for Findex = 1:numel(F)
    for bcoord_i=1:8
        for bcoord_j=1:8
            waveletImpact{Findex, bcoord_i, bcoord_j} = imfilter(spatialImpact{bcoord_i, bcoord_j}, F{Findex}, 'full');
        end
    end
end

%% Create reference cover wavelet coefficients (LH, HL, HH)

% precover
padSize = max([size(F{1})'; size(F{2})']);
PC_SPATIAL_PADDED = padarray(PC_SPATIAL, [padSize padSize], 'symmetric'); % pad image

R_PC = cell(size(F));
for i=1:numel(F)
    R_PC{i} = imfilter(PC_SPATIAL_PADDED, F{i});
end

[k, l] = size(PC_SPATIAL);

rho = zeros(k, l);
C_xi = cell(3, 1);
S_xi = cell(3, 1);
%% Computation of costs
for row = 1:k
    for col = 1:l
        sub_e = e(row, col);
        modRow = mod(row-1, 8)+1;
        modCol = mod(col-1, 8)+1;        
        
        subRows = row-modRow-6+padSize:row-modRow+16+padSize;
        subCols = col-modCol-6+padSize:col-modCol+16+padSize;
     
        for fIndex = 1:3
            % compute residual
            R_PC_sub = R_PC{fIndex}(subRows, subCols);
            % get differences between precover and cover, stego
            wavCoverStegoDiff = waveletImpact{fIndex, modRow, modCol};
            % compute suitability
            
            C_xi{fIndex} = abs(wavCoverStegoDiff.*sub_e)./ (abs(R_PC_sub)+sgm);
            S_xi{fIndex} = abs(wavCoverStegoDiff.*(sub_e-sgn_e(row, col))) ./ (abs(R_PC_sub)+sgm);
        end
        C_rho = C_xi{1} + C_xi{2} + C_xi{3};
        S_rho = S_xi{1} + S_xi{2} + S_xi{3};
        rho(row, col) = sum(S_rho(:)) - sum(C_rho(:));
    end
end

rho = rho + 10^(-4);
rho(rho > wetConst) = wetConst;
rho(isnan(rho)) = wetConst;    
rho((DCT_rounded > 1022)  & (e > 0)) = wetConst;
rho((DCT_rounded < -1022) & (e < 0)) = wetConst;

% Avoid 04 coefficients with e=0.5
maxCostMat = false(size(rho));
maxCostMat(1:8:end, 1:8:end) = true;
maxCostMat(5:8:end, 1:8:end) = true;
maxCostMat(1:8:end, 5:8:end) = true;
maxCostMat(5:8:end, 5:8:end) = true;
rho(maxCostMat & (abs(e)>0.4999)) = wetConst;

save('rho_m.mat', 'rho');

%%-------------------------------
%%      END costs
%%-------------------------------

%% Compute message lenght for each run
nzAC = nnz(DCT_rounded)-nnz(DCT_rounded(1:8:end,1:8:end)); % number of nonzero AC DCT coefficients
totalMessageLength = round(payload*nzAC);

%% Embedding
% permutes path 
perm = randperm(numel(DCT_rounded));
    
[LSBs] = EmbeddingSimulator(DCT_rounded(perm), rho(perm)', totalMessageLength);

LSBs(perm) = LSBs;                           % inverse permutation
LSBs = reshape(LSBs, size(DCT_rounded));  % reshape LSB into image-sized matrix

% Create stego coefficients
temp = mod(DCT_rounded, 2);
S_COEFFS = zeros(size(DCT_rounded));
S_COEFFS(temp == LSBs) = DCT_rounded(temp == LSBs);
S_COEFFS(temp ~= LSBs) = DCT_rounded(temp ~= LSBs) + change(temp ~= LSBs);

S_STRUCT = C_STRUCT;
S_STRUCT.coef_arrays{1} = S_COEFFS;

% Save cover
jpeg_write(C_STRUCT, coverPath);

% Save stego
jpeg_write(S_STRUCT, stegoPath);

function [LSBs] = EmbeddingSimulator(x, rho, m)
       
    rho = rho';
    x = double(x);
    n = numel(x);
    
    lambda = calc_lambda(rho, m, n);
    pChange = 1 - (double(1)./(1+exp(-lambda.*rho)));
    
    randChange = rand(size(x));
    flippedPixels = (randChange < pChange);
    LSBs = mod(x + flippedPixels, 2);
    
    function lambda = calc_lambda(rho, message_length, n)

        l3 = 1e+3;
        m3 = double(message_length + 1);
        iterations = 0;
        while m3 > message_length
            l3 = l3 * 2;
            p = double(1)./(1 + exp(-l3 .* rho));
            m3 = binary_entropyf(p);
            iterations = iterations + 1;
            if (iterations > 10)
                lambda = l3;
                return;
            end
        end        
        
        l1 = 0; 
        m1 = double(n);        
        lambda = 0;
        
        alpha = double(message_length)/n;
        % limit search to 30 iterations
        % and require that relative payload embedded is roughly within 1/1000 of the required relative payload        
        while  (double(m1-m3)/n > alpha/1000.0 ) && (iterations<30)
            lambda = l1+(l3-l1)/2; 
            p = double(1)./(1+exp(-lambda.*rho));
            m2 = binary_entropyf(p);
    		if m2 < message_length
                l3 = lambda;
    			m3 = m2;
            else
    			l1 = lambda;
    			m1 = m2;
            end
    		iterations = iterations + 1;
        end
    end
    
    function Hb = binary_entropyf(p)
        p = p(:);
        Hb = (-p.*log2(p))-((1-p).*log2(1-p));
        Hb(isnan(Hb)) = 0;
        Hb = sum(Hb);
    end

end

end
