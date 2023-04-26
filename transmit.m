clc;
close all;

%% Number of bits
N=2^7;

%% Text to binary
A = imread('images.jpeg');
%imshow(A);
B = imresize(A,[14 NaN]);
shape = size(B);
temp_x = shape(2);
imshow(B);
B_x = B(:,:,1);
B_y = B(:,:,2);
B_z = B(:,:,3);
array = [];
for i = 1:temp_x
    array = [array; B_x(:,i)];
end
for i = 1:temp_x
    array = [array; B_y(:,i)];
end
for i = 1:temp_x
    array = [array; B_z(:,i)];
end

b = dec2bin(array,8);
b_t=transpose(b);
txt_to_bin=b_t(:)-'0';

length_to_be = length(txt_to_bin);

goldseq = comm.GoldSequence('FirstPolynomial','x^5+x^2+1', ...
    'SecondPolynomial','x^5+x^4+x^3+x^2+1', ...
    'FirstInitialConditions',[0 0 0 0 1], ...
    'SecondInitialConditions',[0 0 0 0 1], ...
    'Index',4,'SamplesPerFrame',length_to_be);
encrypt_var = goldseq();

txt_to_bin_encrypted = xor(txt_to_bin, encrypt_var);


%% Adding syncronization bits
h_pn = commsrc.pn('GenPoly', [7 6 0], 'InitialStates', [0 0 0 0 0 0 1], 'NumBitsOut', N);
syncronization_bits = generate(h_pn);
messageLength = dec2bin(length(txt_to_bin), 16);
messageLengthTransposed = transpose(messageLength);
text_to_bin_messageLength = messageLengthTransposed(:) - '0';

%% [sync bits, message length, data bits]\
tx_text = [syncronization_bits; text_to_bin_messageLength; txt_to_bin_encrypted];

%% Differential Encoding
diff_bits = zeros(length(tx_text), 1);
diff_bits(1) = xor(0, tx_text(1));
for i=2:length(diff_bits)
    diff_bits(i) = xor(diff_bits(i-1), tx_text(i)); % XOR ing privious bit to current bit
end

%% Bit Sequence to Signal (Upsampling)
bits_I = 2*(diff_bits - 0.5);
sig_ipt_I = upsample(bits_I, 8);
FIR_coeff75 = Hps75.Numerator;
sig_RRC75_I = conv(FIR_coeff75,sig_ipt_I);
sig_RRC75_I = sig_RRC75_I/max(sig_RRC75_I);

csvwrite('Non_Coherent_Cosine_075_text_ref_PRC.csv', tx_text);
csvwrite('Non_Coherent_Cosine_075_text_data_PRC.csv', sig_RRC75_I);

Trigger_interval=(200*10^3*8)/length(sig_RRC75_I)
Burst_frequency=1/Trigger_interval
 
