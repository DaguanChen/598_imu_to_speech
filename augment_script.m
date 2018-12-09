rng(0);

% clips = load('amod.mat');
 %clips = load('daniel.mat');
clips = load('zhijian.mat');

%test
%clip = clips.a_1_1;
% fs = 44100;
% soundsc(clip,fs);
% pause(1);
% soundsc(add_noise(clip,0.001),fs);
% pause(1);
% soundsc(dynamic_range_compression(clip,0.03,0.5),fs);
% pause(1);
% soundsc(change_speed(clip, 2),fs);
% pause(1);
% soundsc(shift_pitch(clip, 0.8),fs);

fields = fieldnames(clips);
for i = 1:numel(fields)
  noised.(fields{i}) = add_noise(clips.(fields{i}),0.0005);
  fast.(fields{i}) = change_speed(clips.(fields{i}), 0.8);
  slow.(fields{i}) = change_speed(clips.(fields{i}), 1.3);
  high_pitch.(fields{i}) = shift_pitch(clips.(fields{i}), 1.1);
  low_pitch.(fields{i}) = shift_pitch(clips.(fields{i}), 0.9);
end

% save('amod_noised.mat', '-struct', 'noised');
% save('amod_fast.mat', '-struct', 'fast');
% save('amod_slow.mat', '-struct', 'slow');
% save('amod_high_pitch.mat', '-struct', 'high_pitch');
% save('amod_low_pitch.mat', '-struct', 'low_pitch');

% save('daniel_noised.mat', '-struct', 'noised');
% save('daniel_fast.mat', '-struct', 'fast');
% save('daniel_slow.mat', '-struct', 'slow');
% save('daniel_high_pitch.mat', '-struct', 'high_pitch');
% save('daniel_low_pitch.mat', '-struct', 'low_pitch');
% % 
save('zhijian_noised.mat', '-struct', 'noised');
save('zhijian_fast.mat', '-struct', 'fast');
save('zhijian_slow.mat', '-struct', 'slow');
save('zhijian_high_pitch.mat', '-struct', 'high_pitch');
save('zhijian_low_pitch.mat', '-struct', 'low_pitch');

%Augment functions
function y_noisy = add_noise(y, noise_amplitude)
    wgn = noise_amplitude*randn(length(y),1);
    y_noisy = y + wgn;
    for i = 1:length(y)
        y_noisy(i) = min(max(y_noisy(i),-1), 1);
    end
end

function y_compressed = dynamic_range_compression(y, amp_threshold, compress_ratio)
    y_compressed = y;
    for i = 1:length(y_compressed)
        if abs(y_compressed(i)) > amp_threshold
            y_compressed(i) = sign(y_compressed(i))*(amp_threshold + ...
                (abs(y_compressed(i)) - amp_threshold)*compress_ratio);
        end
    end
end

function y_spedup = change_speed(y_orig, ratio)
    WindowLen = 256;
    AnalysisLen = 64;
    SynthesisLen = round(AnalysisLen*ratio);
    Hopratio = SynthesisLen/AnalysisLen;
    win = dsp.Window('Hanning', 'Sampling', 'Periodic');
    
    yprevwin = zeros(WindowLen-SynthesisLen,1);
    gain = 1/(WindowLen*sum(hanning(WindowLen,'periodic').^2)/SynthesisLen);
    unwrapdata = 2*pi*AnalysisLen*(0:WindowLen-1)'/WindowLen;
    yangle = zeros(WindowLen,1);
    firsttime = true;

    y_spedup = [];
    frame = 0;
    while frame*AnalysisLen+WindowLen < length(y_orig)
        y = y_orig(frame*AnalysisLen+1:frame*AnalysisLen+WindowLen);
        yfft = fft(y);

        % Convert complex FFT data to magnitude and phase.
        ymag       = abs(yfft);
        yprevangle = yangle;
        yangle     = angle(yfft);

        % Synthesis Phase Calculation
        % The synthesis phase is calculated by computing the phase increments
        % between successive frequency transforms, unwrapping them, and scaling
        % them by the ratio between the analysis and synthesis hop sizes.
        yunwrap = (yangle - yprevangle) - unwrapdata;
        yunwrap = yunwrap - round(yunwrap/(2*pi))*2*pi;
        yunwrap = (yunwrap + unwrapdata) * Hopratio;
        if firsttime
            ysangle = yangle;
            firsttime = false;
        else
            ysangle = ysangle + yunwrap;
        end

        % Convert magnitude and phase to complex numbers.
        ys = ymag .* complex(cos(ysangle), sin(ysangle));

        % IST-FFT
        release(win);
        ywin  = win(ifft(ys));    % Windowed IFFT

        % Overlap-add operation
        olapadd  = [ywin(1:end-SynthesisLen,:) + yprevwin; ...
                    ywin(end-SynthesisLen+1:end,:)];
        yistfft  = olapadd(1:SynthesisLen,:);
        yprevwin = olapadd(SynthesisLen+1:end,:);

        % Compensate for the scaling that was introduced by the overlap-add
        % operation
        yistfft = yistfft * gain;
        y_spedup = [y_spedup; real(yistfft)];
        frame = frame+1;
    end
end

function y_pitch_shift = shift_pitch(y, ratio)
    y_spedup = change_speed(y, ratio);
    y_pitch_shift = resample(y_spedup,44100,round(44100*ratio));
end
