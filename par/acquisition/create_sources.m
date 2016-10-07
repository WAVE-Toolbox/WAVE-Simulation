clearvars; close all;

OUTPUT_FILENAME='sources.mtx';

%% Requiered Parameters
X=[50 60]; % Coordinates in X (Grid points)
Y=[50 60]; % Coordinates in Y (Grid points)
Z=[50 60]; % Coordinates in Z (Grid points)
SOURCE_TYPE=[1 1]; % Source Type (1=P-Source)
WAVELET_TYPE=[1 1]; % Wavelet Type (1=Synthetic)

%% Optional Parameters
WAVELET_SHAPE=[1 2]; % Wavelet Shape (1=Ricker, 2=FGaussian)
FC=[5 4]; % Center Frequency in Hz
AMP=[5 2]; % Amplitude
TShift=[0 0]; % Time shift in s


%% Write to file

% Create Matrix
SOURCE_FILE=[X' Y' Z' SOURCE_TYPE' WAVELET_TYPE' WAVELET_SHAPE' FC' AMP' TShift'];

% Write mtx file
write_mtx(OUTPUT_FILENAME,SOURCE_FILE);