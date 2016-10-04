clearvars; close all;

OUTPUT_FILENAME='sources.mtx';

%% Requiered Parameters
X=[50]; % Coordinates in X (Grid points)
Y=[50]; % Coordinates in Y (Grid points)
Z=[50]; % Coordinates in Z (Grid points)
SOURCE_TYPE=[1]; % Source Type (1=P-Source)
WAVELET_TYPE=[1]; % Wavelet Type (1=Synthetic)

%% Optional Parameters
WAVELET_SHAPE=[1]; % Wavelet Shape (1=Ricker)
FC=[5]; % Center Frequency in Hz
AMP=[5]; % Amplitude
TShift=[0]; % Time shift in s


%% Write to file

% Create Matrix
SOURCE_FILE=[X' Y' Z' SOURCE_TYPE' WAVELET_TYPE' WAVELET_SHAPE' FC' AMP' TShift'];

% Write mtx file
write_mtx(OUTPUT_FILENAME,SOURCE_FILE);