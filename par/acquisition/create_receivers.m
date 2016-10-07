clearvars; close all;

OUTPUT_FILENAME='receiver.mtx';

%% Requiered Parameters
X=[70 75]; % Coordinates in X (Grid points)
Y=[70 75]; % Coordinates in Y (Grid points)
Z=[70 75]; % Coordinates in Z (Grid points)
RECEIVER_TYPE=[1 1]; % RECEIVER Type (1=P-RECEIVER)

%% Write to file

% Create Matrix
RECEIVER_FILE=[X' Y' Z' RECEIVER_TYPE'];

% Write mtx file
write_mtx(OUTPUT_FILENAME,RECEIVER_FILE);