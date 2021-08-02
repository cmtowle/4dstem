function structDP = getDPgraphene(file)

% Clarissa Bhargava - 2021 April


% Input variables
atomsAbs = abs(readmatrix(file));


% Set size of beam intersection
spot = 40; % Angstroms

% Set cell dimension for each diffraction pattern to be taken
cellDimMini = [spot spot 1];


% Get the first section of the atoms based on x-coordinate
xBool = atomsAbs(:,1) < spot;
xFiltered = atomsAbs(xBool,:);

% Filter based on y-coordinate
yBool = xFiltered(:,2) < spot;
xyFiltered = xFiltered(yBool,:);


% Perform simulation
emdSTEM = PRISM01(xyFiltered,cellDimMini);
emdSTEM = PRISMmultislice(emdSTEM);


% Adjust diffraction pattern scaling
% structDP = fftshift(emdSTEM.MULTI4D);
% plotImage(structDP,'absolute',jetBlack, [0.00000 0.00001]);


