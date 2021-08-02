function [array4dstem] = getDPgrapheneArray(file,stepsize)

% Clarissa Bhargava - 2021 April


% Input variables
atomsAbs = abs(readmatrix(file));
step = stepsize;


% Set size of beam intersection
spot = 80.0; % Angstroms
rad = spot/2;

% Set cell dimension for each diffraction pattern to be taken
cellDimMini = [spot spot 1];

% Get center points of subdivided simulation areas

centers = [];
array4dstem = [];

for xCoord = rad:step:ceil(max(atomsAbs(:,1)))-rad
    
    for yCoord = rad:step:ceil(max(atomsAbs(:,2)))-rad
    
    centers = [centers ; [xCoord yCoord]];
    
    % Get the first section of the atoms based on x-coordinate
    xBoolmin = atomsAbs(:,1) > xCoord - rad;
    xFilteredmin = atomsAbs(xBoolmin,:);
    xBoolmax = xFilteredmin(:,1) < xCoord + rad;
    xFilteredminmax = xFilteredmin(xBoolmax,:);

    % Filter based on y-coordinate
    yBoolmin = xFilteredminmax(:,2) > yCoord - rad;
    yFilteredmin = xFilteredminmax(yBoolmin,:);
    yBoolmax = yFilteredmin(:,2) < yCoord + rad;
    xyFiltered = yFilteredmin(yBoolmax,:);


    % Perform simulation
    emdSTEM = PRISM01(xyFiltered,cellDimMini);
    emdSTEM = PRISMmultislice(emdSTEM);
    
    % Save to array
    structDP = fftshift(emdSTEM.MULTI4D);
    array4dstem = cat(3, array4dstem, structDP);
    
    
    % Name file and save image
%     fileName = strcat('atomsBig3_',num2str(xCoord),'x',num2str(yCoord),'y');
%     plotImage(structDP,'absolute',jetBlack,[0.00000 0.000008],fileName);
    % disp(fileName)
    
    end
    
end

% disp(centers)

end

% Adjust diffraction pattern scaling
% structDP = fftshift(emdSTEM.MULTI4D);
% plotImage(structDP,'absolute',jetBlack, [0.00000 0.00001]);

