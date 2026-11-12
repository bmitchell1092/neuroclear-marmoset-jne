% --- Function to create a custom red-white-blue colormap ---
function cmap = customRedBlue(m)
    if nargin < 1
        m = size(get(gcf, 'colormap'), 1);
    end
    
    % Define the custom blue and red colors
    blue = [0 0.4470 0.7410];   % Custom blue
    red = [0.9500 0.2250 0.0980]; % Custom red
    white = [1 1 1];  % White for the middle
    
    % Interpolate the colormap
    cmap = [linspace(blue(1), white(1), m/2)', linspace(blue(2), white(2), m/2)', linspace(blue(3), white(3), m/2)'; ...  % Blue to white
            linspace(white(1), red(1), m/2)', linspace(white(2), red(2), m/2)', linspace(white(3), red(3), m/2)'];     % White to red
end