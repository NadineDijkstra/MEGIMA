function [frameHz, pixPerDeg, calibrationFile, screenDistance] = get_monitor_info(location)

% [frameHz, pixPerDeg, calibrationFile] = getMonitorInfo(location)
%
% Arguments:
%   location:   'debug'         = Sander's desktop
%               'behavioural'   = Behavioural lab 1
%               'trio'          = Trio scanner
%               'random'        = Random computer
% By SEB, Dec 2013

switch location
    case 'debug'
        % Sander's desktop computer (DCCN583)
        % Screen width = aa cm
        % Screen height = bb cm
        % Distance from eye to screen = 60 cm
        % Screen resolution = 1680 x 1050 pixels
        % vf_horizontal = atan(aa/2/cc)*180/pi*2 = dd deg, pixPerDeg = 1680/dd = x deg
        % vf_vertical  = atan(bb/2/cc)*180/pi*2 = ee deg, pixPerDeg = 1050/ee = y deg
        
        screenDistance = 60;
        pixPerDeg = 29;             % pixPerDeg = (x + y)/2
        frameHz = 60;               % frameHz = Screen('FrameRate', scrnum)
        calibrationFile = [];
        
    case 'behavioural'
        % Behavioural lab 1
        % Screen width: 38.5 cm
        % Screen height: 29 cm
        % Distance from eye to screen: 80 cm
        % Screen resolution = 1024 x 768 pixels
        % vf_horizontal = atan(38.5/2/80)*180/pi*2 = 27.0592 deg, pixPerDeg = 1024/31.0852 = 37.8430
        % vf_vertical  = atan(29/2/80)*180/pi*2 = 20.5467 deg, pixPerDeg = 768/23.5137 = 37.3783
        
        screenDistance = 80;
        pixPerDeg = 37.6;           % pixPerDeg = (x + y)/2
        frameHz = 60;               % frameHz = Screen('FrameRate', scrnum)
        calibrationFile = [];
        
    case 'mri'
        % Projected image, screen of Donders Trio
        % Screen width: 38.5 cm
        % Screen height: 29 cm
        % Distance from eye to screen: 80 cm
        % Screen resolution = 1024 x 768 pixels
        % visualField_horizontal = atan(38.5/2/80)*180/pi*2 = 27.0592 deg, pixPerDeg = 1024/31.0852 = 37.8430
        % visualField_vertical  = atan(29/2/80)*180/pi*2 = 20.5467 deg, pixPerDeg = 768/23.5137 = 37.3783

        screenDistance = 80;
        pixPerDeg = 40;           % pixPerDeg = (x + y)/2
        frameHz = 60;               % frameHz = Screen('FrameRate', scrnum)
        calibrationFile = [];
     
    case 'meg'
        % Projected image, screen of Donders MEG
        % Screen width: 38.5 cm
        % Screen height: 29 cm
        % Distance from eye to screen: 80 cm
        % Screen resolution = 1024 x 768 pixels
        % visualField_horizontal = atan(38.5/2/80)*180/pi*2 = 27.0592 deg, pixPerDeg = 1024/31.0852 = 37.8430
        % visualField_vertical  = atan(29/2/80)*180/pi*2 = 20.5467 deg, pixPerDeg = 768/23.5137 = 37.3783

        screenResolution = [1024 768];
        screenDistance = .8;
        screenWidth = .492;
        pixPerDeg = (tan(pi/180)*screenDistance*screenResolution(1))/screenWidth;
        frameHz = 60;               % frameHz = Screen('FrameRate', scrnum)
        calibrationFile = [];
       
%         screenWidth = 49.2;
%         distanceToScreen = 80;
        
        
    case 'random'
        % Any random computer
        screenWidth = input('What is the screen width in centimeters?');
        screenHeight = input('What is the screen height in centimeters?');
        screenDistance = input('What is the distance to the screen in centimeters?');
        screenNumber = input('What is the screen number you want to use? Default is 0');
        [screenResolution(1), screenResolution(2)] = Screen('WindowSize', screenNumber);
        visualField_horizontal = atan(screenWidth/2/distanceToScreen)*180/pi*2;    % in degrees
        horizontal_pixPerDeg = screenResolution(1)/visualField_horizontal;         % in degrees
        visualField_vertical = atan(screenHeight/2/distanceToScreen)*180/pi*2;     % in degrees
        vertical_pixPerDeg = screenResolution(2)/visualField_vertical;             % in degrees
        
        pixPerDeg = (horizontal_pixPerDeg + vertical_pixPerDeg)/2;
        frameHz = Screen('FrameRate', screenNumber);
        calibrationFile = [];
        
    otherwise
        error('Invalid location entered');
end






