% SetUpCL15-7v Ultrasound Tensiometry Code
% 
% This code is originally based on the  RunSetUpL11_5vFlashHFR_acquire.m - 
% "High Frame Rate" acquisition for ultrafast imaging code with a Vantage 
% 64LE system and was written for Vantage 4.4.0. 
% It has been modified by Darryl Thelen and Lauren Welte. 
% Please submit a report on GitHub if there are any issues.
% 
% It is intended to be used with the Phillips CL15-7V transducer.
% 
% The code is designed to function with several other scripts to manage the
% processing and visualization of ultrasound-tensiometry code. 
% 
% The Event structure controls the sequence of events. There are two
% "groups" of events. Initially, there is plane wave imaging to allow for
% placement of the ultrasound probe on the tissue of interest. The system
% will capture "reducedFrames" number of frames in a loop, until the system
% is told to enter the tensiometry mode. Once the "start" button is
% pressed, it will enter the tensiometry mode, which will send a trig out
% signal (to initiate the tap with an external set-up), and then will
% capture "tapSamples" number of plane wave imaging frames after each tap, 
% at the "sampleRate" selected (recommended 20,000Hz). It will repeat the
% taps for "numTaps", with a controlled tap frequency of "tapRate".
% 
% Note: there is too much data being transferred to produce the ultrasound
% image in MATLAB. It is normal for it to "freeze" while collecting data.
% 
% After the taps are complete, you will have the option to save out the
% data. Depending on the imaging parameters (requested depth of image and
% resolution), this can take quite a long time (several minutes or more).
% The command window will tell you once the data have been saved.
% 
% The data are saved in a .mat file. The RF data (which contains the first
% frame of each tap to permit image reconstruction), the IQ data (which
% contains all of the measured plane wave imaging, and the parameters for
% the collection, including TX, TW, Trans, P, PData, CollectionParams,
% Resource, Receive, SeqControl.
% 
% After a collection is complete, it is possible to verify the quality of
% the collected data using the "Tap Check". This does a very shortened
% processing pipeline to get a rough idea of how well the propagating waves
% are being measured. It will open up a new GUI. Select the number of taps
% (which will be evenly spread across the total number of taps e.g. if you 
% collect 20 total taps, and you select 10 taps, you will see tap 1, 3, 5
% 7...17,19 ). Then, click "select ROI" and on the ultrasound image, drag 
% from the most superficial point to the deepest point between which you'd
% like to check the wave propagation. The colourful lines across the image
% show which depths will be used for calculating the wave propagation.
% Adjust the number and location to suit your needs. Then, click "compute".
% Once the calculations are complete, click on a grey area away from a
% button/set of axes. You should now be able to use the arrow keys (left
% and right) to move between specific taps and see the wave propagation. Up
% and down arrow keys will allow you to move between depths. You can also
% use the buttons on the bottom of the app. 
% 
% 

clear all
%% Change here: Parameters Relevant to Ultrasound - Tensiometry Code ------

% imaging parameters
startDepthMM = 6;% mm
endDepthMM = 14;% mm
resDesiredMM = 0.1; % mm of desired image resolution (between 0.05 and 1 mm) 

% between taps
numTaps = 100; % total number of taps, must be even
tapRate = 50; % Hz (between taps)

% within taps
tapSamples = 75; % Samples after a tap event
sampleRate = 20000; % Hz (within ultrasound frames for a tap)
%% You should not have to modify past here, unless adapting the code for a different system.

% Add an error check for current working folder.
if ~contains(pwd,'Vantage')
    error('Make sure the current working directory is the Vantage folder.')
end

simulateMode = 0;   % set to 0 to acquire data using Vantage 64LE hardware, 1 for simulation mode

% for saving:
CollectionParams.slowFrameRate = tapRate;
CollectionParams.fastFrameRate = sampleRate;


% --- Initiate time tags ----------------------------------
TimeTagEna = 0; % initial state of time tag function% LW: I suspect this is obsolete - test to make sure

% initiate the start of a trial
flag_tap = 0; % no taps completed; changes to 1 once it has tapped

%% ------- Specify Trans structure array ---------------------------------
Trans.name = 'CL15-7','035C';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);    % L11-5v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;      % set maximum high voltage limit for pulser supply.

%% ------- System parameters -------------------------------------------------
filename = 'CL15-7vFlashHFR_acq'; % used to launch VSX automatically
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 64;   % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = simulateMode;
waveLengthMM = 1000*Resource.Parameters.speedOfSound/(10^6*Trans.frequency);

%% Define imaging attributes

P.startDepth = round(startDepthMM/waveLengthMM);  % Acquisition depth in wavelengths
P.endDepth = round(endDepthMM/waveLengthMM);   % This should preferrably be a multiple of 128 samples.
P.numAcqs = tapSamples;      % no. of Acquisitions in a Receive frame (this is a "superframe")
P.numFrames = numTaps;      % no. of Receive frames (real-time images are produced 1 per frame)


% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, round(resDesiredMM/waveLengthMM,2)];% [ x y z] automatically converts z depth to wavelengths
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.
if PData(1).Size(1) * tapSamples * tapRate * numTaps*10^-6 > 150 %
    error('You may be requesting a VERY large matrix because of the resolution and/or number of taps requested. Verify that these values are correct.')
end
% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

%% Specify Resources.
reducedFrames = 10; % in the non-collection time
tframe=zeros(100000,1);
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 1280*tapSamples;   % 4096 - this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = numTaps;       % number of 'taps'

Resource.RcvBuffer(2).datatype = 'int16';% add a second buffer for non collection
Resource.RcvBuffer(2).rowsPerFrame = 1280*tapSamples;   % 4096 - this size allows for maximum range
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = reducedFrames;       % number of plane wave imaging frames

% specify the tensiometry data buffers - I/Q data
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).pagesPerFrame = tapSamples;  % set up the complex buffer
Resource.InterBuffer(1).numFrames = numTaps;  % number of taps

% specify the initial imaging buffer used to align the probe with the
% tissue. Use interBuffer to convert to Image Buffer
Resource.InterBuffer(2).datatype = 'complex'; % set up the complex buffer
Resource.InterBuffer(2).numFrames = 1;  % used to move data over to the image buffer

Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 1;

Resource.DisplayWindow(1).Title = [Trans.name,' HFR Tensiometry'];
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);
Resource.Parameters.initializeOnly = 0;

%% Specify waveform parameters: TX, TW, TGC

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
TX(1).waveform = 1;            % use 1st TW structure.
TX(1).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(1).focus = 0;
TX(1).Steer = [0.0,0.0];       % theta, alpha = 0.
TX(1).Apod = ones(1,Trans.numelements);
TX(1).Delay = computeTXDelays(TX(1));

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,298,395,489,618,727,921,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 1), 1,numTaps*tapSamples+reducedFrames); 

% - Set event specific Receive attributes for the tensiometry section
for i = 1:numTaps
    for j = 1:tapSamples
        % -- Acquisitions 
        rcvNum = tapSamples*(i-1) + j;
        Receive(rcvNum).Apod(1:2:Trans.numelements)=1; % 64LE has only 64 receive channels
        Receive(rcvNum).framenum = i;
        Receive(rcvNum).acqNum = j;

    end
end

for i = numTaps*tapSamples+1 : numTaps*tapSamples + reducedFrames % for the non-collection frames; indexed as after the tensiometry frames
   
        fr2 = i-numTaps*tapSamples;
        rcvNum = i;
        
        Receive(rcvNum).Apod(1:2:Trans.numelements)=1; % 64LE has only 64 receive channels
        Receive(rcvNum).bufnum = 2;
        Receive(rcvNum).framenum = fr2; % resets the frame number for the second buffer

end


%% Specify Recon structure arrays.
Recon = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:tapSamples), 1, numTaps  + reducedFrames ); % one set of recon for IQ transfer, one for a single image

           
for i = 1:numTaps
    
%     Recon(i).rcvBufFrame = i;
    Recon(i).IntBufDest = [1, i]; % for each tap, add to the next increment of the interbuffer
    Recon(i).ImgBufDest = [0,0]; % no image to recreate here
    
    k = (i-1)*tapSamples + 1;
    Recon(i).RINums = k : k+tapSamples-1;
    
end


for i = numTaps+1:numTaps+reducedFrames
    Recon(i).RINums = numTaps*tapSamples + i - numTaps;
    Recon(i).IntBufDest = [2,1]; % send to the complex interbuffer (the second)
%     Note that the structure is initially set up with the image buffer as
%     [1,1] and the tap loop above removes it. No need to re-specify in
%     this loop.
end


%% Define ReconInfo structures.  -> one for every tapsample going to IQ buffer 
ReconInfo = repmat(struct('mode', 'replaceIQ', ...
                   'txnum', 1, ...
                   'rcvnum', 1, ...   
                   'regionnum', 1,...
                   'pagenum',1),1,tapSamples*numTaps + reducedFrames); % one for every IQ frame transferred and an extra set (reducedFrames) for the image processing
                  
% Reconstruction information for the tensiometry. Define the pages and
% receive structures to reference for each tap/sample.
for i = 1:numTaps
   for j = 1:tapSamples
       
       k = (i-1)*tapSamples + j;
       
       ReconInfo(k).pagenum = j;
       ReconInfo(k).rcvnum = k;
   end
end

% Reconstruction information for the image - hence change to
% "replaceIntensity" instead of "replace IQ" for tensiometry
 nRIImage = tapSamples*numTaps+1 ;
for i = nRIImage : tapSamples*numTaps + reducedFrames
   
    ReconInfo(i).mode = 'replaceIntensity';
    ReconInfo(i).rcvnum = i ;
end



%% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMode','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

                     
% external processing function to process time tag data
Process(2).classname = 'External';
Process(2).method = 'readTimeTag';
Process(2).Parameters = {'srcbuffer','receive',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1, ...
                         'dstbuffer','none'};


% EF3 is external function for tap check
Process(3).classname = 'External';
Process(3).method = 'tapCheckFunction';
Process(3).Parameters = {'srcbuffer','inter',... % name of buffer to process.
    'srcbufnum',1,...
    'srcframenum',0,...
    'dstbuffer','none'};
% 
Process(4).classname = 'External';
Process(4).method = 'saveRFData';
Process(4).Parameters = {'srcbuffer','receive',... % name of buffer to process.
    'srcbufnum',1,...
    'srcframenum',0,...
    'dstbuffer','none'};
                   
Process(5).classname = 'External';
Process(5).method = 'saveIQData';
Process(5).Parameters = {'srcbuffer','inter',... % name of buffer to process.
    'srcbufnum',1,...
    'srcframenum',0,...
    'dstbuffer','none'};

%% Specify SeqControl structure arrays.

SeqControl(1).command = 'sync';%'jump'; % jump back to start.

SeqControl(2).command = 'timeToNextAcq';  % time between tap samples
SeqControl(2).argument = 1/sampleRate*1E6; % 50 usecs corresponds to a 20 kHz frame rate
FR=1/(SeqControl(2).argument*1e-6);

SeqControl(3).command = 'triggerOut';
SeqControl(3).argument = 0;

SeqControl(4).command = 'returnToMatlab';

SeqControl(5).command = 'timeToNextAcq';  % time between taps - the time taken to collect the samples
SeqControl(5).argument = 1/tapRate*1E6 - SeqControl(2).argument*tapSamples;  % 5 msec

SeqControl(6).command = 'timeToNextAcq';  % time between imaging frames
SeqControl(6).argument = 1/(tapRate)*1E6; % 5 msec

SFR=1/(tapSamples*SeqControl(2).argument*1e-6+SeqControl(5).argument*1e-6);

nsc = length(SeqControl)+1; % nsc is count of SeqControl objects
n = 1; % n is count of Events
% ntaps = 0; % number of taps


%% Starting EVENT sequence: initial imaging, no tap

nStart = 1;

for i = 1 : reducedFrames
   
    
        Event(n).info = 'Acquire RF';
        Event(n).tx = 1;
        Event(n).rcv = nRIImage+i-1; %tapSamples*(i-1) + j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 6;
        n = n+1;
      
    % Set last acquisitions SeqControl for transferToHost.
    % Do reconstruction and processing for last sub frame
     Event(n-1).seqControl = [6,nsc,4];
     Event(n-1).recon = 0;
     Event(n-1).process = 0;
    SeqControl(nsc).command = 'transferToHost'; 
    nsc = nsc + 1;

    
    Event(n).info = 'reconstruct image'; % noop between frames for frame rate control
    Event(n).tx = 0; % no transmit
    Event(n).rcv = 0; % no rcv
    Event(n).recon = numTaps + i;% nReconImage; % 
    Event(n).process = 1; % external processing function
    Event(n).seqControl = 0; %
    n=n+1;

 
end

    
Event(n).info = 'Jump back';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = nsc;
SeqControl(nsc).command = 'jump';
SeqControl(nsc).argument = nStart; % this forces it to loop back and start this reduced memory load

nsc = nsc+1;
n = n+1;
%% Collect the ultrasound-tensiometry

nCollect = n;


for i = 1:numTaps 
    
%     ntaps=ntaps+1;
    
    for j = 1:tapSamples
    
        Event(n).info = 'Acquire RF';
        Event(n).tx = 1;
        Event(n).rcv = tapSamples*(i-1) + j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        if (j==1)
            Event(n).seqControl = [3,2];
        end
        n = n+1;
        

    end
    
    % Set last acquisitions SeqControl for transferToHost.
    % Do reconstruction and processing for last sub frame
     Event(n-1).seqControl = [5,nsc,4];
     Event(n-1).recon = 0;
     Event(n-1).process = 0;
    SeqControl(nsc).command = 'transferToHost'; % transfer all acqs in one super frame
    nsc = nsc + 1;

    
    Event(n).info = 'reconstruct image'; % noop between frames for frame rate control
    Event(n).tx = 0; % no transmit
    Event(n).rcv = 0; % no rcv
    Event(n).recon = i;%nReconImage; % 
    Event(n).process = 1; % external processing function
    Event(n).seqControl = 0; %
    n=n+1;

end

Event(n).info = 'save RF data'; %
Event(n).tx = 0; % no transmit
Event(n).rcv = 0; % no rcv
Event(n).recon = 0;%nReconImage; %
Event(n).process = [4]; % external processing function for the tap
Event(n).seqControl = 0; %
n=n+1;

Event(n).info = 'save IQ data'; %
Event(n).tx = 0; % no transmit
Event(n).rcv = 0; % no rcv
Event(n).recon = 0;%nReconImage; %
Event(n).process = [5]; % external processing function for the tap
Event(n).seqControl = 0; %
n=n+1;

Event(n).info = 'Jump back';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = nsc;
SeqControl(nsc).command = 'jump';
SeqControl(nsc).argument = nStart; % this forces it to loop back and start this reduced memory load

nsc = nsc+1;
n = n+1;
%% Events to check the quality of the tap

nCheckTap = n;

Event(n).info = 'check tap'; %
Event(n).tx = 0; % no transmit
Event(n).rcv = 0; % no rcv
Event(n).recon = 0;% no recon
Event(n).process = 3; % external processing function for the tap
Event(n).seqControl = 0; %
n=n+1;

Event(n).info = 'Jump back';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = nsc;
SeqControl(nsc).command = 'jump';
SeqControl(nsc).argument = nStart; % this forces it to loop back and start this reduced memory load

nsc = nsc+1;
n = n+1;


%% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Range Change
% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% - Time Tag
UI(3).Control =  {'UserB5','Style','VsSlider','Label','Time Tag',...
                  'SliderMinMaxVal',[0,2,TimeTagEna],...
                  'SliderStep',[0.5, 0.5],'ValueFormat','%1.0f'};
UI(3).Callback = text2cell('%TImeTagCallback');

UI(4).Control = {'Style','pushbutton',...
    'String','Start',...
    'Units','normalized',...
    'Position',[0.6875,0.5000,0.12,0.04],...
    'FontUnits','normalized',...
    'FontSize',0.5,...
    'BackgroundColor',[0,0.8,0.2],...
    'Callback',{@startTrial}};
UI(4).Callback = text2cell('%startTrial');

% save button
UI(5).Control = {'Style','pushbutton',...
    'String','Save',...
    'Units','normalized',...
    'Position',[0.6875+0.12,0.5000,0.12,0.04],...
    'FontUnits','normalized',...
    'FontSize',0.5,...
    'BackgroundColor',[0.5,0.1,0.1],...
    'Callback',{@saveButton}};
UI(5).Callback = text2cell('%saveButton');

% TapCheck for  quality
UI(6).Control = {'Style','pushbutton',...
    'String','Tap Check',...
    'Units','normalized',...
    'Position',[0.6875,0.400,0.20,0.04],...
    'FontUnits','normalized',...
    'FontSize',0.5,...
    'BackgroundColor',[0,0.5,0.5],...
    'Callback',{@TCUI}};
UI(6).Callback =text2cell('%TCUI');


%% External Functions

EF(1).Function = text2cell('%tapCheckFunction'); % this will create the function written below and make it into a function
EF(2).Function = text2cell('%saveRFData');
EF(3).Function = text2cell('%saveIQData');


%% Save all the structures to a .mat file.
% and invoke VSX automatically
filename=['MatFiles/',filename];
save(filename);

disp([ mfilename ': NOTE -- Running VSX automatically!']), disp(' ')

VSX
commandwindow  % just makes the Command window active to show printout


return


%% **** Callback routines to be converted by text2cell function. ****
%SensCutoffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutoffCallback

%RangeChangeCallback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        P.endDepth = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);

evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback


%TImeTagCallback
import com.verasonics.hal.hardware.*
TimeTagEna = round(UIValue);
VDAS = evalin('base', 'VDAS');
switch TimeTagEna
    case 0
        if VDAS % can't execute this command if HW is not present
            % disable time tag
            rc = Hardware.enableAcquisitionTimeTagging(false);
            if ~rc
                error('Error from enableAcqTimeTagging')
            end
        end
        tagstr = 'off';
    case 1
        if VDAS
            % enable time tag
            rc = Hardware.enableAcquisitionTimeTagging(true);
            if ~rc
                error('Error from enableAcqTimeTagging')
            end
        end
        tagstr = 'on';
    case 2
        if VDAS
            % enable time tag and reset counter
            rc = Hardware.enableAcquisitionTimeTagging(true);
            if ~rc
                error('Error from enableAcqTimeTagging')
            end
            rc = Hardware.setTimeTaggingAttributes(false, true);
            if ~rc
                error('Error from setTimeTaggingAttributes')
            end
        end
        tagstr = 'on, reset';
end
% display at the GUI slider value
h = findobj('Tag', 'UserB5Edit');
set(h,'String', tagstr);
assignin('base', 'TimeTagEna', TimeTagEna);
return
%TImeTagCallback


%SaveMyRFData
    % Extract the time stamps associated with each of the superframes
    % Read RcvData buffer 1, and extract the time tag value from the first two
    % samples from channel 1 in each frame (note this would be the time tag
    % from only the first acquisition of a frame that included multiple
    % acquisitions).  Create an array of time stamp values in seconds for all
    % frames in the receive buffer
    numFrames = size(RcvData{1}, 3);  % number of frames in Rcv buffer 1 in Matlab Workspace
    % create a column vector for time tag value from first acquisition of each frame
    TimeTagValues = zeros(numFrames, 1);
    W = zeros(2, 1); % W will be the two time tag words after conversion to unsigned integer values
    for nn = 1:numFrames
        for in=1:2
            % read first and second sample from column 1 of the frame, and
            % convert to double
            W(in) = double(RcvData{1}(ii, 1, nn));
            if W(ii) < 0
                % translate 2's complement negative values to their unsigned integer equivalents
                W(ii) = W(ii) + 65536;
            end
        end
        % first sample is 16 LSB's of time tag value and second sample is 16
        % MSb's so scale and sum to get actual time tag value
        TimeTagValues(nn, 1) = W(1) + 65536 * W(2);
    end
    % the 32 bit time tag counter increments every 25 usec, so we have to scale
    % by 25 * 1e-6 to convert to a value in seconds
    TimeTagValues = TimeTagValues/4e4;
    
%     % Now extract all the RF data
% 
    k=1;
    myRF=zeros(tapsamples,Resource.Parameters.numRcvChannels,numTaps,'int16');
    for i=1:numTaps
        for j=1:tapSamples
            j1=Receive(k).startSample;
            j2=Receive(k).endSample;
            myRF(:,:,j)=RcvData{1}(j1:j2,1:2:end,i); % harcoding that data is collected from every other element
            k=k+1;
        end
    end
% % display at the GUI slider value
% h = findobj('Tag', 'UserB5Edit');
% set(h,'String', tagstr);
% assignin('base', 'TimeTagEna', TimeTagEna);
return
%SaveMyRFData

%startTrial
startTrial(varargin)

% change the start event
nCollect = evalin('base','nCollect'); % send it to the sync start
%     disp('Waiting for sync. ');
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',nCollect};
evalin('base',['Resource.Parameters.startEvent =',num2str(nCollect),';']);
assignin('base','Control',Control);
evalin('base','flag_tap = 1;')
%startTrial


%TCUI
TCUI(varargin)

% this is the UI callback (i.e. button is clicked, send it to the
% processing external function that ACTUALLY checks the tap)
% change the start event
nCheckTap = evalin('base','nCheckTap'); % send it to the sync start
%     disp('Waiting for sync. ');
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',nCheckTap};
evalin('base',['Resource.Parameters.startEvent =',num2str(nCheckTap),';']);
assignin('base','Control',Control);
%TCUI


%tapCheckFunction
tapCheckFunction(varargin)

w = 1:size(varargin{1},2);
IData = squeeze(varargin{1}(:,w,1,:,:));
QData = squeeze(varargin{2}(:,w,1,:,:));
% getParamsTapCheck(IData,QData)
Trans = evalin('base','Trans');
SeqControl = evalin('base','SeqControl');
Resource = evalin('base','Resource');
wx = 0:Trans.spacingMm:(Trans.numelements-1)*Trans.spacingMm;
CollectionParams = evalin('base','CollectionParams');

Fprobe = Trans.frequency;
Ffast = CollectionParams.fastFrameRate;
Fslow = CollectionParams.slowFrameRate;

tapCheck_change(IData, QData, Fprobe, Ffast, Fslow, wx)   
%tapCheckFunction

%saveButton
saveButton(varargin)


Receive = evalin('base','Receive');
Resource = evalin('base','Resource');
numTaps = evalin('base','numTaps');

TX = evalin('base','TX');
TW = evalin('base','TW');
Trans = evalin('base','Trans');
P = evalin('base','P');
PData = evalin('base','PData');
CollectionParams = evalin('base','CollectionParams');
SeqControl = evalin('base','SeqControl');
RF = evalin('base','RF');
I = evalin('base','IData');
Q = evalin('base','QData');
[file,path]  = uiputfile('G:\My Drive\NIMBLE\Projects\ultraFAST\Data\23032022RightAnkle\BundleTest\*.mat');
save (fullfile(path,file), 'RF','I','Q','TX','TW','Trans','P','PData','CollectionParams','Resource','Receive','SeqControl','-v7.3');
fprintf('File saved to: %s\n',fullfile(path,file));

%saveButton


%saveRFData
saveRFData(varargin)
% size(varargin)

    disp('RF Processing!')
%     msg = 'hi'
%     evalin('base','pass_msg = msg')
  
   flag_tap = evalin('base','flag_tap');
%     if flag_tap == 1
    Receive = evalin('base','Receive');
    Resource = evalin('base','Resource');
    numTaps = evalin('base','numTaps');

    j1=Receive(1).startSample;
    j2=Receive(1).endSample;
    % NOTE: we are only saving one RF image per tap
    RF=zeros(j2-j1+1,Resource.Parameters.numRcvChannels,numTaps,'int16');
    for i=1:numTaps
        for j=1%:tapSamples
            j1=Receive(j).startSample;
            j2=Receive(j).endSample;
            RF(:,:,i)=varargin{1}(j1:j2,1:2:end,i); % hardcoding that data is collected from every other element
        end
    end
    assignin('base','RF',RF)
    
    disp('RF Processing Complete!')
%     else
%         warndlg('No data has been collected.')
%     end
%saveRFData

%saveIQData
saveIQData(varargin)
disp('IQ Processing!')
% if flag_tap == 1
    IData = squeeze(varargin{1}(:,:,1,:,:));
    QData = squeeze(varargin{2}(:,:,1,:,:));

    assignin('base','IData',IData)
    assignin('base','QData',QData)
    disp('IQ Processing complete!')
% end
% I = evalin('base','IData');
% I = evalin('base','QData');
%saveIQData