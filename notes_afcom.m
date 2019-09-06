function [ myscreen ] = notes_afcom( varargin )
%%
% Step-by-step reading of afcom with my comments
%
% varargin: Variable length input argument list.
% Allows any number of arguments to a function.
% The variable varargin is a cell array containing the optional arguments to the function.
% varargin must be declared as the last input argument and collects all the inputs from that point onwards.
%
%
% Shih-Wei Wu

global stimulus fixStimulus

stimulus = struct; 
fixStimulus = struct;

stimulus.rotSpd = 90; %what's this?

%% Initialize Variables

% add arguments later
% why define these now and still calls getArgs that will define most of these?
scan = 0;
plots = 0;
noeye = 0;
debug = 0;
replay = 0;
powerwheel = 0;
run = 0; 
eyewindow=0; 
mouse=0; 
practice=0; 
practiceType=-1;
cue=0;

getArgs(varargin,{'scan=0','cue=1','plots=0','noeye=0','powerwheel=1','eyewindow=3','practice=0','practiceType=-1','debug=0','replay=0','run=0','build=0','mouse=0'});
% getArgs(varargin,<validVars>,<verbose=1>)
% passed in varargin, creates variables in the calling function determined by the arguments
% validVars is a cell array of names that will be checked against to see if
% the variables that are being set are valid variable names. With a
% validVars, it will set defaults, and complain if a variable outside of
% the list is set. 
%
% So
% scan,cue,plots,noeye,powerwheel,eyewindow,practice,practiceType,debug,reply,run,build,mouse
% are the variables you can specify in the input arguments. If the
% variables are not specified in the input arguments, they will be set by
% the default specified here.
%
%
stimulus.scan = scan; % scanning experiment? default is no.
stimulus.plots = plots; %plot performance after experiment? default is no.
stimulus.noeye = noeye; %no eye tracking? default is no, which means there is eye tracking?
stimulus.cue = cue; % cue = 1 means direction, cue = 2 means color
stimulus.practice = practice; %practice? default is no.
stimulus.practiceType = practiceType; %practiceType? default is -1. what is this?
stimulus.mousedebug = mouse; %mouse? default is no (not using mouse). 
stimulus.powerwheel = powerwheel; %powerwheel? default is yes (turning the knob to make a response)
stimulus.eyewindow = eyewindow; %eyewindow? not sure what this is for. default is 3.
stimulus.debug = debug; %debug? default is no. 
% debug is a MATLAB built-in function. so this is problem.
stimulus.replay = replay; %replay? not sure what it is for. default is no.
stimulus.overrideRun = run; %run? not sure what this is for. default is no. 
% run is a MATLAB built-in function. so problem.

clear localizer invisible scan noeye task test2 attend build eyewindow mouse practice powerwheel cue

if stimulus.scan
    warning('Disabling eyewindow');
    stimulus.eyewindow=0;
end

%% Open Old Stimfile
if ~stimulus.replay %replay default is 0. so by default we could go into this if statement. 
    stimulus.counter = 1;
    
    if ~isempty(mglGetSID) && isdir(sprintf('~/data/afcom/%s',mglGetSID))
        % Directory exists, check for a stimfile
        % mglGetSID: Gets whatever SID has been set by mglSetSID
        % mglSetSID: Sets the subject ID which is used by initScreen
        % and mglTaskLog. Subject ID is then retrieved with mglGetSID
        % 
        % So there should be a call to mglSetSID before this?
        
        files = dir(sprintf('~/data/afcom/%s/1*mat',mglGetSID));
        
        if length(files) >= 1
            fname = files(end).name;
            
            s = load(sprintf('~/data/afcom/%s/%s',mglGetSID,fname));
            % copy staircases and run numbers
            stimulus.counter = s.stimulus.counter + 1;
            stimulus.colors = s.stimulus.colors;
            stimulus.colorwheel = s.stimulus.colorwheel;
            clear s;
            disp(sprintf('(afcom) Data file: %s loaded.',fname));
        else
            warning('(afcom) Unable to load previous data files. If this is *not* the first run there is something wrong.');
        end
    end
end

%% Display run info
stimulus.counter = -1;
if ~stimulus.replay
    disp('*************************');
    disp(sprintf('(afcom) This is scan #%i',stimulus.counter));
    disp('*************************');
end

%% Setup Screen
% what does initScreen do?
% > it will open up a screen window
% > it will return a structural array called myscreen
if stimulus.replay
    myscreen = initScreen('replayScreen');
else
    myscreen = initScreen('VPixx');
end
% set background to black
myscreen.background = 0;

%% Initialize Stimulus
% what does this mean?
% this is just definining myscreen.stimulusNames{1} and
% stimulus.responseKeys
if ~stimulus.replay
    myscreen.stimulusNames{1} = 'stimulus';
    
    if stimulus.powerwheel
        stimulus.responseKeys = 1;
    else
        stimulus.responseKeys = [1 2]; % left right
    end
else
end

%% load the calib
% Would I need this? this seems like for calibrating the screen? (gamma
% calibration)
if isfield(myscreen,'calibFullFilename')
    calib = load(fullfile(myscreen.calibFullFilename));
    stimulus.calib = calib.calib;
else
    stimulus.calib = []; % need this so that mglLab2rgb doesn't fail
end

%% Colors
% Define stimulus color in rgb.
if ~isfield(stimulus,'colors')
    stimulus.colors.white = [0.8 0.8 0.8]; stimulus.colors.red = [0.8 0 0];
    stimulus.colors.green = [0 1 0]; stimulus.colors.black = [0 0 0];
end

% available range of color/direction increments
stimulus.theta_ = pi/64; % increment size
stimulus.thetas = 0:stimulus.theta_:(2*pi);
    
if 1 %~isfield(stimulus,'colorwheel')
    % get the lab space rgb values
    stimulus.backgroundLab = rgb2lab([0.5 0.5 0.5]);
    % rgb2lab Convert RGB to CIE 1976 L*a*b*
    % need backgroundLab for definining stimulus.colorwheel.rgb
    
    % setup color wheel
    stimulus.colorwheel.acenter = 0;%stimulus.backgroundLab(2);
    stimulus.colorwheel.bcenter = 0;%stimulus.backgroundLab(3);

    % compute the ranges around 0 and pi for the colorwheel

    D = 60;
    stimulus.colorwheel.distanceLab = D;

    stimulus.colorwheel.rgb = zeros(length(stimulus.thetas),3);
    for ti = 1:length(stimulus.thetas)
        theta = stimulus.thetas(ti);
        rgb = ang2rgb(theta);
        % ang2rgb is the helper function (subfunction) at the end of afcom
        % basically it uses ang to define a and b in L*a*b and stimulus.backgroundLab(1) to define the first
        % element of L*a*b. Then it uses mglLab2rgb to transform L*a*b to
        % rgb.

        stimulus.colorwheel.rgb(ti,:) = rgb;
    end

    % if any values are outside RGB space just cut them off
    stimulus.colorwheel.rgb(stimulus.colorwheel.rgb<0) = 0;
    stimulus.colorwheel.rgb(stimulus.colorwheel.rgb>1) = 1;
end

stimulus.colors.mean = [1 1 1]*mean(stimulus.colorwheel.rgb(:));
% stimulus.color.wheel.rgb is a n-by-3 matrix where n is the number of
% color angles (theta) that specifies the rgb of different colors defined by these color angles. 
% This takes the mean of all these rgb values and make it a 1-by-3 vector for mean color. 

%% Sizes
stimulus.fixWidth = 0.5; %fixation width?
stimulus.targetWidth = 10; %will be used to define area to show moving dots
stimulus.patchEcc = 8; %eccentricity?


%% Setup dot patches and stencils

% there will be 12 possible locations where we can show dots. We will use 6
% patches of dots 
% How do you express patches of dots?
%
stimulus.dotScale = 0.3;
stimulus.cueScale = 0.1;

stimulus.dotDirs = [0.5 1.5 0.5 1.5]*pi; % when cue=1 we use these to set the dot directions
% what does 0.5 1.5 mean? 0.5*pi is 90-deg (moving upward), 1.5*pi is 270-deg (moving downward) 
% 4 patches of dots and that's why there are 4 elements in [0.5 1.5 0.5
% 1.5]?
stimulus.dotColors = [0.5 1.5 0.5 1.5]*pi; % when cue=2 we use these to set the color
% horizontal directions:
% stimulus.dotDirs = [0 1 0 1]*pi;
 
dots = struct;
% information about dots, create a a struct array for it

dots.density = 0.2; %some scalar value for dots density. not sure how it is defined.
dots.speed = 3.5; %some scalar value for dots speed. not sure how it is defined
dots.maxAlive = myscreen.framesPerSecond/4; %don't know what it is used for.
dots.maxX = stimulus.targetWidth; %used to define x,y position of dots
dots.maxY = stimulus.targetWidth;

stimulus.dotThetas = [0 0 pi pi]; %don't know what this is for.

for di = 1:4
    stimulus.patches{di} = struct;
    % create a new struct arracy -- stimulus.patches -- under stimulus
    
    % patch dots
    stimulus.patches{di}.dots = initDots(dots);
    % initDots is a helper function (subfunction).
    % initDots:
    % define a patch of dots -- their locations, color, coherence level.
    % output a struct array containing those information.
    % why does dots number have to be even number? because dots will be split into half
    % define dots number based on how big the area is and the dots density
    % set; define dots (x,y) position by randomly placing them within the
    % area (a square defined by maxX, maxY); define coherence level of the
    % dots (set to 1 here); 
    
    % color
    if stimulus.cue==1
        % in this case, all dots in a patch are moving in the same direction. dots either move upward or downward.
        % all dots have the same color (mean of all possible colors)
        stimulus.patches{di}.color = stimulus.colors.mean; %patch color is mean color
        stimulus.patches{di}.dots.dir = stimulus.dotDirs(di); %dots direction is either 0.5*pi (90-deg) or 1.5*pi (270-deg)
    else
        % in this case, all dots are moving horizontally; color is either
        % 0.5*pi or 1.5*pi (in a-b space in L*a*b)
        stimulus.patches{di}.color = ang2rgb(stimulus.dotColors(di));
        stimulus.patches{di}.dots.dir = 0;
    end
    
    % location
    stimulus.patches{di}.theta = stimulus.dotThetas(di);
    stimulus.patches{di}.ecc = stimulus.patchEcc;
    stimulus.patches{di}.xcenter = stimulus.patches{di}.ecc * cos(stimulus.patches{di}.theta); %patch center x 
    stimulus.patches{di}.ycenter = stimulus.patches{di}.ecc * sin(stimulus.patches{di}.theta); %patch center y
end

% stencils are used to control drawing only to specific parts of the screen
% the dots will be drawn on those stencils
mglClearScreen(0); %clear buffer
mglStencilCreateBegin(1); %begin drawing to stencil.
% draw two oval stencils; why using first and third patch location? because
% first and second have the same location, third and fourth pathes have the
% same location.
for di = [1 3]
    %draw a filled oval to stencil that has size targetWidth
    mglFillOval(stimulus.patches{di}.xcenter,stimulus.patches{di}.ycenter,[stimulus.targetWidth, stimulus.targetWidth],[1 1 1]);
end
mglFillOval(0,0,[stimulus.targetWidth stimulus.targetWidth]/4,[1 1 1]); %another stencil at (0,0)
mglFlush; %swap front and back buffer
mglStencilCreateEnd;
% Now that you have created the stencils, later on when you want to draw
% you should call mglStencilSelect(1) before the drawing and
% mglStencilSelect(0) after flushing. This is called in
% drawCue,drawStim,drawTarget

%% Create the cue patch
% dots is a struct array already defined above for dot patches. Here fields below dots are
% redefined for the cue patch. So cue would also show moving dots.
dots.maxX = stimulus.targetWidth/4;
dots.maxY = stimulus.targetWidth/4;
dots.density = 2;
dots.dotScale = 3;
dots.maxAlive = 1000;
stimulus.cueDots = initDots(dots); %initialize dots according to the specifications above. initDots a subfunction
% mostly on locations of these cue dots.
% there are stimulus dots and there are cue dots. Stimulus dots are what appears on the side (left or right)
% cue dots serve as cue and appears at center.

%% Setup Probe Task
% what is it used for? 
task{1}{1} = struct;
% task waits for fixation on first segment
% below looks like an index for order. iti is first (1), fix is second (2) ...
stimulus.seg.iti = 1; 
stimulus.seg.fix = 2;
stimulus.seg.cue = 3;
stimulus.seg.isi = 4;
stimulus.seg.stim = 5;
stimulus.seg.delay = 6;
stimulus.seg.resp = 7;
stimulus.seg.feedback = 8;

% Define min and maximum value (time length in seconds) for different
% segments of a trial
task{1}{1}.segmin = [0 inf 0.5 0.75 inf 1 4 0.75]; %8 elements corresponding to the 8 segments of a trial
task{1}{1}.segmax = [2 inf 0.5 0.75 inf 1 4 0.75];

if stimulus.noeye
    task{1}{1}.segmin(stimulus.seg.fix) = 0;
    task{1}{1}.segmax(stimulus.seg.fix) = 0;
end

if stimulus.practice==1
    task{1}{1}.segmin(stimulus.seg.cue) = 1;
    task{1}{1}.segmax(stimulus.seg.cue) = 1;
    task{1}{1}.segmin(stimulus.seg.isi) = 1;
    task{1}{1}.segmax(stimulus.seg.isi) = 1;
    task{1}{1}.segmin(stimulus.seg.resp) = 6;
    task{1}{1}.segmax(stimulus.seg.resp) = 6;
    task{1}{1}.segmin(stimulus.seg.feedback) = 1.5;
    task{1}{1}.segmax(stimulus.seg.feedback) = 1.5;
elseif stimulus.practice==2
    task{1}{1}.segmin(stimulus.seg.cue) = 1;
    task{1}{1}.segmax(stimulus.seg.cue) = 1;
    task{1}{1}.segmin(stimulus.seg.isi) = 1;
    task{1}{1}.segmax(stimulus.seg.isi) = 1;
end

task{1}{1}.waitForBacktick = 1;

task{1}{1}.getResponse = zeros(1,length(task{1}{1}.segmin));
task{1}{1}.getResponse(stimulus.seg.resp) = 1; 
% task{1}{1}.getResponse = [0 0 0 0 0 0 1 0]
% a vector containing 1s and 0s used as index for whether to detect a response

if stimulus.scan==1
    task{1}{1}.numTrials = Inf; %why set number of trials to Inf if scan=1?
else
    task{1}{1}.numTrials = 40;
end

task{1}{1}.random = 1; %what is this for?

task{1}{1}.parameter.trialType = [1 1 1 2 2 2 0 0 3 4]; % 1 = spatial, 2 = feature, 0 = no cue, 3 = exact cue (1+2), 4 = target only
task{1}{1}.parameter.target = [1 2 3 4]; % which patch is the target
task{1}{1}.parameter.duration = [0.25 1.0]; % bump to 0.25/0.50/1.00 for full task? 

if stimulus.practice==1
    task{1}{1}.parameter.duration = 1.0;
end

if stimulus.practiceType>=0
    task{1}{1}.parameter.trialType= stimulus.practiceType;
end

task{1}{1}.parameter.cue = stimulus.cue; % which cue condition, 1=direction cues, 2=color cues

if ~stimulus.replay && stimulus.scan
    task{1}{1}.synchToVol = zeros(1,length(task{1}{1}.segmin));
    task{1}{1}.synchToVol(end) = 1;
end

% feature target; all are nan for now.
task{1}{1}.randVars.calculated.dead = nan;
task{1}{1}.randVars.calculated.targetAngle = nan; % angle of the target
task{1}{1}.randVars.calculated.distractorAngle = nan; % angle of the other thing you had to attend
task{1}{1}.randVars.calculated.distractor = nan;
task{1}{1}.randVars.calculated.featdist = nan; % number of the matched feature
task{1}{1}.randVars.calculated.sidedist = nan; % number of the matched side
task{1}{1}.randVars.calculated.distdist = nan; % number of the one you ignored (not matched side or matched feature)
task{1}{1}.randVars.calculated.angle1 = nan;
task{1}{1}.randVars.calculated.angle2 = nan;
task{1}{1}.randVars.calculated.angle3 = nan;
task{1}{1}.randVars.calculated.angle4 = nan;
task{1}{1}.randVars.calculated.respAngle = nan;
task{1}{1}.randVars.calculated.respDistance = nan;
task{1}{1}.randVars.calculated.distDistance = nan;
task{1}{1}.randVars.calculated.cwOffset = nan; % colorwheel offset rotation

%% Mouse movement storage data

% average reaction time is ~300, but the matrix will get filled with zeros
% (bad) if we don't pre-fill it with nan
stimulus.data.mouseTrack = nan(min(task{1}{1}.numTrials,50),500); % a n-by-50 matrix of NaN where n is number of trials (40 or 50)
stimulus.data.mouseTick = 1;


%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,@endTrialCallback,[]);
    % initTask is a mgl function that initializes task for stimuli programs
    %
    % initTask takes input task{1}{phaseNum} and myscreen and calls a bunch of functions that take these two as inputs
    % and return an update on task{1}{phaseNum} and myscreen.
    %
    % 1. startSegmentCallback: does nothing.
    % 2. screenUpdateCallback: does all the drawing and get response  
    % 3. getResponseCallback: if task.thistrial.gotResponse==0 (didn't get
    % a response?), then task = jumpSegment(task) (jump to the next segment).
    % 4. startTrialCallback: set mouse position, get mouse position, cur dots direction, distractors, targetAngle, distractorAngle ... 
    % 5. endTrialCallback: if task.thistrial.dead, return, end. 
    % 
    % are the callback functions listed here. They have to be in this
    % particular order. See initTask. In initTask, the sixth callback
    % function is startBlockCallback, the seventh callback function is
    % randCallback
    % 
    %
    %
    % These functions are defined at the end as subfunctions. 
    % All these functions take task{1}{phaseNum} and myscreen as inputs
    %
    %
    % > what does "callbacks" mean?
    % Usually use together with function: callback function
    % e.g. plot(x,y,'ButtonDownFcn',@lineCallback)
    % lineCallBack is a callback function for the lines created by the plot
    % function. Use the @ operator to assign the function handle to the
    % ButtonDownFcn property of each line created by plot
    % 
    % function lineCallback(src,~)
    % src.Color = 'red';
    % end
    %
    % > what does "function handles" mean?
    %     FUNHANDLE = @FUNCTION_NAME returns a handle to the named function,
    %      FUNCTION_NAME. A function handle is a MATLAB value that provides a
    %      means of calling a function indirectly. You can pass function
    %      handles in calls to other functions (which are often called function
    %      functions). You can also store function handles in data structures for
    %      later use (for example, as Handle Graphics callbacks). A function
    %      handle is one of the standard MATLAB data types. Its class is
    %      'function_handle'.
    %
    % To create a function handle, precede the function mane with an @ sign.
    % e.g. create a handle named f as follows:
    % f = @myfunction
    % You call a function using a handle the same way you call the function
    % directly.
    % e.g. 
    % function y = computeSquare(x)
    % y=x.^2;
    % end
    %
    % f = @computeSquare
    % a=4;
    % b = f(a)
    % 
    % function handles are variables you can pass to other functions
    % e.g. q = integral(f,0,1)
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~stimulus.replay && ~stimulus.noeye
    myscreen = eyeCalibDisp(myscreen);
    % eyeCalibDisp is a mgl function that runs eye calibration routine
    
    % let the user know
    disp(sprintf('(afcom) Starting run number: %i.',stimulus.counter));
end


%% Main Task Loop

phaseNum = 1;
% Again, only one phase. 
% What does the above sentence mean? only one phase?
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
    % updateTask is mgl function that runs experimental tasks. It updates
    % task in running stimulus programs. 
    % It updates task{1} and phaseNum, both are used to evaluate the while
    % loop.
    %
    % So what is the rationale here? 
    % >> in initTask you call functions that takes care of evertyhing you ever need 
    % (including drawing stimuli and detecting a response) in this experiment. 
    % >> initTask does the above by calling the functions needed to perform
    % these things "if" these functions are called:
    %   >> startSegmentCallback: does nothing.
    %   >> trialResponseCallback: 
    %   >> screenUpdateCallback: refreshes the screen. draw stimuli
    %   depending on which segment of a trial you are in. e.g. case
    %   stimulus.seg.stim, drawStim and drawFix.
    %   >> endTrialCallback
    %   >> startTrialCallback
    %   >> startBlockCallback
    %   >> randCallback
    %   ...
    % These functions are subfunctions at afcom_wei.
    %
    % > in updateTask you show them? 
    % > in updateTask you detect response?
    %
    % so all the drawings are done in updateTask?
    %
    % I see drawFix, drawBorder, drawAllBorders, drawCue, drawStim, drawTarget, afPoints, drawPicker, drawResp, mglGluPartialDisk_ 
    % But where are these functions called?
    % drawStim, drawFix, drawCue, drawTarget, drawPicker, drawAllBorders
    % called in screenUpdateCallback. screenUpdateCallback is called by
    % initTask.
    %
    % flip screen
    myscreen = tickScreen(myscreen,task);
    % why flip screen? to show the drawing on screen.
    % tickScreen calls mglFlush to "flip" screen (swap front and back
    % buffer)
    % front buffer: what is being shown on screen
    % back buffer: where you are drawing 
    % front buffer is displayed on screen and you draw to the back buffer,
    % then you swap them when you are done drawing so the image in the back
    % buffer is shown. This is so the user never sees a partially drawn
    % image as you are updating the screen -- they see the last complete
    % image, i.e. the front buffer.
end

% task ended
mglClearScreen;
mglTextSet([],32,stimulus.colors.white);
% get count
mglTextDraw('Please wait',[0 0]);
mglFlush %swap front and back buffer (waits for one frame tick)
myscreen.flushMode = 1;

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);
% endTask is mgl function that packages all variables into myscreen and
% reports any thrown errors
%
if ~stimulus.replay && stimulus.plots
    disp('(afcom) Displaying plots');
    dispInfo(stimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%
function dispInfo()
% This function shows information about subjects' performance after they
% finish the experiment.
% This is called at the end of the experiment:
% when stimulus.replay is false and when stimulus.plots is
% true: dispInfo(stimulus)
%
% This function doesn't have any input. But when it is called, it does have stimulus
% as input. This is weird because MATLAB should complain there are too many input arguments 
% in this case.
%
% This function uses stimulus, myscreen and task. I guess all these are saved in data file.
% 

files = dir(fullfile('~/data/afcom/',mglGetSID,'*.mat'));
% mglGetSID: gets whatever SID has been set by mglSetSID
% mglSetSID: sets the subject ID which is used by initScreen and mglTaskLog
% what does it mean to have more than one file for a subject? each file is
% a block?

maxTrackLength = 0;
for fi = 1:length(files)
    load(fullfile('~/data/afcom/',mglGetSID,files(fi).name));
    % load subject's data file. what is in it? myscreen and task? and stimulus
    
    exp = getTaskParameters(myscreen,task); %need myscreen and task as inputs. must be from the file loaded...
    % get the task parameters, reaction times etc, out of the screen and
    % task variables
    
    % Setting up cell arrays e,mt,maxTrackLength
    e{fi} = exp{1};
    mt{fi} = stimulus.data.mouseTrack(1:e{fi}.nTrials,:);
    mt{fi}(mt{fi}==0) = nan;
    maxTrackLength = max(maxTrackLength,size(mt{fi},2));
end

% clear duration
% warning('adding duration = 1 if missing');
% for ei = 1:length(e)
%     if ~isfield(e{ei}.parameter,'duration')
%         e{ei}.parameter.duration = ones(size(e{ei}.parameter.trialType));
%     end
% end

%% concatenate all trials
% This probably means defining a bunch of variables each a 1xn vector where
% each n indicates the number of runs
pvars = {'target','trialType','cue','duration'}; %pvars probably means task parameters 
rvars = {'dead','targetAngle','distractorAngle','angle1','angle2','angle3',...
    'angle4','respAngle','respDistance','distDistance'}; %rvars probably means "random" parameters (so pvars is fixed?)
runs = [];

% Define variables based on above by creating an empty matrix for each
for pii = 1:length(pvars)
    eval(sprintf('%s = [];',pvars{pii}));
end
% what this will do is to create empty matrix
% target = [] trialType=[] cue=[] duration=[] ...
for ri = 1:length(rvars)
    eval(sprintf('%s = [];',rvars{ri}));
end
% dead = [] targetAngle=[] distractorAngle=[] angle1=[] angle2=[] ...
% angle3=[],...

runcount = [0 0];
for run = 1:length(e)
    if e{run}.nTrials>0
        runs = [runs ones(1,e{run}.nTrials)];
        runcount(e{run}.parameter.cue(1)) = runcount(e{run}.parameter.cue(1)) + 1;
        for pii = 1:length(pvars)
            eval(sprintf('%s = [%s e{run}.parameter.%s];',pvars{pii},pvars{pii},pvars{pii}));
        end
        for ri = 1:length(rvars)
            eval(sprintf('%s = [%s e{run}.randVars.%s];',rvars{ri},rvars{ri},rvars{ri}));
        end
    end
end
% runs = [40 40 40 40 40] if you have 5 runs each having 40 trials
% runcount = [3 2] probably means there are two cues (therefore 2 columns); each column corresponds to 
% the number of runs for each cue condition.
% target = [target e{}.parameter.target]
% trialType = [trialType e{}.parameter.trialType]
% cue = [cue e{}.parameter.cue]
% ... like that

eval('dur = duration;');
% Why do this?

%% concatenate mouse tracks
amt = nan(length(target),maxTrackLength); %mt for mouse track
start = 1;
for run = 1:length(e) %length of e should be number of runs
    stop = (start+e{run}.nTrials-1); %where to stop
    amt(start:stop,1:size(mt{run},2)) = mt{run}; %mt{run} is a matrix based on stimulus.data.mouseTrack
    start = stop + 1;
end

%% go backward through mouseTracks and fix jumps
% What are jumps? Why is this an issue?
% assume that you end near zero, so if you jump -pi you need to -pi the
% earlier section, etc
amt = fliplr(amt); %flip array in left/right direction
% so ros is trial, column is time? by flipping, last column becomes the
% first.
for ai = 1:size(amt,1) %row means trials?
    track = amt(ai,:); 
    dtrack = diff(track); %difference between adjacent points.
    posidx = find(dtrack>5); %find difference greater than 5. returns their index number.
    negidx = find(dtrack<-5); %find difference smaller than -5
    for pii = 1:length(posidx)
        idx = posidx(pii)+1; % e.g. posidx(1)=5, then idx=5+1=6
        track(idx:end) = track(idx:end)-2*pi; %don't understand. you are updating track using this equation each time there is a difference greater than 5. 
    end
    for nii = 1:length(negidx)
        idx = negidx(nii)+1;
        track(idx:end) = track(idx:end)+2*pi;
    end
    dtrack = diff(track); %why define dtrack again here? it does nothing.
    amt(ai,:) = track;
end
amt = fliplr(amt);

%% create one giant matrix, but just of a few variables that matter
data = [cue' runs' trialType' respDistance' dur']; %these variables are defined in this subfunction. respDistance = e{run}.randVars.respDistance
keepIdxs = ~any(isnan(data(:,4)),2);
data = data(keepIdxs,:);
amt = amt(keepIdxs,:);

disp(sprintf('Total trials: %i',size(data,1)));

%% print out information
disp(sprintf('Runs so far: %i cue direction (cue=1), %i cue color (cue=2)',runcount(1),runcount(2)));

%% plot mousetracks
% step 1: rotate mousetracks so that they are relative to the target
amt_ = amt - repmat(targetAngle(keepIdxs)',1,size(amt,2));
% test plot the average mousetrack
figure; hold on
plot(amt_','-k');
hline(0,'--r');
xlabel('Time from response window start');
ylabel('Rotation (rad)');
drawPublishAxis;

%% plot
% split data by difficulty
edata = data(data(:,5)==1,:);
hdata = data(data(:,5)==0.25,:);

dispInfoFigures(edata,'easy');
dispInfoFigures(hdata,'hard');
% this is a subfunction.

%%
function dispInfoFigures(data,diff)

% build one figure for each task
titles = {'Cue direction: ','Cue color: '};
bins = pi/32:pi/16:pi;
blabels = {};
for bi = 0:(length(bins)-1)
    blabels{bi+1} = sprintf('%i/16',bi);
end

cmap = brewermap(5,'Dark2');

for cue = 1:2
    disp(sprintf('%s cue %s',diff,titles{cue}));
    cdata = data(data(:,1)==cue,:);
    
    disp(sprintf('Trials of: %s so far %i',titles{cue},size(cdata,1)));
    
    all = cdata(cdata(:,3)==0,:);
    disp(sprintf('Type all: %i',size(all,1)));
    spatial = cdata(cdata(:,3)==1,:);
    disp(sprintf('Type spatial: %i',size(spatial,1)));
    feature = cdata(cdata(:,3)==2,:);
    disp(sprintf('Type feature: %i',size(feature,1)));
    target = cdata(cdata(:,3)==3,:);
    disp(sprintf('Type target: %i',size(target,1)));
    base = cdata(cdata(:,3)==4,:);
    disp(sprintf('Type baseline: %i',size(base,1)));

    figure;
    
    group = {'all','spatial','feature','target','base'};
    legends = {'All','Spatial','Feature','Target','Baseline'};
    
    for s = 1:5
        cdat = eval(sprintf('%s(:,4)',group{s}));
        his = hist(cdat,bins);
        his = his/sum(his);
        
        subplot(5,1,s); hold on
        b = bar(bins,his,pi/8);
        set(b,'FaceColor',cmap(s,:),'EdgeColor','w');
        vline(nanmedian(cdat),'--k');
        legend(legends{s});
        ylabel('Proportion (%)');
        xlabel('Response distance from target (target=0');
        set(gca,'XTick',bins,'XTickLabel',blabels);
        drawPublishAxis;
    end
end

%%


%%
function [task, myscreen] = startTrialCallback(task,myscreen)
global stimulus

% swap seglen in
task.thistrial.seglen(stimulus.seg.stim) = task.thistrial.duration;

if stimulus.powerwheel
    mglSetMousePosition(myscreen.screenWidth/2+rand*2*pi*stimulus.rotSpd-pi*stimulus.rotSpd,myscreen.screenHeight/2,1);
else
    mglSetMousePosition(myscreen.screenWidth/2,myscreen.screenHeight/2,2);
end

% get the current mouse position:
mInfo = mglGetMouse(myscreen.screenNumber);
stimulus.live.mouseStart = -mInfo.x/stimulus.rotSpd;

if stimulus.cue==1
    stimulus.cueDots.dir = stimulus.patches{task.thistrial.target}.dots.dir;
else
    stimulus.cueDots.dir = 0; % doesn't matter, dots are incoherent
end

if (task.thistrial.trialType==0) || (task.thistrial.trialType==4)
    task.thistrial.distractor = nan; % no cue, so everything is a potential distractor
else
    if task.thistrial.trialType==2
        distractors = [3 4 1 2];
    else % this could be trial type 1 or 3, but in both cases this is a spatial distractor
        distractors = [2 1 4 3];
    end
    task.thistrial.distractor = distractors(task.thistrial.target);
end

% set the angles of the patches
for di = 1:length(stimulus.patches)
    ctheta = randsample(stimulus.thetas,1);
    
    if di==task.thistrial.target
        task.thistrial.targetAngle = ctheta;
    end
    if di==task.thistrial.distractor
        task.thistrial.distractorAngle = ctheta;
    end
    
    if stimulus.cue==1
        % direction cue, so set the colors to be different
        stimulus.patches{di}.color = ang2rgb(ctheta);
    else
        % color cue, so set the directions to be different
        stimulus.patches{di}.dots.dir = ctheta;
    end
    
    task.thistrial.(sprintf('angle%i',di)) = ctheta;
end

% colorwheel random rotation
task.thistrial.cwOffset = rand*2*pi;
task.thistrial.respAngle = -task.thistrial.cwOffset;

if stimulus.cue==1
    trialTypes = {'nocue','spatial','direction','target','baseline'};
else
    trialTypes = {'nocue','spatial','color','target','baseline'};
end
disp(sprintf('(afcom) Starting trial %i. Attending %s',task.trialnum,trialTypes{task.thistrial.trialType+1}));

% eye tracking 
task.thistrial.dead = 0;
stimulus.live.eyeCount=0;
stimulus.live.fixCount = 0;

% mouse tracking
stimulus.data.mouseTick = 1;
%%


function [task, myscreen] = endTrialCallback(task,myscreen)
%%
if task.thistrial.dead, return; end
% if this is true, then exit function.

respType = {'timeout','click','multiclick','multiclick','multiclick'};
if isnan(task.thistrial.respDistance)
    %if respDistance is not a number
    task.thistrial.respDistance = angdist(task.thistrial.respAngle,task.thistrial.targetAngle);
    task.thistrial.distDistance = angdist(task.thistrial.respAngle,task.thistrial.distractorAngle);
    disp(sprintf('Recorded: %s. angle of %1.2f true %1.2f: %1.2f distance',respType{task.thistrial.gotResponse+1},task.thistrial.respAngle,task.thistrial.targetAngle,task.thistrial.respDistance));
end

function d = angdist(t1,t2)
%%
d = acos(cos(t1)*cos(t2)+sin(t1)*sin(t2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = startSegmentCallback(task,myscreen)
%%
% global stimulus


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Drawing functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function drawCue(task)
%%
global stimulus

switch task.thistrial.trialType
    case 0
        cues = 0;
    case 1
        cues = 1;
    case 2
        cues = 2;
    case 3
        cues = [1 2];
    case 4
        cues = [1 2];
end

if any(cues==0)
    % draw lines to both sides
    dotDirs = unique(stimulus.dotThetas);
    for di = 1:length(dotDirs)
        x = 1.5*stimulus.fixWidth * cos(dotDirs(di));
        y = 1.5*stimulus.fixWidth * sin(dotDirs(di));
        mglLines2(x,y,2*x,2*y,2,stimulus.colors.white);    
    end
end

if any(cues==1)
    % spatial - draw lines to attended locations
        
    % draw the line from fixWidth to 2*fixWidth
    x = 1.5*stimulus.fixWidth * cos(stimulus.patches{task.thistrial.target}.theta);
    y = 1.5*stimulus.fixWidth * sin(stimulus.patches{task.thistrial.target}.theta);
    mglLines2(x,y,2*x,2*y,2,stimulus.colors.white);
end
if any(cues==2)
    % feature - draw the motion direction or the color

    if stimulus.cue==1
        coherence = 1;
        color = stimulus.colors.white;
    elseif stimulus.cue==2
        coherence = 0;
        color = ang2rgb(stimulus.dotColors(task.thistrial.target));
    end
    % cue dots version
    stimulus.cueDots = updateDots(stimulus.cueDots,coherence,false);
    
    mglStencilSelect(1); %before you draw, call this to allow drawing on stencils set earlier.
    afPoints(stimulus.cueDots.x-stimulus.cueDots.maxX/2,stimulus.cueDots.y-stimulus.cueDots.maxY/2,stimulus.cueScale,color);
    mglStencilSelect(0);
end

function drawStim(task,stimSeg)
%%
global stimulus

mglStencilSelect(1); %before you draw, call this to allow drawing on stencils set earlier.
% update and collapse x/y coordinates for drawing
n = stimulus.patches{1}.dots.n;

x = nan(1,n*length(stimulus.patches));
y = x;
r = stimulus.colors.mean(1)*ones(1,n*length(stimulus.patches));
g = r;
b = r;

for di = 1:length(stimulus.patches)
    if task.thistrial.trialType~=4 || ~stimSeg || (stimSeg && (di==task.thistrial.target))
        if stimulus.cue==1 || (stimulus.cue==2 && stimSeg)
            % if this is the actual stim seg and using motion, update
            % coherently
            stimulus.patches{di}.dots = updateDots(stimulus.patches{di}.dots,1,false);
        else
            % otherwise use incoherent motion
            stimulus.patches{di}.dots = updateDots(stimulus.patches{di}.dots,0,false);
        end

        offX = stimulus.patches{di}.xcenter - stimulus.patches{di}.dots.maxX/2;
        offY = stimulus.patches{di}.ycenter - stimulus.patches{di}.dots.maxY/2;

        x(((di-1)*n+1):(di*n)) = offX + stimulus.patches{di}.dots.x;
        y(((di-1)*n+1):(di*n)) = offY + stimulus.patches{di}.dots.y;
        if stimulus.cue==2 || (stimulus.cue==1 && stimSeg) % if this is the actual stimulus segment and we are using color, show the colors
            r(((di-1)*n+1):(di*n)) = stimulus.patches{di}.color(1);
            g(((di-1)*n+1):(di*n)) = stimulus.patches{di}.color(2);
            b(((di-1)*n+1):(di*n)) = stimulus.patches{di}.color(3);
        end
    end
end

drop = isnan(x);
if any(drop)
    x = x(~drop);
    y = y(~drop);
    r = r(~drop);
    g = g(~drop);
    b = b(~drop);
end

% randomly sort x/y/r/g/b so that overlapping patches render correctly
perm = randperm(length(x));
x = x(perm); y = y(perm); r = r(perm); g = g(perm); b = b(perm);

% Draw dots; afPoints is a subfunction that draws dots one at a time.
afPoints(x,y,stimulus.dotScale,[r' g' b']);

% draw all the dots at once
% mglPoints2c(x,y,stimulus.dotScale,r,g,b);

% Drawing dots one at a time vs. drawing all dots at once: which one is
% preferred?

mglStencilSelect(0);

function afPoints(x,y,scale,c)
%%
% draw the dots one at a time with mglGluDisk
% mglGluDisk: draw disk(s) at location x,y.
% 
% This is called in drawTarget,drawCue,drawStim.
% drawStim is where moving dots are drawn.

cFlag = size(c,1)==1;
for di = 1:length(x)
    if cFlag
        mglGluDisk(x(di),y(di),scale,c);
        % scale: size. c: color
    else
        mglGluDisk(x(di),y(di),scale,c(di,:));
    end
end

function drawTarget(task)
%%
global stimulus

if stimulus.cue==1
    stimulus.patches{task.thistrial.target}.dots = updateDots(stimulus.patches{task.thistrial.target}.dots,1,false);
    color = stimulus.colors.mean;
else
    % if we we cued color set the coherence to zero so that there's no
    % direction information
    stimulus.patches{task.thistrial.target}.dots = updateDots(stimulus.patches{task.thistrial.target}.dots,0,false);
    color = ang2rgb(stimulus.dotColors(task.thistrial.target));
end
% compute the offset position
offX = stimulus.patches{task.thistrial.target}.xcenter - stimulus.patches{task.thistrial.target}.dots.maxX/2;
offY = stimulus.patches{task.thistrial.target}.ycenter - stimulus.patches{task.thistrial.target}.dots.maxY/2;

% draw the actual points
mglStencilSelect(1); %before you draw, call this to allow drawing on stencils set earlier.
afPoints(stimulus.patches{task.thistrial.target}.dots.x + offX,stimulus.patches{task.thistrial.target}.dots.y + offY,stimulus.dotScale,color);
mglStencilSelect(0);

function drawPicker(task)
%%
global stimulus

if stimulus.cue==1
    % When we cue spatial/direction we need to draw the color picker
    for ti = 1:length(stimulus.thetas)
        theta = stimulus.thetas(ti) + task.thistrial.cwOffset;
        mglGluPartialDisk_(0,0,1,1.25,180/pi*(theta-stimulus.theta_/2),180/pi*stimulus.theta_,stimulus.colorwheel.rgb(ti,:));
    end
    % Also draw a little marker to indicate the current rotation
    mglGluPartialDisk_(0,0,1,1.25,180/pi*(task.thistrial.respAngle+task.thistrial.cwOffset)-2.5,5,[0.75 0.75 0.75]);
else
    % Don't rotate the marker using cwOffset
    mglGluPartialDisk_(0,0,1,1.25,180/pi*(task.thistrial.respAngle)-2.5,5,[0.75 0.75 0.75]);
end

function drawResp(angle)
%%
global stimulus
% Draw the chosen color as a backgroundc ircle
if stimulus.cue==1
    mglFillOval(0,0,stimulus.fixWidth*[1 1],ang2rgb(angle));
else
    mglGluPartialDisk_(0,0,1,1.25,180/pi*angle-2.5,5,[0.75 0.75 0.75]);
end

function mglGluPartialDisk_(x,y,isize,osize,sangle,sweep,color)
%%
% just a wrapper around mglGluPartialDisk which converst from REAL angles
% to MGL angles. I absolutely hate this aspect of MGL which I assume is
% inherited from OpenGL...
sangle = 90-sangle; % this sets 0 to be vertical and all coordinates go clockwise
mglGluPartialDisk(x,y,isize,osize,sangle,sweep,color);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
%%
global stimulus

mglClearScreen();

if task.thistrial.dead && mglGetSecs(stimulus.live.deadTime)>1
    task = jumpSegment(task,inf);
end

% skip screen updates if you are already dead
if task.thistrial.dead
    if task.thistrial.dead
        mglTextSet([],32,stimulus.colors.red);
        mglTextDraw('Eye Movement Detected',[0 0]);
    end
    return
end

if (task.thistrial.thisseg==stimulus.seg.resp)
    if stimulus.powerwheel
        mInfo = mglGetMouse(myscreen.screenNumber);
        task.thistrial.respAngle = -(mInfo.x-myscreen.screenWidth/2)/stimulus.rotSpd;
    else
        mInfo = mglGetMouse(myscreen.screenNumber);
        degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
        if stimulus.cue==1
            task.thistrial.respAngle = atan2(degy,degx) - task.thistrial.cwOffset;
        else
            task.thistrial.respAngle = atan2(degy,degx);
        end
    end
    task.thistrial.respAngle = mod(task.thistrial.respAngle,2*pi);
    
    stimulus.data.mouseTrack(task.trialnum,stimulus.data.mouseTick) = task.thistrial.respAngle;
    stimulus.data.mouseTick = stimulus.data.mouseTick + 1;
    
    % note that respAngle is stored in *real* angles -- so that it
    % corresponds correctly to the direction task. This means that when you
    % transform into visual space you need to flip into MGL angles, see
    % mglGluDiskAnnulus_ which does this step
end

switch task.thistrial.thisseg
        
    case stimulus.seg.iti
        drawStim(task,false);
        drawFix(task,stimulus.colors.white);
    case stimulus.seg.fix % same as for ITI
        drawStim(task,false);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.cue
        % fixation
        drawStim(task,false);
        drawCue(task);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.isi
        drawStim(task,false);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.stim
        drawStim(task,true);
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.delay
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.resp
        drawTarget(task);
        drawPicker(task);
        if stimulus.cue==1
            % only draw the chosen color at fixation if we're doing cued
            % direction
            drawResp(task.thistrial.respAngle);
        end
        drawFix(task,stimulus.colors.white);
        
    case stimulus.seg.feedback
        drawResp(task.thistrial.targetAngle);
        drawFix(task,stimulus.colors.white);
        
end

drawAllBorders(stimulus.patches,stimulus.targetWidth/2);

% do eye position tracking, but only during some segments
if (~stimulus.noeye) && (stimulus.eyewindow>0) && any(task.thistrial.thisseg==[stimulus.seg.fix stimulus.seg.cue stimulus.seg.stim])
    % check eye pos

    % mouse version for testing with no eyetracker
    if stimulus.mousedebug
        mInfo = mglGetMouse(myscreen.screenNumber);
        degx = (mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
        degy = (mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;

        pos = [degx, degy];
    else
        [pos,~] = mglEyelinkGetCurrentEyePos;
    end
    % compute distance
    dist = hypot(pos(1),pos(2));

    if task.thistrial.thisseg==stimulus.seg.fix
        if stimulus.live.fixCount > stimulus.eyeFrames
            task = jumpSegment(task);
        elseif ~any(isnan(pos))
            if dist < stimulus.eyewindow
                stimulus.live.fixCount = stimulus.live.fixCount + 1;
            else
                stimulus.live.fixCount = 0;
            end
        end
    else
        % Eye movement detection code
        if (~stimulus.noeye) && (stimulus.eyewindow>0) && ~task.thistrial.dead
            if ~any(isnan(pos))
                if dist > stimulus.eyewindow && stimulus.live.eyeCount > stimulus.eyeFrames
                    disp('Eye movement detected!!!!');
                    stimulus.live.deadTime = mglGetSecs;
                    task.thistrial.dead = 1;
                    return
                elseif dist > stimulus.eyewindow
                    stimulus.live.eyeCount = stimulus.live.eyeCount + 1;
                end
            end
        end
    end

end

function [task, myscreen] = getResponseCallback(task, myscreen)
%%
if task.thistrial.dead, return; end

if task.thistrial.gotResponse==0
    % if you don't get a response in this trial
    % jump to the feedback segment
    task = jumpSegment(task);
end

function drawFix(task,color)
%%
global stimulus;

if task.thistrial.dead
    mglGluDisk(0,0,[1 1],stimulus.colors.red,60);
else
    mglFixationCross(stimulus.fixWidth,1,color);
end

function drawBorder(x,y,r,c)
%%
for t = 0:pi/4:2*pi
    mglGluPartialDisk(x,y,r,r+0.05,(180/pi)*(t-pi/16),22.5,c);
end

function drawAllBorders(locations,r)
%%
% draw the borders
for li = 1:length(locations)
    drawBorder(locations{li}.xcenter,locations{li}.ycenter,r,[0.05 0.05 0.05]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rgb = ang2rgb(ang)
%%
global stimulus

a = stimulus.colorwheel.distanceLab*cos(ang)+stimulus.colorwheel.acenter;
b = stimulus.colorwheel.distanceLab*sin(ang)+stimulus.colorwheel.bcenter;

rgb = mglLab2rgb([stimulus.backgroundLab(1) a b],stimulus.calib);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for horizontal motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDots(dots)
%%
dots.dir = 0;

area = dots.maxX*dots.maxY;

dots.n = round(area * dots.density);

% make a some points
% dots.n = 500*dots.density;
% make sure it's an even number. Why dots have to be even number? because
% you want to set them in half
dots.n = dots.n + mod(dots.n,2);

% set half to white and half to black
dots.con = repmat([1 2],1,dots.n/2);

dots.x = rand(1,dots.n)*dots.maxX;
dots.y = rand(1,dots.n)*dots.maxY;

% Why replace dots? Because if you don't then peripheral overlapping dot
% patches will rival!!

dots.alive = randi(dots.maxAlive,1,dots.n); % set to random up to 200 ms
% don't understand what this does
% maxAlive = 60/4 = 15 (when 60Hz is frame rate). don't know why divided by
% 4. randomly from 1 to 15 to create a n-by-n matrix where n is the number
% of dots.

dots.xdisp = dots.x;
dots.ydisp = dots.y;

% set incoherent dots to 0
dots.coherency = 1;
dots.incoherent = rand(1,dots.n) > dots.coherency;
dots.incoherentn = sum(dots.incoherent);
dots.coherent = ~dots.incoherent;

% track time
dots.time = mglGetSecs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for horizontal motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDots(dots,coherence,repick)
%%
% update dots location (frame to frame)
%
elapsed = mglGetSecs-dots.time;
dots.time = mglGetSecs;

dots.alive = dots.alive+1;
rIdx = dots.alive>dots.maxAlive;
replace = sum(rIdx);
dots.x(rIdx) = rand(1,replace)*dots.maxX;
dots.y(rIdx) = rand(1,replace)*dots.maxY;
dots.alive(rIdx) = 0;

% get the coherent and incoherent dots
if repick
    dots.incoherent = rand(1,dots.n) > coherence;
    dots.incoherentn = sum(dots.incoherent);
    dots.coherent = ~dots.incoherent;
    dots.coherency = coherence;
elseif dots.coherency ~= coherence
    cohDiff = coherence - dots.coherency;
    numDots = round(abs(cohDiff) * dots.n); % actual number of dots to flip
    if numDots > dots.n, numDots = dots.n; end
    if cohDiff > 0
        % we need to add more coherent dots
        flipDots = [zeros(1,numDots) ones(1,sum(dots.incoherent)-numDots)];
        dots.incoherent(dots.incoherent) = flipDots(randperm(length(flipDots)));
    else
        % we need to add more incoherent dots
        flipDots = [ones(1,numDots) zeros(1,sum(dots.coherent)-numDots)];
        dots.incoherent(dots.coherent) = flipDots(randperm(length(flipDots)));
    end
    dots.incoherentn = sum(dots.incoherent);
    dots.coherent = ~dots.incoherent;
    dots.coherency = sum(dots.coherent)/dots.n;
end
dots.coherentn = dots.n-dots.incoherentn;

vectorLength = dots.speed*elapsed;

% move coherent dots
dots.x(dots.coherent) = dots.x(dots.coherent) + vectorLength * cos(dots.dir);
dots.y(dots.coherent) = dots.y(dots.coherent) + vectorLength * sin(dots.dir);

dots.x(dots.incoherent) = dots.x(dots.incoherent) + vectorLength * cos(rand(1,dots.incoherentn)*2*pi);
dots.y(dots.incoherent) = dots.y(dots.incoherent) + vectorLength * sin(rand(1,dots.incoherentn)*2*pi);

offscreen = dots.x > dots.maxX;
dots.x(offscreen) = mod(dots.x(offscreen),dots.maxX);
offscreen = dots.x < 0;
dots.x(offscreen) = mod(dots.x(offscreen),dots.maxX);

offscreen = dots.y > dots.maxY;
dots.y(offscreen) = mod(dots.y(offscreen),dots.maxY);
offscreen = dots.y < 0;
dots.y(offscreen) = mod(dots.y(offscreen),dots.maxY);








