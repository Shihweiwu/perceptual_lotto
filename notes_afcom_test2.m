function [ myscreen ] = notes_afcom_test2
global stimulus fixStimulus
stimulus = struct; 
fixStimulus = struct;
stimulus.rotSpd = 90;
scan = 0;
plots = 0;
noeye = 1;
debug = 0;
replay = 0;
powerwheel = 1;
run = 0; 
eyewindow=3; 
mouse=0; 
practice=0; 
practiceType=-1;
cue=1;

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
stimulus.replay = replay; %replay? not sure what it is for. default is no.
stimulus.overrideRun = run;

clear localizer invisible scan noeye task test2 attend build eyewindow mouse practice powerwheel cue

if stimulus.scan
    warning('Disabling eyewindow');
    stimulus.eyewindow=0;
end
if ~stimulus.replay %replay default is 0. so by default we could go into this if statement. 
    stimulus.counter = 1;
    
    if ~isempty(mglGetSID) && isdir(sprintf('~/data/afcom/%s',mglGetSID))
        
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
stimulus.counter = -1;
if ~stimulus.replay
    disp('*************************');
    disp(sprintf('(afcom) This is scan #%i',stimulus.counter));
    disp('*************************');
end

if stimulus.replay
    myscreen = initScreen('replayScreen');
else
    myscreen = initScreen('VPixx');
end
% set background to black
myscreen.background = 0;

if ~stimulus.replay
    myscreen.stimulusNames{1} = 'stimulus';
    
    if stimulus.powerwheel
        stimulus.responseKeys = 1;
    else
        stimulus.responseKeys = [1 2]; % left right
    end
else
end
if isfield(myscreen,'calibFullFilename')
    calib = load(fullfile(myscreen.calibFullFilename));
    stimulus.calib = calib.calib;
else
    stimulus.calib = []; % need this so that mglLab2rgb doesn't fail
end

if ~isfield(stimulus,'colors')
    stimulus.colors.white = [0.8 0.8 0.8]; stimulus.colors.red = [0.8 0 0];
    stimulus.colors.green = [0 1 0]; stimulus.colors.black = [0 0 0];
end
stimulus.theta_ = pi/64; % increment size
stimulus.thetas = 0:stimulus.theta_:(2*pi);

%%%
stimulus.fixWidth = 0.5; %fixation width?
stimulus.targetWidth = 10; %will be used to define area to show moving dots
%%%
stimulus.patchEcc = 8; %eccentricity?
stimulus.dotScale = 0.3;
stimulus.cueScale = 0.1;
stimulus.dotDirs = [0.5 1.5 0.5 1.5]*pi;
stimulus.dotColors = [0.5 1.5 0.5 1.5]*pi;

dots = struct;
dots.density = 0.1; %some scalar value for dots density. not sure how it is defined.
dots.speed = 10; %some scalar value for dots speed. not sure how it is defined
dots.maxAlive = 1000;%myscreen.framesPerSecond/4; %don't know what it is used for.
dots.maxX = stimulus.targetWidth; %used to define x,y position of dots
dots.maxY = stimulus.targetWidth;
stimulus.dotThetas = [0 0 pi pi];
for di = 1:1
    stimulus.patches{di} = struct;
    stimulus.patches{di}.dots = initDots(dots);
    if stimulus.cue==1
        stimulus.patches{di}.dots.dir = stimulus.dotDirs(di); %dots direction is either 0.5*pi (90-deg) or 1.5*pi (270-deg)
    else
        stimulus.patches{di}.dots.dir = 0;
    end
    
    % location
    stimulus.patches{di}.theta = stimulus.dotThetas(di);
    stimulus.patches{di}.ecc = stimulus.patchEcc;
    stimulus.patches{di}.xcenter = 0;%stimulus.patches{di}.ecc * cos(stimulus.patches{di}.theta); %patch center x 
    stimulus.patches{di}.ycenter = 0;%stimulus.patches{di}.ecc * sin(stimulus.patches{di}.theta); %patch center y
end
mglClearScreen(0); %clear buffer
mglStencilCreateBegin(1);
for di = 1
    %draw a filled oval to stencil that has size targetWidth
    mglFillOval(stimulus.patches{di}.xcenter,stimulus.patches{di}.ycenter,[stimulus.targetWidth, stimulus.targetWidth],[1 1 1]);
end
task{1}{1} = struct;

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
task{1}{1}.getResponse(stimulus.seg.resp) = 1; %a vector contains 1s and 0s used as index for whether to detect a response
if stimulus.scan==1
    task{1}{1}.numTrials = 1; %why set number of trials to Inf if scan=1?
else
    task{1}{1}.numTrials = 1;
end

task{1}{1}.random = 1; %what is this for?
task{1}{1}.parameter.trialType = [1 1 1 2 2 2 0 0 3 4]; 
task{1}{1}.parameter.target = [1 1 1 1]; % which patch is the target
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
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,@endTrialCallback,[]);
end
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
    myscreen = tickScreen(myscreen,task);
end
mglClearScreen;
mglTextSet([],32,stimulus.colors.white);
% get count
mglTextDraw('Please wait',[0 0]);
mglFlush %swap front and back buffer (waits for one frame tick)
myscreen.flushMode = 1;
myscreen = endTask(myscreen,task);

function [task, myscreen] = startTrialCallback(task,myscreen)
global stimulus

stimulus.final_location_theta = (randi([0,360]));
pd=makedist('Normal','mu',stimulus.final_location_theta,'sigma',30);
t=truncate(pd,0,360);
r=random(t,1,length(stimulus.patches{1}.dots.x));
stimulus.final_location_x=(stimulus.patches{1}.dots.maxX/2)+(stimulus.patches{1}.dots.maxX/2)*cosd(r);
stimulus.final_location_y=(stimulus.patches{1}.dots.maxY/2)+(stimulus.patches{1}.dots.maxY/2)*sind(r);
% swap seglen in
task.thistrial.seglen(stimulus.seg.stim) = task.thistrial.duration;

% get the current mouse position:
mInfo = mglGetMouse(myscreen.screenNumber);
stimulus.live.mouseStart = -mInfo.x/stimulus.rotSpd;

for di = 1:length(stimulus.patches)
    task.thistrial.(sprintf('angle%i',di)) = stimulus.final_location_theta;
end

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
function [task, myscreen] = endTrialCallback(task,myscreen)
%%
if task.thistrial.dead, return; end
% if this is true, then exit function.

respType = {'click','click','multiclick','multiclick','multiclick'};
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

function drawStim(task,stimSeg)
%%
global stimulus

mglStencilSelect(1); %before you draw, call this to allow drawing on stencils set earlier.
% update and collapse x/y coordinates for drawing
n = stimulus.patches{1}.dots.n;

x = nan(1,n*length(stimulus.patches));
y = x;
r = stimulus.colors.white(1)*ones(1,n*length(stimulus.patches));
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
            r(((di-1)*n+1):(di*n)) = stimulus.colors.white(1);
            g(((di-1)*n+1):(di*n)) = stimulus.colors.white(2);
            b(((di-1)*n+1):(di*n)) = stimulus.colors.white(3);
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

function drawResp(angle)
%%
global stimulus
% Draw the chosen color as a backgroundc ircle
if stimulus.cue==1
    mglFillOval(0,0,stimulus.fixWidth*[1 1],[0.5,0.5,0.5]);
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
    mInfo = mglGetMouse;
    a=(mInfo.x-myscreen.screenWidth/2)*myscreen.imageWidth/myscreen.screenWidth;
    b=(mInfo.y-myscreen.screenHeight/2)*myscreen.imageHeight/myscreen.screenHeight;
    mglFillOval(a,b,[2,2],[0.5 0.5 0.5]);

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
        %drawTarget(task);
        %drawPicker(task);
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
function dots = initDots(dots)

%%
dots.dir = 0;
area = (dots.maxX)*(dots.maxY)*pi;
dots.n = round(area * dots.density);
dots.n = dots.n + mod(dots.n,2);
theta=randi([0,360],1,dots.n);


dots.x = (dots.maxX/2)+rand(1,dots.n)*(dots.maxX/2).*cosd(theta);
dots.y = (dots.maxY/2)+rand(1,dots.n)*(dots.maxY/2).*sind(theta);

dots.alive = dots.maxAlive;%randi(dots.maxAlive,1,dots.n); % set to random up to 200 ms

dots.xdisp = dots.x;
dots.ydisp = dots.y;

% track time
dots.time = mglGetSecs;
function dots = updateDots(dots,coherence,repick)
%%
% update dots location (frame to frame)
%
global stimulus
elapsed = mglGetSecs-dots.time;
dots.time = mglGetSecs;

dots.alive = dots.alive+1;

vectorLength = dots.speed*elapsed;
angle_togo=atan2(stimulus.final_location_y-dots.y,stimulus.final_location_x-dots.x);
angle_trans=wrapTo360(rad2deg(angle_togo));

dots.x = dots.x + vectorLength * cosd(angle_trans);
dots.y = dots.y + vectorLength * sind(angle_trans);

dis = sqrt(((dots.maxX/2)-dots.x).^2+((dots.maxY/2)-dots.y).^2);

if ~isempty(find(dis>=(dots.maxX/2)))       
 offscreen = (dis>=(dots.maxX/2));
 theta=randi([0,360],1,length(find(offscreen==1)));
 dots.x(offscreen) = (dots.maxX/2)+rand(1,length(find(offscreen==1)))*(dots.maxX/2).*cosd(theta);
 dots.y(offscreen) = (dots.maxY/2)+rand(1,length(find(offscreen==1)))*(dots.maxY/2).*sind(theta);
end

