function myscreen = testExperiment_dotCatch_wResp(screenName)
%
% Modified from testExperiment.
% 2019.05.04: SWW. Change the task such that
% 1. There is no staircase task (previously task{1}).
% 2. Task 2 has two phases -- randomly moving dots (phase 1) followed by a fixation (phase 2).
% Previously phase 2 is randomly moving dots with 0 coherence.
%


% check arguments
if ~any(nargin == [0 1])
    help testExperiment
    return
end

% initalize the screen
if nargin >= 1
    myscreen = initScreen(screenName);
else
    myscreen = initScreen;
end
myscreen.saveData = 1;
myscreen.datadir = cd;
myscreen.subjectID = 'sub01';

%Hide cursor
mglDisplayCursor(0);

% set the first task to be fixation task
global fixStimulus
fixStimulus.pospix = [myscreen.screenWidth myscreen.screenHeight]/2;
[task{1},myscreen] = fixation(myscreen);

% the second task will have 5 segments -- (1)fixation,(2)moving
% dots, (3)response, (4)confirmation, (5)feedback 
task{2}{1}.segmin = [1 3 5 1 1];
task{2}{1}.segmax = [1 3 5 1 1];
% These two parameters, mu and keppa, are the two parameters of the von
% Mises distribution
task{2}{1}.parameter.mu = (0:30:360)*(pi/180); %mean parameter in radian
task{2}{1}.parameter.keppa = [1 10 50 100 200]; %dispersion parameter. greater value means smaller SD.
task{2}{1}.random = 1; %randomize the order in which you present different combinations

%Set initial position of the mouse randomly
% task{2}{1}.randVars.initAngledeg=randi([0,359],[numel(task{2}{1}.randVars.myRandomDir),1]); % from steeve's code taskDotDir.
task{2}{1}.randVars.initAngledeg=0:1:359;

%Get response and reaction time
%response segment is set at "3" or "4" for mouse event
task{2}{1}.getResponse=[0 0 1 0 0];
task{2}{1}.numTrials=20; %for now.
%task{2}{1}.numTrials=numel(task{2}{1}.randVars.myRandomDir); %from Steeve's code taskDotDir.

%Store data
task{2}{1}.randVars.calculated.prodcoor=[nan nan];

% initialize our task
for phaseNum = 1:length(task{2})
    [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback,@getResponseCallBack);
end

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen); %this is to have the stimulus variable saved at the end of experiment
stimulus = initDots(stimulus,myscreen); %subfunction see below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task{2})) && ~myscreen.userHitEsc
    % update the dots
    [task{2} myscreen phaseNum] = updateTask(task{2},myscreen,phaseNum);
    
    %update the fixation task
    [task{1},myscreen]=updateTask(task{1},myscreen,1);
    
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);
clear global stimulus;
clear global fixStimulus;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;
global fixStimulus;
global Inicoord;

if (task.thistrial.thisseg == 2)
    % get (mu,keppa) of the von Mises of this trial
    stimulus.dots.mu = task.thistrial.mu; % mean parameter of von Mises
    %mu=stimulus.dots.mu;
    stimulus.dots.keppa = task.thistrial.keppa; %dispersion parameter of von Mises
    
    % get initial position
    stimulus.dots.x = stimulus.dots.width*(rand(1,stimulus.dots.n) - 0.5);
    stimulus.dots.y = stimulus.dots.height*(rand(1,stimulus.dots.n) - 0.5);
    
    % get pass-through location
    % sample from von Mises
    theta = slCirc_vmrnd2(stimulus.dots.mu,stimulus.dots.keppa,1,stimulus.dots.n); %sample from von Mises
    % figure(1);hold on;
    % hist(theta)
    
    r=0.5;
    % calculate pass-through location
    x_passthru = stimulus.dots.xcenter + r*stimulus.dots.rmax*cos(theta);
    y_passthru = stimulus.dots.ycenter + r*stimulus.dots.rmax*sin(theta);
    stimulus.dots.x_passthru = x_passthru;
    stimulus.dots.y_passthru = y_passthru;
    
    % size(x_passthru)
    % size(stimulus.dots.x)
    
    % calculate motion direction
    d_x = x_passthru - stimulus.dots.x;
    d_y = y_passthru - stimulus.dots.y;
    dir_motion = atan(d_y./d_x);
    index_x_st0 = d_x<0; index_x_gt0 = d_x>0;
    index_y_st0 = d_y<0; index_y_gt0 = d_y>0;
    index_Q2 = index_x_st0 + index_y_gt0 ==2;
    index_Q3 = index_x_st0 + index_y_st0 ==2;
    index_Q4 = index_x_gt0 + index_y_st0 ==2;
    dir_motion(index_Q2) = dir_motion(index_Q2)+pi;
    dir_motion(index_Q3) = dir_motion(index_Q3)+pi;
    dir_motion(index_Q4) = dir_motion(index_Q4)+2*pi;
    
    % direction of motion; each dot has its own direction
    stimulus.dots.dir = dir_motion;
    
    clear theta x_pass y_pass d_x d_y d dir_motion index_Q2 index_Q3 index_Q4 index_x_st0 index_x_gt0 index_y_st0 index_y_gt0
end

%response (from Steeve's taskDotDir)
if (task.thistrial.thisseg == 3)

    %Position mouse in the center of the screen
    %mglSetMousePosition(fixStimulus.pospix(1), fixStimulus.pospix(2), myscreen.screenNumber);
    %Calculate angle coordinates on aperture relative to center
    [coord]=polar2cartesian((225+270)/2,stimulus.dots.rmax);
    Inicoord.x.angle.onap=coord.x;
    Inicoord.y.angle.onap=coord.y;

    %Calculate pixels coordinates on aperture relative to center
    Inicoord.x.pix.onap=Inicoord.x.angle.onap * mglGetParam('xDeviceToPixels');
    Inicoord.y.pix.onap=Inicoord.y.angle.onap * mglGetParam('yDeviceToPixels');

    %Calculate coordinates on aperture relative to screen's root.
    Inicoord2root.x.pix.onap=Inicoord.x.pix.onap + fixStimulus.pospix(1);
    Inicoord2root.y.pix.onap=Inicoord.y.pix.onap + fixStimulus.pospix(2);

    %Position mouse
    mglSetMousePosition(Inicoord2root.x.pix.onap, Inicoord2root.y.pix.onap, myscreen.screenNumber);
end

% confirmation (from Steeve's taskDotDir)
if (task.thistrial.thisseg == 4)

    %Get mouse position
    mi=mglGetMouse(myscreen.screenNumber); %get pixel positions.
    mouseinfo.x.pix=mi.x;
    mouseinfo.y.pix=mi.y;

    %Check if subject confirmed his choice: ("space bar is down")
    %if you want to you use mouse click to enter choice
    %mouseinfo.buttons=mi.buttons;

    %if you want to you use keyboard button '1' to enter choice
    mouseinfo.buttons=mglGetKeys(19);

    %Position mouse on aperture
    [mouseinfo]=pos2ap(myscreen,fixStimulus,mouseinfo,stimulus);

    %Interface PowerMate
    [mouseinfo]=PowerMate2ap(mouseinfo,stimulus,task);

    %Collect response
    task.thistrial.prodcoor=[mouseinfo.x.angle.onap mouseinfo.y.angle.onap];%(visual angle)
    [~,deg] = SLcart2polar(task.thistrial.prodcoor);

    %print trial number with the key and direction choosen.
%     fprintf('%i %1.0f %1.0f \n',...
%         task.trialnum,task.thistrial.myRandomDir,deg)
    %then jump to feedback
    task=jumpSegment(task,4);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)

global stimulus
global fixStimulus
mglClearScreen;

% segment 2 is to display moving dots
if (task.thistrial.thisseg == 2)
    stimulus = updateDots(stimulus,myscreen);
end
% else
%     mglFixationCross(1,1,[1 1 1]);
% end

% segment 5 is to display a single moving dot and see if the basket
% subjects place catches it 
if (task.thistrial.thisseg == 5)
    stimulus = shootDots(stimulus,myscreen);
end

% segment 3: response (from Steeve's taskDotDir); subjects can move the
% mouse freely
if (task.thistrial.thisseg == 3)

    % Position mouse on aperture
    mi=mglGetMouse(myscreen.screenNumber); %(pixels)
    mouseinfo.x.pix=mi.x;
    mouseinfo.y.pix=mi.y;
    mouseinfo=pos2ap(myscreen,fixStimulus,mouseinfo,stimulus);

    % Interface PowerMate
    [mouseinfo]=PowerMate2ap(mouseinfo,stimulus,task);

    % draw a white arrow (radius)
    mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
        mouseinfo.x.angle.onap,...
        mouseinfo.y.angle.onap,2,[1 1 1],1);
end

% segment 4: confirmation (from Steeve's taskDotDir); subjects need to
% confirm his choice on basket placement
if (task.thistrial.thisseg == 4)

    %Position mouse on aperture
    %deliver pixel positions, not angle
    mi=mglGetMouse(myscreen.screenNumber);
    mouseinfo.x.pix=mi.x;
    mouseinfo.y.pix=mi.y;
    mouseinfo=pos2ap(myscreen,fixStimulus,mouseinfo,stimulus);

    %Interface PowerMate
    [mouseinfo]=PowerMate2ap(mouseinfo,stimulus,task);

    %Confirm subjects' choice by drawing a red arrow
    mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
        mouseinfo.x.angle.onap,...
        mouseinfo.y.angle.onap,2,[1 0 0],1);
end

% get response (from Steeve's taskDotDir)
function [task,myscreen]=getResponseCallBack(task,myscreen)

%When subject chose, confirm (jump to segment 4 in updateScreenCallback)
task=jumpSegment(task,4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen)

% convert the passed in parameters to real units
if ~isfield(stimulus,'dots') || ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax = 15;,end
if ~isfield(stimulus.dots,'xcenter'), stimulus.dots.xcenter = 0;,end
if ~isfield(stimulus.dots,'ycenter'), stimulus.dots.ycenter = 0;,end
if ~isfield(stimulus.dots,'dotsize'), stimulus.dots.dotsize = 3;,end
if ~isfield(stimulus.dots,'density'), stimulus.dots.density = 3;,end
if ~isfield(stimulus.dots,'coherence'), stimulus.dots.coherence = 1;,end
if ~isfield(stimulus.dots,'speed'), stimulus.dots.speed = 5;,end
if ~isfield(stimulus.dots,'dir'), stimulus.dots.dir = 0;,end

% actually a square patch of dots that get stenciled
% so calculate width and height
stimulus.dots.width = stimulus.dots.rmax*2;
stimulus.dots.height = stimulus.dots.rmax*2;

% get the number of dots
stimulus.dots.n = round(stimulus.dots.width*stimulus.dots.height*stimulus.dots.density);

% get max and min points for dots
stimulus.dots.xmin = -stimulus.dots.width/2;
stimulus.dots.xmax = stimulus.dots.width/2;
stimulus.dots.ymin = -stimulus.dots.height/2;
stimulus.dots.ymax = stimulus.dots.height/2;

% get the step size
stimulus.dots.stepsize = stimulus.dots.speed/myscreen.framesPerSecond;

% create stencil
mglClearScreen;
mglStencilCreateBegin(1);
% and draw that oval
mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax]/2,[1 1 1],60); %original code
mglStencilCreateEnd;
mglClearScreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update dot positions and draw them to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateDots(stimulus,myscreen)

% When dots go off the patch, place new dots on the patch.
bin_xOut = zeros(1,stimulus.dots.n);
bin_yOut = zeros(1,stimulus.dots.n);
r=0.5;
% Identify which dots are out in the x direction
if sum(stimulus.dots.x < r*stimulus.dots.xmin)>0 || sum(stimulus.dots.x > r*stimulus.dots.xmax)>0
    index_xOut = [find(stimulus.dots.x < r*stimulus.dots.xmin) find(stimulus.dots.x > r*stimulus.dots.xmax)];
    bin_xOut(index_xOut)=1;
end
% Identify which dots are out in the y direction
if sum(stimulus.dots.y < r*stimulus.dots.ymin)>0 || sum(stimulus.dots.y > r*stimulus.dots.ymax)>0
    index_yOut = [find(stimulus.dots.y < r*stimulus.dots.ymin) find(stimulus.dots.y > r*stimulus.dots.ymax)];
    bin_yOut(index_yOut)=1;
end
% either out in x or y are out
% test=(bin_xOut + bin_yOut);
% figure(1);hold on;
% hist(test);
index_out = find((bin_xOut + bin_yOut) >0);
n_out = length(index_out);

if n_out > 0
    % For dots that go off the patch, place new dots at random location to replace them
    stimulus.dots.x(index_out) = stimulus.dots.xcenter + stimulus.dots.width*(rand(1,n_out)-0.5);
    stimulus.dots.y(index_out) = stimulus.dots.ycenter + stimulus.dots.height*(rand(1,n_out)-0.5);
    
    % For new dots, sample from von Mises to determine pass-through
    % location
    theta_new = slCirc_vmrnd2(stimulus.dots.mu,stimulus.dots.keppa,1,n_out); %sample from von Mises
    
    % calculate pass-through location
    x_passthru = stimulus.dots.xcenter + r*stimulus.dots.rmax*cos(theta_new);
    y_passthru = stimulus.dots.ycenter + r*stimulus.dots.rmax*sin(theta_new);
    stimulus.dots.x_passthru(index_out) = x_passthru;
    stimulus.dots.y_passthru(index_out) = y_passthru;
    
    % calculate motion direction
    d_x = x_passthru - stimulus.dots.x(index_out);
    d_y = y_passthru - stimulus.dots.y(index_out);
    dir_new = atan(d_y./d_x);
    
    index_x_st0 = d_x<0; index_x_gt0 = d_x>0;
    index_y_st0 = d_y<0; index_y_gt0 = d_y>0;
    index_Q2 = index_x_st0 + index_y_gt0==2; 
    index_Q3 = index_x_st0 + index_y_st0==2; 
    index_Q4 = index_x_gt0 + index_y_st0==2; 
    dir_new(index_Q2) = dir_new(index_Q2)+pi;
    dir_new(index_Q3) = dir_new(index_Q3)+pi;
    dir_new(index_Q4) = dir_new(index_Q4)+2*pi;
    
    % update direction of motion for new dots
    stimulus.dots.dir(index_out) = dir_new;
    
    clear theta_new x_pass y_pass d_x d_y d dir_new index_Q2 index_Q3 index_Q4 index_x_st0 index_x_gt0 index_y_st0 index_y_gt0
end

% get the dots step: each dot has its own (xstep,ystep)
stimulus.dots.xstep = cos(stimulus.dots.dir).*stimulus.dots.stepsize;
stimulus.dots.ystep = sin(stimulus.dots.dir).*stimulus.dots.stepsize;

% move dots
stimulus.dots.x = stimulus.dots.x+stimulus.dots.xstep;
stimulus.dots.y = stimulus.dots.y+stimulus.dots.ystep;

% % randomwalk rule
% thisdir = rand(1,sum(~stimulus.dots.coherent))*2*pi; % the direction of non-coherently moving dots
% stimulus.dots.x(~stimulus.dots.coherent) = stimulus.dots.x(~stimulus.dots.coherent)+cos(thisdir)*stimulus.dots.stepsize;
% stimulus.dots.y(~stimulus.dots.coherent) = stimulus.dots.y(~stimulus.dots.coherent)+sin(thisdir)*stimulus.dots.stepsize;

% movshon noise
%stimulus.dots.x(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.width;
%stimulus.dots.y(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.height;

% draw the dots
mglStencilSelect(1);
mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize,[1 1 1]);
mglStencilSelect(0);

function stimulus = shootDots(stimulus,myscreen)
% sample one dot from the designated von Mises and let it move and see if the
% basket catches it.
%

% pick one dot and put it at random location
stimulus.dots.xshoot = stimulus.dots.xcenter + stimulus.dots.width*(rand-0.5);
stimulus.dots.yshoot = stimulus.dots.ycenter + stimulus.dots.height*(rand-0.5);

% Sample from von Mises to determine its pass-through location
theta = slCirc_vmrnd2(stimulus.dots.mu,stimulus.dots.keppa,1,1); %sample from von Mises

% calculate pass-through location
r=0.5;
x_passthru = stimulus.dots.xcenter + r*stimulus.dots.rmax*cos(theta);
y_passthru = stimulus.dots.ycenter + r*stimulus.dots.rmax*sin(theta);
stimulus.dots.xshoot_passthru = x_passthru;
stimulus.dots.yshoot_passthru = y_passthru;

% calculate motion direction
d_x = x_passthru - stimulus.dots.xshoot;
d_y = y_passthru - stimulus.dots.yshoot;
dir_shoot = atan(d_y./d_x);

if d_x<0 && d_y>0 %quadrant 2
    dir_shoot = dir_shoot+pi;
elseif d_x<0 && d_y<0 %quadrant 3
    dir_shoot = dir_shoot+pi;
elseif d_x>0 && d_y<0 %quadrant 4
    dir_shoot = dir_shoot+2*pi;
end

% update direction of motion for new dots
stimulus.dots.dir_shoot = dir_shoot;

function [mouseinfo]=pos2ap(myscreen,fixStimulus,mouseinfo,stimulus)
%function to convert from mouse's to aperture's coordinates
% From Steeve's taskDotDir

%and the mouse coordinates relative to the center of the target screen
mouseinfo.x.pix=mouseinfo.x.pix-fixStimulus.pospix(1); %(in pix)
mouseinfo.y.pix=mouseinfo.y.pix-fixStimulus.pospix(2);

%convert pixels 2 angles (mglLines works with angles)
mouseinfo.x.angle.screen=mouseinfo.x.pix*mglGetParam('xPixelsToDevice');
mouseinfo.y.angle.screen=mouseinfo.y.pix*mglGetParam('yPixelsToDevice');

%calculate the transformation parameter that position the mouse on the aperture
%calculate length of the arrow (pythagoras theorem)
arrowsz2=mouseinfo.x.angle.screen^2 + mouseinfo.y.angle.screen^2;%(angle)

%calculate transformation parameter
transpara=stimulus.dots.rmax^2/arrowsz2;

%transform actual coordinates to on-aperture coordinates.
mouseinfo.x.angle.onap=sqrt(mouseinfo.x.angle.screen^2*transpara)*sign(mouseinfo.x.angle.screen);
mouseinfo.y.angle.onap=sqrt(mouseinfo.y.angle.screen^2*transpara)*sign(mouseinfo.y.angle.screen);

function [mouseinfo]=PowerMate2ap(mouseinfo,stimulus,task)
%function to convert from PowerMate's to aperture's coordinates
% From Steeve's taskDotDir

global Inicoord;
%Set speed, initial position and size of arrow.
wheelspeed=0.008;
%fixed
%initAngle=(225+270)/2*pi/180;
%random
initAngle.rad=task.thistrial.initAngledeg*pi/180;

r=stimulus.dots.rmax;

%Calculate arrow's coordinates on aperture
mouseinfo.x.angle.onap= r * cos(initAngle.rad + wheelspeed*(mouseinfo.x.pix - Inicoord.x.pix.onap)); %(-) --> goes left when turn left, vvs.
mouseinfo.y.angle.onap= r * sin(initAngle.rad + wheelspeed*(mouseinfo.x.pix - Inicoord.x.pix.onap));

function [coord]=polar2cartesian(theta,r)
%function to convert from polar to cartesian coordinates
% From Steeve's taskDotDir

%theta is an angle in degree
%r is the radius of the unit circle
%Coord are in visual angle
%Record angle in degree
theta2.deg=theta;
%Convert from degree to radian
theta2.rad=theta2.deg*pi/180;
%Calculate visual angles coordinates
coord.x=r*cos(theta2.rad);
coord.y=r*sin(theta2.rad);



