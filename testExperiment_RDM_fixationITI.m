function myscreen = testExperiment_RDM_fixationITI(screenName)
%
% Modified from testExperiment.
% 2019.05.04: SWW. Change the task such that
% 1. There is no staircase task (previously task{1}).
% 2. Task 2 has two phases -- randomly moving dots (phase 1) followed by a fixation (phase 2).
% Previously phase 2 is randomly moving dots with 0 coherence.
%
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

% set the first task to be the fixation staircase task
% [task{1} myscreen] = fixStairInitTask(myscreen);

% set our task to have two phases. 
% one starts out with dots moving for incohrently for 10 seconds
task{1}{1}.waitForBacktick = 1;
task{1}{1}.seglen = 10;
task{1}{1}.numBlocks = 1;
task{1}{1}.parameter.dir = 0;
task{1}{1}.parameter.coherence = 0;

% the second phase has 2 second bursts of directions, followed by 
% a top-up period of the same direction
task{1}{2}.segmin = [2 6];
task{1}{2}.segmax = [2 10];
task{1}{2}.parameter.dir = 0:60:360; % initTask will define task.thistrial.dir, which will be used to define stimulus.dots.dir in startSegmentCallback 
task{1}{2}.parameter.coherence = 0:0.2:1; %initTask will define task.thistrial.dir, which will be used to define stimulus.dots.coherence in startSegmentCallback
task{1}{2}.random = 1; %randomize the order in which you present different directions

% initialize our task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback);
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
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  % update the dots
  [task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);

  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);
clear global stimulus;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;
if (task.thistrial.thisseg == 1)
  stimulus.dots.coherence = task.thistrial.coherence;
else
  stimulus.dots.coherence = 0.5;
end
stimulus.dots.dir = task.thistrial.dir; %direction of coherently moving dots.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)

global stimulus
mglClearScreen;
if (task.thistrial.thisseg == 1)
    stimulus = updateDots(stimulus,myscreen);
else
    mglFixationCross(1,1,[1 1 1]);
end

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

% get initial position
stimulus.dots.x = rand(1,stimulus.dots.n)*stimulus.dots.width;
stimulus.dots.y = rand(1,stimulus.dots.n)*stimulus.dots.height;

% get the step size
stimulus.dots.stepsize = stimulus.dots.speed/myscreen.framesPerSecond;

% create stencil
mglClearScreen;
mglStencilCreateBegin(1);
% and draw that oval
mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax]/2,[1 1 1],60);
mglStencilCreateEnd;
mglClearScreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update dot positions and draw them to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateDots(stimulus,myscreen)

% get the dots step
stimulus.dots.xstep = cos(pi*stimulus.dots.dir/180)*stimulus.dots.stepsize;
stimulus.dots.ystep = sin(pi*stimulus.dots.dir/180)*stimulus.dots.stepsize;

% pick a random set of dots
stimulus.dots.coherent = rand(1,stimulus.dots.n) < stimulus.dots.coherence;

% now move those dots in the right direction
stimulus.dots.x(stimulus.dots.coherent) = stimulus.dots.x(stimulus.dots.coherent)+stimulus.dots.xstep;
stimulus.dots.y(stimulus.dots.coherent) = stimulus.dots.y(stimulus.dots.coherent)+stimulus.dots.ystep;

% randomwalk rule
thisdir = rand(1,sum(~stimulus.dots.coherent))*2*pi; % the direction of non-coherently moving dots
stimulus.dots.x(~stimulus.dots.coherent) = stimulus.dots.x(~stimulus.dots.coherent)+cos(thisdir)*stimulus.dots.stepsize;
stimulus.dots.y(~stimulus.dots.coherent) = stimulus.dots.y(~stimulus.dots.coherent)+sin(thisdir)*stimulus.dots.stepsize;

% movshon noise
%stimulus.dots.x(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.width;
%stimulus.dots.y(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.height;

% make sure we haven't gone off the patch
% do the dots separately for left and right hand side
stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin) = stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin)+stimulus.dots.width;
stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax) = stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax)-stimulus.dots.width;
stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin) = stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin)+stimulus.dots.height;
stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax) = stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax)-stimulus.dots.height;

% draw the dots
mglStencilSelect(1);
mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize,[1 1 1]);
mglStencilSelect(0);




