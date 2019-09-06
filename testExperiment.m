% testExperiment.m
%
%        $Id$
%      usage: testExperiment
%         by: justin gardner
%       date: 09/07/06
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: example program to show how to use the task structure
%
function myscreen = testExperiment(screenName)

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
[task{1} myscreen] = fixStairInitTask(myscreen);

% set our task to have two phases. 
% one starts out with dots moving for incohrently for 10 seconds
task{2}{1}.waitForBacktick = 1;
task{2}{1}.seglen = 10;
task{2}{1}.numBlocks = 1;
task{2}{1}.parameter.dir = 0;
task{2}{1}.parameter.coherence = 0;

% the second phase has 2 second bursts of directions, followed by 
% a top-up period of the same direction
task{2}{2}.segmin = [2 6];
task{2}{2}.segmax = [2 10];
task{2}{2}.parameter.dir = 0:60:360;
task{2}{2}.parameter.coherence = 1;
task{2}{2}.random = 1;

% initialize our task
for phaseNum = 1:length(task{2})
  [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback);
end

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
stimulus = initDots(stimulus,myscreen);

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
  % update the fixation task
  [task{1} myscreen] = updateTask(task{1},myscreen,1);
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
  stimulus.dots.coherence = 0;
end
stimulus.dots.dir = task.thistrial.dir;
% coh=stimulus.dots.coherence
% dd=stimulus.dots.dir

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)

global stimulus
mglClearScreen;
stimulus = updateDots(stimulus,myscreen);

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

% get initial position: original definition
stimulus.dots.x = rand(1,stimulus.dots.n)*stimulus.dots.width;
stimulus.dots.y = rand(1,stimulus.dots.n)*stimulus.dots.height;
%
% it shouldn't matter where you set these dots to be -- as long as they are
% scattered evenly within the boundary defined by the width and height of
% the patch.
% This is because their position will be updated almost immediately in "updateDots" (see the function "updateDots" below) 
% where if the dots fall outside the boundary [xmin xmax ymin ymax] they will be redrawn at
% the opposite side of boundary where they fall off. 
% Try the following:
% stimulus.dots.x = (rand(1,stimulus.dots.n))*stimulus.dots.width+10000;
% stimulus.dots.y = (rand(1,stimulus.dots.n))*stimulus.dots.height;
% Or
% stimulus.dots.x = (rand(1,stimulus.dots.n)-0.5)*stimulus.dots.width;
% stimulus.dots.y = (rand(1,stimulus.dots.n)-0.5)*stimulus.dots.height;
% You will notice no difference in the dots. They look just like using the original definition. 
%
% However, if they are set too densely, like the following:
% stimulus.dots.x = rand(1,stimulus.dots.n)*stimulus.dots.width*0.1;
% stimulus.dots.y = rand(1,stimulus.dots.n)*stimulus.dots.height*0.1;
% there will be problems

% get the step size
stimulus.dots.stepsize = stimulus.dots.speed/myscreen.framesPerSecond;

% create stencil
mglClearScreen;
mglStencilCreateBegin(1);
% and draw that oval: original code
mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax]/2,[1 1 1],60);
% mglGluDisk(x,y,size,...) size is the radius.
% I know size here means radius because when size is greater than rmax --
% which is half of height and width (see definition above) then I start
% seeing non-circular aperture. Try the following:
% mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax
% stimulus.dots.rmax]*1,[1 1 1],60); % still look like a circular aperture
% mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax
% stimulus.dots.rmax]*1.2,[1 1 1],60); % no longer a circular aperture
% mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax
% stimulus.dots.rmax]*1.5,[1 1 1],60); % a square aperture.
%
% So the original code set the aperture to be half the size of the maximum possible size for the stimuli to look like 
% a circular aperture.
% 
% What is not clear is why set to be half the maximum possible size?
% 
% Another question is the color of the patch. I change color but it doesn't
% make any difference. So why set color? what does it do?
%
% mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax]/2,[1 0 0],60);
% mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax]/2,[0 0 0],60);

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
thisdir = rand(1,sum(~stimulus.dots.coherent))*2*pi;
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




