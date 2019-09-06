function myscreen = testExperiment_dotCatch_fromCtr(screenName)
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

% set the first task to be fixation task
global fixStimulus
fixStimulus.pospix = [myscreen.screenWidth myscreen.screenHeight]/2;
[task{1},myscreen] = fixation(myscreen);

% the second task will have two segments -- moving dots and fixation
task{2}{1}.segmin = [2 1];
task{2}{1}.segmax = [2 1];
% These two parameters, mu and keppa, are the two parameters of the von
% Mises distribution
task{2}{1}.parameter.mu = (0:30:360)*(pi/180); %mean parameter in radian
task{2}{1}.parameter.keppa = [1 10 50 100 200]; %dispersion parameter. greater value means smaller SD.
task{2}{1}.random = 1; %randomize the order in which you present different combinations

% initialize our task
for phaseNum = 1:length(task{2})
    [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

% get (mu,keppa) of the von Mises of this trial
stimulus.dots.mu = task.thistrial.mu; % mean parameter of von Mises
%mu=stimulus.dots.mu;
stimulus.dots.keppa = task.thistrial.keppa; %dispersion parameter of von Mises

% get initial position
stimulus.dots.x = stimulus.dots.xcenter*ones(1,stimulus.dots.n)+40*stimulus.dots.stepsize*(rand(1,stimulus.dots.n)-0.5);
stimulus.dots.y = stimulus.dots.ycenter*ones(1,stimulus.dots.n)+40*stimulus.dots.stepsize*(rand(1,stimulus.dots.n)-0.5);

% get pass-through location
% sample from von Mises
theta = slCirc_vmrnd2(stimulus.dots.mu,stimulus.dots.keppa,1,stimulus.dots.n); %sample from von Mises
% figure(1);hold on;
% hist(theta)

% calculate pass-through location
r=1;
x_pass = r*stimulus.dots.rmax*cos(theta);
y_pass = r*stimulus.dots.rmax*sin(theta);

% size(x_pass)
% size(stimulus.dots.x)

% calculate motion direction
d_x = x_pass - stimulus.dots.x;
d_y = y_pass - stimulus.dots.y;
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
if ~isfield(stimulus.dots,'density'), stimulus.dots.density = 1;,end
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
mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax]/2,[1 1 1],60);
mglStencilCreateEnd;
mglClearScreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update dot positions and draw them to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateDots(stimulus,myscreen)

% When dots go off the patch, place new dots on the patch.
bin_xOut = zeros(1,stimulus.dots.n);
bin_yOut = zeros(1,stimulus.dots.n);
% Identify which dots are out in the x direction
r=1;
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
index_out = find((bin_xOut + bin_yOut) >0);
n_out = length(index_out);

if n_out > 0
    % For dots that go off the patch, place new dots at random location to replace them
    stimulus.dots.x(index_out) = stimulus.dots.xcenter*ones(1,n_out)+40*stimulus.dots.stepsize*(rand(1,n_out)-0.5);
    stimulus.dots.y(index_out) = stimulus.dots.ycenter*ones(1,n_out)+40*stimulus.dots.stepsize*(rand(1,n_out)-0.5);
    
    % For new dots, sample from von Mises to determine pass-through
    % location
    theta_new = slCirc_vmrnd2(stimulus.dots.mu,stimulus.dots.keppa,1,n_out); %sample from von Mises
    
    % calculate pass-through location
    x_pass = stimulus.dots.xcenter+r*stimulus.dots.rmax*cos(theta_new);
    y_pass = stimulus.dots.xcenter+r*stimulus.dots.rmax*sin(theta_new);
    
    % calculate motion direction
    d_x = x_pass - stimulus.dots.x(index_out);
    d_y = y_pass - stimulus.dots.y(index_out);
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




