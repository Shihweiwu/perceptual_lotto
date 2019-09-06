function myscreen = testExperiment_dotCatch_limitLife_2p2(screenName)
%
% Modified from testExperiment_dotCatch_limitLife_2p1. dot-catching task.
% > No redraw. set stimulus.dots.nFramesToReDraw to very high
% > When a dot moves out, replace it at random location.
%
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
task{2}{1}.segmin = [10 1];
task{2}{1}.segmax = [10 1];
% These two parameters, mu and keppa, are the two parameters of the von
% Mises distribution
task{2}{1}.parameter.mu = (0:30:360)*(pi/180); %mean parameter in radian
task{2}{1}.parameter.keppa = [1 10 50 100 200]; %dispersion parameter. greater value means smaller SD.
task{2}{1}.random = 1; %randomize the order in which you present different combinations
% task{2}{1}.parameter.coherence = [0.05 0.1 0.2 0.4];
task{2}{1}.parameter.coherence = [0.2];


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

global stimulus

% get (mu,keppa) of the von Mises of this trial
stimulus.dots.mu = task.thistrial.mu; % mean parameter of von Mises
stimulus.dots.keppa = task.thistrial.keppa; %dispersion parameter of von Mises
stimulus.dots.coherence = task.thistrial.coherence; %boundary for adding random angle to the angle selected from von Mises.
stimulus.dots.frameCount = 0;
stimulus.dots.n_VM_dots = round(stimulus.dots.n*stimulus.dots.coherence);
% s=task.thistrial.keppa;
% mu=task.thistrial.mu;
% coh=task.thistrial.coherence;

% get initial position: randomly place the dots within a circular aperture.
% This will get rid of the "edgy" perception due to placing them within a sqaure aperture.
r_grid = unifrnd(0,stimulus.dots.rmax,1,stimulus.dots.n);
theta_start = unifrnd(0,2*pi,1,stimulus.dots.n);
stimulus.dots.x = stimulus.dots.xcenter + r_grid.*cos(theta_start);
stimulus.dots.y = stimulus.dots.ycenter + r_grid.*sin(theta_start);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)

global stimulus
mglClearScreen;
if (task.thistrial.thisseg == 1)
    stimulus = updateDots_2(stimulus,myscreen);
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

% % get initial position
% stimulus.dots.x = rand(1,stimulus.dots.n)*stimulus.dots.width;
% stimulus.dots.y = rand(1,stimulus.dots.n)*stimulus.dots.height;
% The above code randomly places all dots within [0 width] in x and [0
% height] in y.

% get initial position: randomly place the dots within a circular aperture.
% This will get rid of the "edgy" perception due to placing them within a sqaure aperture.
r_grid = unifrnd(0,stimulus.dots.rmax,1,stimulus.dots.n);
theta_start = unifrnd(0,2*pi,1,stimulus.dots.n);
stimulus.dots.x = stimulus.dots.xcenter + r_grid.*cos(theta_start);
stimulus.dots.y = stimulus.dots.ycenter + r_grid.*sin(theta_start);
% size(stimulus.dots.x);
% size(stimulus.dots.y);
% get the step size
stimulus.dots.stepsize = stimulus.dots.speed/myscreen.framesPerSecond;

% initialize frame count
stimulus.dots.frameCount = 0;

% define number of frames for re-drawing all dots
stimulus.dots.nFramesToReDraw = 10^10;


% create stencil
mglClearScreen;
mglStencilCreateBegin(1);
% and draw that oval
mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax],[1 1 1],60); %original code
% mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax]/2,[1 1 1],60); %original code
% mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax]*2,[1 1 1],60); %original code
mglStencilCreateEnd;
mglClearScreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update dot positions and draw them to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateDots_2(stimulus,myscreen)

% coh = stimulus.dots.coherence;
% n_VM_dots = round(stimulus.dots.n*stimulus.dots.coherence);
% stimulus.dots.n_VM_dots = n_VM_dots;

% stimulus.dots.frameCount = stimulus.dots.frameCount+1;
% fc=stimulus.dots.frameCount;
% if stimulus.dots.frameCount>stimulus.dots.nFramesToReDraw
%     stimulus.dots.frameCount=1; % reset framecount
%     %stimulus = reDrawDots(stimulus,myscreen); % resample from von Mises to determine the direction of motion for the von Mises dots.
% end
% if stimulus.dots.frameCount==1
%     % whenever framecount gets reset
%     % randomly choose a subset of dots to be the von Mises dots
%     %rand_index = randsample(stimulus.dots.n,stimulus.dots.n_VM_dots,0)';
%     rand_order = randperm(stimulus.dots.n);
%     index_cohM = rand_order(1:stimulus.dots.n_VM_dots);
%     index_rndM = rand_order(stimulus.dots.n_VM_dots+1:stimulus.dots.n);
%     stimulus.dots.index_cohM = index_cohM;
%     stimulus.dots.index_rndM = index_rndM;
%     %stimulus.dots.coherent = zeros(1,stimulus.dots.n);
%     %stimulus.dots.coherent(rand_index)=1;
%     stimulus = reDrawDots(stimulus,myscreen);
% end

% randomly choose a subset of dots to be the von Mises dots
%rand_index = randsample(stimulus.dots.n,stimulus.dots.n_VM_dots,0)';
rand_order = randperm(stimulus.dots.n);
index_cohM = rand_order(1:stimulus.dots.n_VM_dots);
index_rndM = rand_order(stimulus.dots.n_VM_dots+1:stimulus.dots.n);
stimulus.dots.index_cohM = index_cohM;
stimulus.dots.index_rndM = index_rndM;
%stimulus.dots.coherent = zeros(1,stimulus.dots.n);
%stimulus.dots.coherent(rand_index)=1;
stimulus = reDrawDots(stimulus,myscreen);

% get the dots step: each dot has its own (xstep,ystep)
stimulus.dots.xstep = stimulus.dots.stepsize*cos(stimulus.dots.dir);
stimulus.dots.ystep = stimulus.dots.stepsize*sin(stimulus.dots.dir);

%indx=stimulus.dots.coherent;
%size(stimulus.dots.coherent)
%sum(stimulus.dots.coherent);

% move the von Mises dots in their set directions
%fc=stimulus.dots.frameCount
stimulus.dots.x(stimulus.dots.index_cohM) = stimulus.dots.x(stimulus.dots.index_cohM)+stimulus.dots.xstep;
%stimulus.dots.x(stimulus.dots.coherent) = stimulus.dots.x(stimulus.dots.coherent)+stimulus.dots.xstep;
stimulus.dots.y(stimulus.dots.index_cohM) = stimulus.dots.y(stimulus.dots.index_cohM)+stimulus.dots.ystep;

% move the remaining dots in random directions
thisdir = rand(1,stimulus.dots.n-stimulus.dots.n_VM_dots)*2*pi;
size(thisdir);
size(stimulus.dots.index_rndM);
stimulus.dots.x(stimulus.dots.index_rndM) = stimulus.dots.x(stimulus.dots.index_rndM)+cos(thisdir)*stimulus.dots.stepsize;
stimulus.dots.y(stimulus.dots.index_rndM) = stimulus.dots.y(stimulus.dots.index_rndM)+sin(thisdir)*stimulus.dots.stepsize;

% Calculate distance from center. Replace dots that move out of the
% aperture by randomly placing new dots in the aperture
% index_out = sqrt(stimulus.dots.x.^2 + stimulus.dots.y.^2) > stimulus.dots.rmax/2;
index_out = sqrt(stimulus.dots.x.^2 + stimulus.dots.y.^2) > stimulus.dots.rmax;
n_out = sum(index_out);
if n_out>0
    % r_grid = rand(1,n_out)*(stimulus.dots.rmax/2);
    r_grid = rand(1,n_out)*stimulus.dots.rmax;
    theta_new = rand(1,n_out)*2*pi;
    stimulus.dots.x(index_out) = stimulus.dots.xcenter + r_grid.*cos(theta_new);
    stimulus.dots.y(index_out) = stimulus.dots.ycenter + r_grid.*sin(theta_new);
end

% % When dots go off the patch, place new dots on the patch.
% bin_xOut = zeros(1,stimulus.dots.n);
% bin_yOut = zeros(1,stimulus.dots.n);
% r=1;
% % Identify which dots are out in the x direction
% if sum(stimulus.dots.x < r*stimulus.dots.xmin)>0 || sum(stimulus.dots.x > r*stimulus.dots.xmax)>0
%     index_xOut = [find(stimulus.dots.x < r*stimulus.dots.xmin) find(stimulus.dots.x > r*stimulus.dots.xmax)];
%     bin_xOut(index_xOut)=1;
% end
% % Identify which dots are out in the y direction
% if sum(stimulus.dots.y < r*stimulus.dots.ymin)>0 || sum(stimulus.dots.y > r*stimulus.dots.ymax)>0
%     index_yOut = [find(stimulus.dots.y < r*stimulus.dots.ymin) find(stimulus.dots.y > r*stimulus.dots.ymax)];
%     bin_yOut(index_yOut)=1;
% end
% % either out in x or y are out
% test=(bin_xOut + bin_yOut);
% figure(1);hold on;
% hist(test);
% index_out = find((bin_xOut + bin_yOut) >0);
% n_out = length(index_out);
%
%
% if n_out > 0
%     % For dots that go off the patch, place new dots at random location to replace them
% %     stimulus.dots.x(index_out) = stimulus.dots.xcenter + stimulus.dots.width*(rand(1,n_out)-0.5);
% %     stimulus.dots.y(index_out) = stimulus.dots.ycenter + stimulus.dots.height*(rand(1,n_out)-0.5);
%     stimulus.dots.x(index_out) = stimulus.dots.xcenter + unifrnd(-stimulus.dots.width/2,stimulus.dots.width/2,1,n_out);
%     stimulus.dots.y(index_out) = stimulus.dots.ycenter + unifrnd(-stimulus.dots.height/2,stimulus.dots.height/2,1,n_out);
%
%     % For new dots, sample from von Mises to determine pass-through
%     % location
%     theta_new = slCirc_vmrnd2(stimulus.dots.mu,stimulus.dots.keppa,1,n_out); %sample from von Mises
%
%     % calculate pass-through location
%     x_passthru = stimulus.dots.xcenter + r*stimulus.dots.rmax*cos(theta_new);
%     y_passthru = stimulus.dots.ycenter + r*stimulus.dots.rmax*sin(theta_new);
%     stimulus.dots.x_passthru(index_out) = x_passthru;
%     stimulus.dots.y_passthru(index_out) = y_passthru;
%
%     % calculate motion direction
%     d_x = x_passthru - stimulus.dots.x(index_out);
%     d_y = y_passthru - stimulus.dots.y(index_out);
%     dir_new = atan(d_y./d_x);
%
%     index_x_st0 = d_x<0; index_x_gt0 = d_x>0;
%     index_y_st0 = d_y<0; index_y_gt0 = d_y>0;
%     index_Q2 = index_x_st0 + index_y_gt0==2;
%     index_Q3 = index_x_st0 + index_y_st0==2;
%     index_Q4 = index_x_gt0 + index_y_st0==2;
%     dir_new(index_Q2) = dir_new(index_Q2)+pi;
%     dir_new(index_Q3) = dir_new(index_Q3)+pi;
%     dir_new(index_Q4) = dir_new(index_Q4)+2*pi;
%
%     % update direction of motion for new dots
%     stimulus.dots.dir(index_out) = dir_new;
%
%     clear theta_new x_pass y_pass d_x d_y d dir_new index_Q2 index_Q3 index_Q4 index_x_st0 index_x_gt0 index_y_st0 index_y_gt0
% end


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


function stimulus = reDrawDots(stimulus,myscreen)
%
% resample dots motion path for updateDots_2 by resampling the direction of
% motion from von Mises.
%

% sample from von Mises
theta = slCirc_vmrnd2(stimulus.dots.mu,stimulus.dots.keppa,1,stimulus.dots.n_VM_dots); %sample from von Mises
% figure(1);hold on;
% hist(theta)

r=1;
% calculate pass-through location
x_passthru = stimulus.dots.xcenter + r*stimulus.dots.rmax*cos(theta);
y_passthru = stimulus.dots.ycenter + r*stimulus.dots.rmax*sin(theta);
stimulus.dots.x_passthru = x_passthru;
stimulus.dots.y_passthru = y_passthru;

% size(x_passthru)
% size(stimulus.dots.x)

% calculate motion direction
d_x = x_passthru - stimulus.dots.x(stimulus.dots.index_cohM);
d_y = y_passthru - stimulus.dots.y(stimulus.dots.index_cohM);
dir_motion = atan(d_y./d_x);
index_x_st0 = d_x<0; index_x_gt0 = d_x>0;
index_y_st0 = d_y<0; index_y_gt0 = d_y>0;
index_Q2 = index_x_st0 + index_y_gt0 ==2;
index_Q3 = index_x_st0 + index_y_st0 ==2;
index_Q4 = index_x_gt0 + index_y_st0 ==2;
% n_Q2=sum(index_Q2)
% n_Q3=sum(index_Q3)
% n_Q4=sum(index_Q4)
dir_motion(index_Q2) = dir_motion(index_Q2)+pi;
dir_motion(index_Q3) = dir_motion(index_Q3)+pi;
dir_motion(index_Q4) = dir_motion(index_Q4)+2*pi;

% direction of motion; each dot has its own direction
stimulus.dots.dir = dir_motion;


% % % Original version of update dots
% % %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function to update dot positions and draw them to screen
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function stimulus = updateDots(stimulus,myscreen)
% %
% % % get the dots step
% % stimulus.dots.xstep = cos(pi*stimulus.dots.dir/180)*stimulus.dots.stepsize;
% % stimulus.dots.ystep = sin(pi*stimulus.dots.dir/180)*stimulus.dots.stepsize;
% %
% % % pick a random set of dots
% % stimulus.dots.coherent = rand(1,stimulus.dots.n) < stimulus.dots.coherence;
% %
% % % now move those dots in the right direction
% % stimulus.dots.x(stimulus.dots.coherent) = stimulus.dots.x(stimulus.dots.coherent)+stimulus.dots.xstep;
% % stimulus.dots.y(stimulus.dots.coherent) = stimulus.dots.y(stimulus.dots.coherent)+stimulus.dots.ystep;
% %
% % % randomwalk rule
% % thisdir = rand(1,sum(~stimulus.dots.coherent))*2*pi;
% % stimulus.dots.x(~stimulus.dots.coherent) = stimulus.dots.x(~stimulus.dots.coherent)+cos(thisdir)*stimulus.dots.stepsize;
% % stimulus.dots.y(~stimulus.dots.coherent) = stimulus.dots.y(~stimulus.dots.coherent)+sin(thisdir)*stimulus.dots.stepsize;
% %
% % % movshon noise
% % %stimulus.dots.x(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.width;
% % %stimulus.dots.y(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.height;
% %
% % % make sure we haven't gone off the patch
% % % do the dots separately for left and right hand side
% % stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin) = stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin)+stimulus.dots.width;
% % stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax) = stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax)-stimulus.dots.width;
% % stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin) = stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin)+stimulus.dots.height;
% % stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax) = stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax)-stimulus.dots.height;
% %
% % % draw the dots
% % mglStencilSelect(1);
% % mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize,[1 1 1]);
% % mglStencilSelect(0);
% %
%
%
%
%
%
%
