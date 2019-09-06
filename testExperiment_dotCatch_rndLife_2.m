function myscreen = testExperiment_dotCatch_rndLife_2(screenName)
%
% 2019.05.30.SWW.Modified from testExperiment_dotCatch.
% > Limit dots lifespan (e.g. 6 frames for 60Hz refresh rate = 100ms) so
% that after every 6 frames we will redraw all dots and resample their
% motion path.
% > Each dot starts with random lifespan [1,n] where n is the maximum lifespan
% (e.g. 6 frames). After a dot's life goes to 0 it gets reset to the max lifespan.
% defined in task{2}{1}.parameter.lifespan.
% > In this version, at each frame, before we plot the dots we identify which dots
% die. We reset those dots -- giving it max lifespan and randomly re-positioning its 
% location -- and we plot all dots. Hence the number of dots plotted would be the same
% across all frames. testExperiment_dotCatch_rndLife is the version that
% plots only the alive dots.
%
% 2019.05.04.SWW.Change the task such that
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
myscreen.saveData = 1;
myscreen.datadir = cd;
myscreen.subjectID = 'sub01';

% set the first task to be fixation task
global fixStimulus
fixStimulus.pospix = [myscreen.screenWidth myscreen.screenHeight]/2;
[task{1},myscreen] = fixation(myscreen);

% the second task will have two segments -- moving dots and fixation
task{2}{1}.segmin = [1.5 1];
task{2}{1}.segmax = [1.5 1];
% These two parameters, mu and keppa, are the two parameters of the von
% Mises distribution
task{2}{1}.parameter.mu = (0:30:360)*(pi/180); %mean parameter in radian
task{2}{1}.parameter.keppa = [1 10 50 100 200]; %dispersion parameter. greater value means smaller SD.
task{2}{1}.random = 1; %randomize the order in which you present different combinations
task{2}{1}.parameter.lifespan = 12; %set dots lifespan in number of frames

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
stimulus.dots.lifespan = task.thistrial.lifespan; 

% initialize dots life by randomly sampling from 1:lifespan
stimulus.dots.life = randsample(stimulus.dots.lifespan,stimulus.dots.n,1);

% get initial position
% stimulus.dots.x = stimulus.dots.width*(rand(1,stimulus.dots.n) - 0.5);
% stimulus.dots.y = stimulus.dots.height*(rand(1,stimulus.dots.n) - 0.5);
%
% get initial position: randomly place the dots within a circular aperture.
% This will get rid of the "edgy" perception due to placing them within a sqaure aperture. 
r_grid = unifrnd(0,stimulus.dots.rmax,1,stimulus.dots.n);
theta_start = unifrnd(0,2*pi,1,stimulus.dots.n);
stimulus.dots.x = stimulus.dots.xcenter + r_grid.*cos(theta_start); 
stimulus.dots.y = stimulus.dots.ycenter + r_grid.*sin(theta_start); 

% sample from von Mises to determine dots exit location
theta_exit = slCirc_vmrnd2(stimulus.dots.mu,stimulus.dots.keppa,1,stimulus.dots.n);
% figure(1);hold on;
% hist(theta)

r=1; % a scalar just in case. keep it to 1.
% calculate dots exit location
x_exit = stimulus.dots.xcenter + r*stimulus.dots.rmax*cos(theta_exit);
y_exit = stimulus.dots.ycenter + r*stimulus.dots.rmax*sin(theta_exit);
stimulus.dots.x_exit = x_exit; %not sure whether to write it to file.
stimulus.dots.y_exit = y_exit;

% calculate dots motion direction
d_x = x_exit - stimulus.dots.x;
d_y = y_exit - stimulus.dots.y;
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
if ~isfield(stimulus.dots,'speed'), stimulus.dots.speed = 10;,end
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

% get the step size
stimulus.dots.stepsize = stimulus.dots.speed/myscreen.framesPerSecond;

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

% count down dots life
stimulus.dots.life = stimulus.dots.life-1;

% examine whether any dot is dead before countdown
stimulus.dots.dead = stimulus.dots.life==0;

% re-position dead dots, give them new life, and determine its motion path
if sum(stimulus.dots.dead)>0
    stimulus.dots.life(stimulus.dots.dead) = stimulus.dots.lifespan;
    stimulus = reDrawDots(stimulus,myscreen);
end

% Draw dots
stimulus.dots.xstep = stimulus.dots.stepsize*cos(stimulus.dots.dir);
stimulus.dots.ystep = stimulus.dots.stepsize*sin(stimulus.dots.dir);
    
% move dots
stimulus.dots.x = stimulus.dots.x+stimulus.dots.xstep;
stimulus.dots.y = stimulus.dots.y+stimulus.dots.ystep;    

% draw dots
mglStencilSelect(1);
mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize,[1 1 1]);
%mglPoints2(stimulus.dots.x(stimulus.dots.dead),stimulus.dots.y(stimulus.dots.dead),stimulus.dots.dotsize,[1 0 0]);
mglStencilSelect(0);



function stimulus = reDrawDots(stimulus,myscreen)
%
% re-position dead dots and sets up their motion path for updateDots_2.
%

n_dead = sum(stimulus.dots.dead); %number of dead dots
index_dead = find(stimulus.dots.dead)';
% size(stimulus.dots.x)

% get initial position: randomly place the dots within a circular aperture.
% This will get rid of the "edgy" perception due to placing them within a sqaure aperture. 
r_grid = unifrnd(0,stimulus.dots.rmax,1,n_dead);
theta_position = unifrnd(0,2*pi,1,n_dead);
% size(theta_position);
% size(stimulus.dots.x(index_dead));
stimulus.dots.x(index_dead) = stimulus.dots.xcenter + r_grid.*cos(theta_position); 
stimulus.dots.y(index_dead) = stimulus.dots.ycenter + r_grid.*sin(theta_position); 

% sample from von Mises to determine the exit location
theta_exit = slCirc_vmrnd2(stimulus.dots.mu,stimulus.dots.keppa,1,n_dead); %sample from von Mises
% figure(1);hold on;
% hist(theta)

r=1; % a scaler just in case. keep it to 1.
% calculate dots exit location
x_exit = stimulus.dots.xcenter + r*stimulus.dots.rmax*cos(theta_exit);
y_exit = stimulus.dots.ycenter + r*stimulus.dots.rmax*sin(theta_exit);
stimulus.dots.x_exit = x_exit;
stimulus.dots.y_exit = y_exit;

% calculate dots motion direction
d_x = x_exit - stimulus.dots.x(index_dead);
d_y = y_exit - stimulus.dots.y(index_dead);
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
stimulus.dots.dir(index_dead) = dir_motion;


% % Original version of updatDots
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function to update dot positions and draw them to screen
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function stimulus = updateDots(stimulus,myscreen)
% 
% % get the dots step
% stimulus.dots.xstep = cos(pi*stimulus.dots.dir/180)*stimulus.dots.stepsize;
% stimulus.dots.ystep = sin(pi*stimulus.dots.dir/180)*stimulus.dots.stepsize;
% 
% % pick a random set of dots
% stimulus.dots.coherent = rand(1,stimulus.dots.n) < stimulus.dots.coherence;
% 
% % now move those dots in the right direction
% stimulus.dots.x(stimulus.dots.coherent) = stimulus.dots.x(stimulus.dots.coherent)+stimulus.dots.xstep;
% stimulus.dots.y(stimulus.dots.coherent) = stimulus.dots.y(stimulus.dots.coherent)+stimulus.dots.ystep;
% 
% % randomwalk rule
% thisdir = rand(1,sum(~stimulus.dots.coherent))*2*pi;
% stimulus.dots.x(~stimulus.dots.coherent) = stimulus.dots.x(~stimulus.dots.coherent)+cos(thisdir)*stimulus.dots.stepsize;
% stimulus.dots.y(~stimulus.dots.coherent) = stimulus.dots.y(~stimulus.dots.coherent)+sin(thisdir)*stimulus.dots.stepsize;
% 
% % movshon noise
% %stimulus.dots.x(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.width;
% %stimulus.dots.y(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.height;
% 
% % make sure we haven't gone off the patch
% % do the dots separately for left and right hand side
% stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin) = stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin)+stimulus.dots.width;
% stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax) = stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax)-stimulus.dots.width;
% stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin) = stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin)+stimulus.dots.height;
% stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax) = stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax)-stimulus.dots.height;
% 
% % draw the dots
% mglStencilSelect(1);
% mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize,[1 1 1]);
% mglStencilSelect(0);
% 






