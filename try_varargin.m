function try_varargin( varargin )

global stimulus fixStimulus

stimulus = struct;
fixStimulus = struct;

stimulus.rotSpd = 90;

%% Initialize Variables

% add arguments later
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

getArgs(varargin,{'scan=0','cue=1','plots=0','noeye=0','powerwheel=1','eyewindow=3','practice=0','practiceType=-1','debug=0','replay=0','run=0','build=0','mouse=0'})

