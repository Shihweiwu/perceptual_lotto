function [stimulus,yourInputs] = test_varargin( varargin )
%
% Understanding varargin and getArgs
%
% e.g. [stimulus,yours]=test_varargin('scan',1,'cue',0)
%
%
% Shih-Wei Wu

% scan = 0;
% plots = 0;
% noeye = 0;
% debug = 1;
% replay = 0;
% powerwheel = 0;
% run = 0; 
% eyewindow=0; 
% mouse=0; 
% practice=0; 
% practiceType=-1;
% cue=0;

yourInputs = getArgs(varargin,{'scan=0','cue=1','plots=0','noeye=0','powerwheel=1','eyewindow=3',...
    'practice=0','practiceType=-1','debug=0','replay=0','run=0','build=0','mouse=0'});


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
% Attempt to execute SCRIPT debug as a function:
% /Applications/MATLAB_R2015a.app/toolbox/matlab/codetools/debug.m
% 
% Error in test_varargin (line 36)
% stimulus.debug = debug; %debug? default is no.
%  

stimulus.replay = replay; %replay? not sure what it is for. default is no.
stimulus.overrideRun = run; %run? not sure what this is for. default is no. 
% Error using run
% Too many output arguments.
% 
% Error in test_varargin (line 36)
% stimulus.overrideRun = run; %run? not sure what this is for. default is no. 
