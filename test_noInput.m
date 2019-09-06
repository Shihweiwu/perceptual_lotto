function test_noInput
%
%
%

global stimulus

stimulus = struct;
stimulus.mousetrack = rand(1,10);

dispInfo(stimulus)

function dispInfo()

mt = stimulus.mousetrack;

