% set up screen
myscreen = initScreen('test');
% init screen parameters

% set background to black
myscreen.background = 0;

mglClearScreen(0);
% mglClearScreen([clearColor], [clearBits])
% Sets the background color and clears the buffer.

mglStencilCreateBegin(1)
% mglStencilCreateBegin(stencilNumber,invert)
% Begin drawing to stencil.Until mglStencilCreateEnd
% is called, all drawing operations will also draw
% to the stencil.

mglFillOval(300,400,[100 100]);
mglStencilCreateEnd;
mglClearScreen;

% now draw some dots, masked by the oval stencil
mglStencilSelect(1);
mglPoints2(rand(1,5000)*500,rand(1,5000)*500);
mglFlush;
mglStencilSelect(0);

% Fill oval
% mglVisualAngleCoordinates(57,[16 12]);
% x = [-1 -4 -3 0 3 4 1];
% y = [-1 -4 -3 0 3 4 1];
% sz = [1 1]; 
% mglFillOval(x, y, sz,  [1 0 0]);
% % mglFillOval(x,y, size, color)
% %x,y - vectors of x and y coordinates of center
% %size - [width height] of oval
% %color - [r g b] triplet for color
% 
% mglFlush();