function animation_TRH(filename,FigHandle,first,pause_t)
% Purpose:  This function is designed to be called repeatedly and save a
%           series of MATLAB figures as an animated .gif
% Author:   Tobias Hasse    tobiack@udel.edu
% Date:     Spring, 2016, edited June 2021

if exist( 'pause_t' , 'var' )
else                % if the frame rate pause_t is not set choos 1/2 second
    pause_t = .5; 
end
figure(FigHandle)
drawnow

frame = getframe(FigHandle);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256,'nodither');% convert rgb image to 'flat' image

if first == 1;                  % make the file with the first frame
    imwrite(imind,cm,filename,'gif','DelayTime',2*pause_t,'Loopcount',inf);
else                            % add frames
    imwrite(imind,cm,filename,'gif','DelayTime',pause_t,...
        'WriteMode','append');
end

end
