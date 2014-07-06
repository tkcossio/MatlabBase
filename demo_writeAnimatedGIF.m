clearvars, close all, clc

figure(1)
set(1, 'Visible', 'off')

for idx = 1:1:100;
    y = phantom(256,256) + randn(256,256) .* 0.5;
    y = y - min(y(:));
    y = y ./ max(y(:));
    imagesc(y), drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if idx==1;
        imwrite(imind,cm,'test.gif','gif','Loopcount',inf,'DelayTime', 0.01);
    else
        imwrite(imind,cm,'test.gif','gif','DelayTime',0.01,'WriteMode','append');
    end
end

close(1)