% Starter script for CSCI 1290 project on Poisson blending. 
% Written by James Hays and Pat Doran.
% imblend.m is the function where you will implement your blending method.
% By default, imblend.m performs a direct composite.

close all
clear variables

data_dir = './data';
out_dir = './results';

%there are four inputs for each compositing operation -- 
% 1. 'source' image. Parts of this image will be inserted into 'target'
% 2. 'mask' image. This binary image, the same size as 'source', specifies
%     which pixels will be copied to 'target'
% 3. 'target' image. This is the destination for the 'source' pixels under
%     the 'mask'
% 4. 'offset' vector. This specifies how much to translate the 'source'
%     pixels when copying them to 'target'. These vectors are hard coded
%     below for the default test cases. They are of the form [y, x] where
%     positive values mean shifts down and to the right, respectively.

offset = cell(13,1);
offset{1} = [ 210  10 ];
offset{2} = [  10  28 ];
offset{3} = [ 140 80 ];
offset{4} = [  -40  90 ];
offset{5} = [  60 100 ];
offset{6} = [ -28  88 ];
offset{7} = [ 100  88 ];
offset{8} = [ -130  -60 ];
offset{9} = [ 55  100 ];
offset{10} = [ 1000  200 ];
offset{11} = [ 1400  650 ];
offset{12} = [ 750  200 ];
offset{13} = [ 10  180 ];
offset{14} = [ 600  180 ];
j=[1 2 3 4 5 6 10 12 13 14];
for z = 10:10%length(offset)
    i =j(z) 
    source = imread(sprintf('%s/source_%02d.jpg',data_dir,i));
    %mask   = imread(sprintf('%s/mask_%02d.jpg',data_dir,i));
    target = imread(sprintf('%s/target_%02d.jpg',data_dir,i));
    
    source = im2double(source);
    %mask = round(im2double(mask));
    target = im2double(target);

    % Interactive mask
    mask = getmask(source);
    %target = [0.2 0.5 0.2 0.2; 0.7 .7 .7 .7; .9 .9 .8 .90];
    %mask = [0 0 0 0 ; 0 1 1 0; 0 0 0 0];
    imwrite(mask,sprintf('%s/mask_%02d.jpg',out_dir,i),'jpg','Quality',95);
    figure, imshow(mask)
    [source, mask, target] = fiximages(source, mask, target, offset{i});
    
%     output = imblend_mixgrad_max(source, mask, target, 1);
%     output = imblend_mixgrad_max(source, mask, target, 1);
    output = imblend(source, mask, target, 1);
    figure(i)
    imshow(output)
    
    imwrite(output,sprintf('%s/gradmax_result_%02d.jpg',out_dir,i),'jpg','Quality',95);
end

