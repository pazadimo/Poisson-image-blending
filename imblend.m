function  output = imblend( source, mask, target, transparent )
%Source, mask, and target are the same size (as long as you do not remove
%the call to fiximages.m). You may want to use a flag for whether or not to
%treat the source object as 'transparent' (e.g. taking the max gradient
%rather than the source gradient).
flag=0;
t = target;
if(sum(mask(1,:,1))>0 || sum(mask(:,1,1))>0 || sum(mask(size(mask,1),:,1))>0 || sum(mask(:,size(mask,2),1))>0)
    S=size(source);
    source2=zeros(S(1)+2,S(2)+2,3);
    source2(2:(S(1)+1), 2:(S(2)+1),:) = source;
    source = source2;
    target2=zeros(S(1)+2,S(2)+2,3);
    target2(2:(S(1)+1), 2:(S(2)+1),:) = target;
    target=target2;
    mask2=zeros(S(1)+2,S(2)+2,3);
    mask2(2:(S(1)+1), 2:(S(2)+1),:) = mask;
    mask =mask2;
    flag=1;
end
Size = size(target);
[j_mask_index,i_mask_index] = find(mask(:,:,1) == 1);

%[j_mask_index,i_mask_index] = find(mask(:,:) == 1);

%A =speye(Size(1)*Size(2))*4;
%sparse(diag(4*ones(1,Size(1)*Size(2))));
indexes = sub2ind(Size(1:2),j_mask_index,i_mask_index);
all_row = ones((Size(1)*Size(2)),1);
all_row(:,1) = 1:(Size(1)*Size(2));
total_ind_row = cat(1,all_row,indexes,indexes,indexes,indexes,indexes);

indexes_up = sub2ind(Size,j_mask_index-1,i_mask_index);
indexes_down = sub2ind(Size,j_mask_index+1,i_mask_index);
indexes_right = sub2ind(Size,j_mask_index,i_mask_index+1);
indexes_left = sub2ind(Size,j_mask_index,i_mask_index-1);
total_ind_col = cat(1,all_row,indexes_up,indexes_down,indexes_right,indexes_left,indexes);

v = ones(length(indexes)*5 , 1)*-1;
v((length(indexes)*4+1):length(indexes)*5 , 1) = v((length(indexes)*4+1):length(indexes)*5 , 1)*-3;
all_row = ones((Size(1)*Size(2)),1);
v = cat(1,all_row,v);
disp(size(total_ind_row));
disp(size(total_ind_col));
disp(size(v));
A= sparse(total_ind_row,total_ind_col,v,Size(1)*Size(2),Size(1)*Size(2));
% A(indexes,indexes_up ) = -1;
% A(indexes,indexes_down ) = -1;
% A(indexes,indexes_right ) = -1;
% A(indexes,indexes_left ) = -1;
sourcered = source(:,:,1);
sourcered = sourcered(:);
sourcegreen = source(:,:,2);
sourcegreen = sourcegreen(:);
sourceblue = source(:,:,3);
sourceblue = sourceblue(:);
bred = target(:,:,1);
bred = bred(:);
bred(indexes) = 4*sourcered(indexes)-sourcered(indexes_up)-sourcered(indexes_down)-sourcered(indexes_left)-sourcered(indexes_right);
bgreen = target(:,:,2);
bgreen = bgreen(:);
bgreen(indexes) = 4*sourcegreen(indexes)-sourcegreen(indexes_up)-sourcegreen(indexes_down)-sourcegreen(indexes_left)-sourcegreen(indexes_right);

bblue = target(:,:,3);
bblue = bblue(:);
bblue(indexes) = 4*sourceblue(indexes)-sourceblue(indexes_up)-sourceblue(indexes_down)-sourceblue(indexes_left)-sourceblue(indexes_right);

%b(indexes) 

outred = A\bred;
outgreen = A\bgreen;
outblue = A\bblue;
out1= reshape(outred,Size(1),Size(2));
out2= reshape(outgreen,Size(1),Size(2));
out3= reshape(outblue,Size(1),Size(2));

if(flag ==0)
   output = ones(Size(1),Size(2),3);
   output(:,:,1)=out1;
   output(:,:,2)=out2;
   output(:,:,3)=out3; 
else
   output = ones(Size(1)-2,Size(2)-2,3);
   output(:,:,1)=out1(2:Size(1)-1,2:Size(2)-1,1);
   output(:,:,2)=out2(2:Size(1)-1,2:Size(2)-1,1);
   output(:,:,3)=out3(2:Size(1)-1,2:Size(2)-1,1); 
end
output = (output * transparent + t*(1 - transparent));
%output = source .* mask + target .* ~mask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As explained on the web page, we solve for output by setting up a large
% system of equations, in matrix form, which specifies the desired value or
% gradient or Laplacian (e.g.
% http://en.wikipedia.org/wiki/Discrete_Laplace_operator)

% The comments here will walk you through a conceptually simple way to set
% up the image blending, although it is not necessarily the most efficient
% formulation. 

% We will set up a system of equations A * x = b, where A has as many rows
% and columns as there are pixels in our images. Thus, a 300x200 image will
% lead to A being 60000 x 60000. 'x' is our output image (a single color
% channel of it) stretched out as a vector. 'b' contains two types of known 
% values:
%  (1) For rows of A which correspond to pixels that are not under the
%      mask, b will simply contain the already known value from 'target' 
%      and the row of A will be a row of an identity matrix. Basically, 
%      this is our system of equations saying "do nothing for the pixels we 
%      already know".
%  (2) For rows of A which correspond to pixels under the mask, we will
%      specify that the gradient (actually the discrete Laplacian) in the
%      output should equal the gradient in 'source', according to the final
%      equation in the webpage:
%         4*x(i,j) - x(i-1, j) - x(i+1, j) - x(i, j-1) - x(i, j+1) = 
%         4*s(i,j) - s(i-1, j) - s(i+1, j) - s(i, j-1) - s(i, j+1)
%      The right hand side are measurements from the source image. The left
%      hand side relates different (mostly) unknown pixels in the output
%      image. At a high level, for these rows in our system of equations we
%      are saying "For this pixel, I don't know its value, but I know that
%      its value relative to its neighbors should be the same as it was in
%      the source image".

% commands you may find useful: 
%   speye - With the simplest formulation, most rows of 'A' will be the
%      same as an identity matrix. So one strategy is to start with a
%      sparse identity matrix from speye and then add the necessary
%      values. This will be somewhat slow.
%   sparse - if you want your code to run quickly, compute the values and
%      indices for the non-zero entries in A and then construct 'A' with a
%      single call to 'sparse'.
%      Matlab documentation on what's going on under the hood with a sparse
%      matrix: www.mathworks.com/help/pdf_doc/otherdocs/simax.pdf
%   reshape - convert x back to an image with a single call.
%   sub2ind and ind2sub - how to find correspondence between rows of A and
%      pixels in the image. It's faster if you simply do the conversion
%      yourself, though.
%   see also find, sort, diff, cat, and spy


