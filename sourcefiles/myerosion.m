function filtered = myerosion(imageIN, valthreshold)



if isempty(imageIN),        error('Empty image');       end;
if isempty(valthreshold),   error('Empty threshold');   end;

[m,n] = size(imageIN);

pts2 = sub2ind([m n],[1 1 1 1 2 2 (m-1) (m-1) m m m m],[1 2 (n-1) n 1 n 1 n 1 2 (n-1) n]);
val2 = imageIN(pts2);

im1med = medfilt2(imageIN,[5 5]);   

diffIm = abs(imageIN - im1med);
diffIm(pts2) = 0;

pts1 = find(diffIm > valthreshold);

imageIN(pts1) = im1med(pts1);
imageIN(pts2) = val2;

filtered = imageIN;

disp(['Number of points corrected with myFilt = ' num2str(numel(pts1)) ' over ' num2str(numel(imageIN)) ' (' num2str(numel(pts1)/numel(imageIN)*100) ' %) \n']);
fprintf('\n\n');