pixels = [ones(1,5000) zeros(1,5000)];

A = reshape(pixels,[100 100]);

A_f = getHalfFourier(A);
B_f = A_f;
B_f.phase = rand(1,length(B_f.phase))*pi*2;
B = reconstructFromHalfFourier(B_f);
B = B + abs(min(B(:)));
B = B * 1 / max(B(:));

figure
subplot(211)
imagesc(A,[0 1])
colormap('gray')
subplot(212)
imagesc(B,[0 1])
colormap('gray')

% we now have an image, A, and a phase scrambled version, B. We are going
% to combine them 'on the circle'.