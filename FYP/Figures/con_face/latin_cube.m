function lat_cube = latin_cube(N)
%LATIN_CUBE Generates a randomized NxNxN latin cube.
%   Given a positive integer specifying size (N), this code will compute
%   an example latin cube, a 3 dimensional extension of the latin square in
%   which each orthogonal slice through the cube is itself a latin square.
%   Note that the cube created is not an example of a /balanced/ latin
%   cube, rather a randomized one.  
%
%   Example: 
%   >>M = latin_cube(3)
%     M(:,:,1) =
% 
%          3 1 2
%          2 3 1
%          1 2 3
% 
%     M(:,:,2) =
% 
%          1 2 3
%          3 1 2
%          2 3 1
% 
%     M(:,:,3) =
% 
%          2 3 1
%          1 2 3
%          3 1 2
%
%   Author: W. Owen Brimijoin
%   MRC Institute of Hearing Research (Scottish Section)
%   owen(at)ihr.gla.ac.uk
%
%   Date: 7th July, 2010
%   with a hearty acknowledgment to Jos van der Geest's latsq.m
%==========================================================================

if floor(N)~=N | N<1
    display('Error: Cube size must be positive integer')
    lat_cube = [];
    return
end

[lat_cube, lat_square] = deal(rem(cumsum([1:N;ones(N-1,N)])-1,N)+1); %create a latin square

for n = 1:N-1
    lat_cube = cat(3,lat_cube,circshift(lat_square,[n 0])); %append shifted versions to the stack
end

lat_cube = lat_cube(randperm(N),randperm(N),randperm(N)); %randomize rows, columns, and stacks.

%the end
