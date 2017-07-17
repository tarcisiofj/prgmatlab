function [Result]=Function_ColorToGray(CitraInput)

%get red channel
R_Channel=CitraInput(:,:,1);

%figure; imshow(R_Channel); title('R_Channel');

%get green channel
G_Channel=CitraInput(:,:,2);

%figure; imshow(G_Channel); title('G_Channel');

%get blue channel
B_Channel=CitraInput(:,:,3);

%figure; imshow(B_Channel); title('B_Channel');

% get gray image with lightness
GrayImage=(max(max(R_Channel, G_Channel), B_Channel) ...
    + min(min(R_Channel, G_Channel), B_Channel))./2;

%figure; imshow(GrayImage); title('GrayImage lightness');

% get gray image with average
GrayImage=(R_Channel+G_Channel+B_Channel)./3;

%figure; imshow(GrayImage); title('GrayImage average');

% get gray image with luminosity
GrayImage=0.2989 * R_Channel + 0.5870 * G_Channel + 0.1140 * B_Channel;

%figure; imshow(GrayImage); title('GrayImage luminosity');
Result=GrayImage;

