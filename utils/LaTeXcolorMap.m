function map = LaTeXcolorMap(M)

map = [0,80,158; % ntnublue
       178,0,0; % myRed (black!30!red)
       59,124,37; % OliveGreen
       149,49,157; % Purple
       247, 158,30; % YellowOrange
       0,172,239; % Cyan
       236,0,140; % Magenta
       182,50,28; % BrickRed
       0,169,154; % JungleGreen
       140,54,140; % Fuchsia
       247,146,29; % BurntOrange
       70,179,221; % SkyBlue
       148,150,152; % Gray
       152,204,112; % YellowGreen
       239,88,160; % VioletRed
       175,114,176; % Orchid
       ]/255;
   
if M > 16
    map = hsv(M);
else
    map = map(1:M,:);
end