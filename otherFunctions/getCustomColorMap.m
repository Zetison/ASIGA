function cmap = getCustomColorMap()

std       = [0      0       0.5625;
             0      0       1;
             0      0.7       0.7;
             1.5*[44 77 32]/255;
             0.95      0.95       0;
             1      0       0;
             0.5      0       0];
      
n = [6 15 15 15 15 6];

column1 = [];
column2 = [];
column3 = [];
for i = 1:6
    column1 = [column1; std(i,1); linspace2(std(i,1), std(i+1,1), n(i))'];
    column2 = [column2; std(i,2); linspace2(std(i,2), std(i+1,2), n(i))'];
    column3 = [column3; std(i,3); linspace2(std(i,3), std(i+1,3), n(i))'];
end
column1 = [column1; std(end,1)];
column2 = [column2; std(end,2)];
column3 = [column3; std(end,3)];

cmap = [column1 column2 column3];


n = 100;
cmap   = [insertUniform4(cmap(:,1), n)  insertUniform4(cmap(:,2), n)  insertUniform4(cmap(:,3), n)]; 