clc
clear 
close all

m = [6
8
10
12
14
16
18
20
22
24
26
28
30
];

N = [14.1562
14.1338
14.1232
14.1175
14.114
14.1117
14.1102
14.1091
14.1083
14.1076
14.1071
14.1068
14.1064

];


figure
hold on; box on; grid on;

plot(m,N,'.b','MarkerSize',24,'LineWidth',2);
plot(m,N,'-b','MarkerSize',24,'LineWidth',2);
ff = gca;
ff.FontSize = 24;
ff.FontName = 'Times New Roman';
xlabel('m')
ylabel('N_{critical}')





