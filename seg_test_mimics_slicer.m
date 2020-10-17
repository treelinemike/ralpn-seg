% restart
close all; clear all; clc;

% stl_mimics_raw = stlread('H:\CT\31584-001\RK_RAW.stl');

stl_mimics_raw = stlread('H:\CT\manual_seg\31584-004\RK_RAW_mimics_unsmoothed.stl');
stl_slicer_raw = stlread('H:\CT\manual_seg\31584-004\RK_RAW_slicer_unsmoothed.stl');

% stl_mimics_raw = stlread('H:\CT\manual_seg\31584-001\RK_RAW_mimics_unsmoothed.stl');
% stl_slicer_raw = stlread('H:\CT\manual_seg\31584-001\RK_RAW_slicer_unsmoothed.stl');

% stl_mimics = stlread('C:/users/f002r5k/Desktop/ralpn_local/31584-001 CT/D19002-RALPN_RK_001.stl');
% stl_slicer = stlread('C:\Users\f002r5k\GitHub\ralpn-seg\manual_seg\31584-001\RK_001_smoothed.stl');

figure;
hold on; grid on;
axis equal;
% plot3(stl_mimics.Points(:,1),stl_mimics.Points(:,2),stl_mimics.Points(:,3)-2.5,'.','MarkerSize',2,'Color',[0.8 0.2 0.2]);
% plot3(stl_slicer.Points(:,1),stl_slicer.Points(:,2),stl_slicer.Points(:,3),'.','MarkerSize',2,'Color',[0.2 0.8 0.2]);
% 
plot3(stl_mimics_raw.Points(:,1),stl_mimics_raw.Points(:,2),stl_mimics_raw.Points(:,3)-2.5,'.','MarkerSize',2,'Color',[0.2 0.2 0.8]);
plot3(stl_slicer_raw.Points(:,1),stl_slicer_raw.Points(:,2),stl_slicer_raw.Points(:,3),'.','MarkerSize',2,'Color',[0.8 0.2 0.8]);

view([90,0]);

% stl_mimics.Points(1,:)-stl_slicer.Points(1,:)
% 
% mean(stl_mimics.Points)-mean(stl_slicer.Points)