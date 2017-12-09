clear all
clc 
close all

max_samples = 50;

comments1 = {'(N', '(AFIB', '(AFIL', '(N', '(AFIB'};
comments2 = {'(AFIB', '(AFIB', '(AFIL', '(N', '(N'};
comments3 = {'(AFIL', '(AFIB', '(AFIL', '(N', '(AFIB'};
comments4 = {'(N', '(AFIB', '(AFIL', '(N', '(AFIL'};

commentPositions1 = [1, 20, 30, 40, 50];
commentPositions2 = [10, 20, 30, 40, 45];

res1 = get_annotation_vector(max_samples, commentPositions1, comments1);
res2 = get_annotation_vector(max_samples, commentPositions1, comments2);
res3 = get_annotation_vector(max_samples, commentPositions2, comments3);
res4 = get_annotation_vector(max_samples, commentPositions2, comments4);