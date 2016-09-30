close all
clear all
clc

varstr = {'x'}
varstr2 = {'x'  'y'}
varstr3 = {'x'  'y'  'zp' }

dataW = struct('x',1:2:4)
dataW2 = struct('x',1:2:4,'y',5:3:15)
dataW3 = struct('x',1:2:4,'y',5:3:15,'zp',10)

%%
writeDataToFile(dataW,varstr,'test.txt')

%%
writeDataToFile(dataW2,varstr2,'test2.txt')

%%
writeDataToFile(dataW3,varstr3,'test3.txt')