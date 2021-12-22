clear all;
clc; 

load('/Users/anette/Documents/MATLAB/Experiment2/Subdata/TimingTest/STARTSubject_TimingTest3_1.mat');

AudCond = find(o.targets.audTargets == 1);

WkA(1,:) = (o.timing.audTarOnset(1,AudCond) - o.timing.wakeup(2,AudCond))*1000;
WkA(2,:) = (o.timing.audTarOnset(2,AudCond) - o.timing.wakeup(4,AudCond))*1000;

WkTr(1,:) = (o.targets.TriggerOnsets(1,:) - o.timing.wakeup(2,:))*1000;
WkTr(2,:) = (o.targets.TriggerOnsets(2,:) - o.timing.wakeup(4,:))*1000;

TrA(1,:) = (o.timing.audTarOnset(1,AudCond) - o.targets.TriggerOnsets(1,AudCond))*1000;
TrA(2,:) = (o.timing.audTarOnset(2,AudCond) - o.targets.TriggerOnsets(2,AudCond))*1000;

TrV(1,:) = (o.timing.flipTarOnset(1,:) - o.targets.TriggerOnsets(1,:))*1000;
TrV(2,:) = (o.timing.flipTarOnset(2,:) - o.targets.TriggerOnsets(2,:))*1000;

SOAa = (o.timing.audTarOnset(2,AudCond) - o.timing.audTarOnset(1,AudCond))*1000;
SOAv = (o.timing.flipTarOnset(2,:) - o.timing.flipTarOnset(1,:))*1000;

