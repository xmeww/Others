close all; 

te = sum(o.response.tooEarly.trial);
fprintf('\n \n number of button presses that came too early: %i \n \n', te)

% auditory conditions
AVcond = sum(o.conditions == 1);
Acond  = sum(o.conditions == 2);

HT(1,1) = sum(o.response.hit(o.conditions == 1) == 1)/AVcond; % AV hits 
HT(2,1) = sum(o.response.hit(o.conditions == 2) == 1)/Acond; % AV hits 

HT(1,2) = sum(o.response.hit(o.conditions == 1) == -2)/AVcond; % AV miss 
HT(2,2) = sum(o.response.hit(o.conditions == 2) == -2)/Acond; % AV miss 

HT(1,3) = sum(o.response.hit(o.conditions == 1) == 0)/AVcond; % AV hits 
HT(2,3) = sum(o.response.hit(o.conditions == 2) == 0)/Acond; % AV hits 

subplot(2,1,1)
bar(HT)
title('auditory conditions')
legend('hit','miss', 'no Response')
xlabel('condition')
ylabel('response in %')
set(gca, 'xtick', 1:2, 'xticklabel', {'AV', 'A'})
axis([0.5 2.5 0 1])

% no auditory conditions
Vcond = sum(o.conditions == 3);
NScond  = sum(o.conditions == 4);

CR(1,1) = sum(o.response.hit(o.conditions == 3) == -1)/Vcond; % AV hits 
CR(2,1) = sum(o.response.hit(o.conditions == 4) == -1)/NScond; % AV hits 

CR(1,2) = sum(o.response.hit(o.conditions == 3) == 2)/Vcond; % AV miss 
CR(2,2) = sum(o.response.hit(o.conditions == 4) == 2)/NScond; % AV miss 

CR(1,3) = sum(o.response.hit(o.conditions == 3) == 0)/Vcond; % AV hits 
CR(2,3) = sum(o.response.hit(o.conditions == 4) == 0)/NScond; % AV hits 

fg = subplot(2,1,2);
bar(CR);
title('No auditory conditions')
legend('Correct rejection','FA', 'no Response')
xlabel('condition')
ylabel('response in %')
set(gca, 'xtick', 1:2, 'xticklabel', {'V', 'NS'})
axis([0.5 2.5 0 1])

waitfor(fg);

% disp('ok')

% % frequencies
% %freqs = Sd.Tarfreqs;
%  freqs = o.sound.targetFreqs;
% 
% 
% for iFreq = 1:length(freqs)
%     AVcond = sum((o.conditions == 1) & (o.targets.tarFreqs == freqs(iFreq)));
%     Acond  = sum((o.conditions == 2) & (o.targets.tarFreqs == freqs(iFreq)));
% 
%     HTf(1,iFreq) = sum(o.response.hit((o.conditions == 1) & (o.targets.tarFreqs == freqs(iFreq))) == 1)/AVcond; %#ok<*SAGROW> % AV hits 
%     HTf(2,iFreq) = sum(o.response.hit((o.conditions == 2) &(o.targets.tarFreqs == freqs(iFreq))) == 1)/Acond; % AV hits 
% end
% 
% figure
% bar(HTf)
% title('% hit per frequency')
% legend(num2str(freqs'))
% xlabel('condition')
% ylabel('response in %')
% set(gca, 'xtick', 1:2, 'xticklabel', {'AV', 'A'})
% %axis([0.5 5.5 0 1])
