%%%% get the specific hypothesis%%%%%%%%%
rulePos(I, H, N, A) :- output(I, H, N), input(I, A, NA), NA >= N.
ruleNeg(I, H, N, A) :- universe(A), output(I, H, N), not input(I, A, _). 


%%%%leave the unsubsumed rules  as the last rules only%%%%%%%%%%
%the head must exist, whereas the body does not.
lastRuleHead(J, H, N) :- output(J, H, N).
lastRulePos(I, H, N, A) :- lastRuleHead(I, H, N), rulePos(I, H, N, A).
lastRuleNeg(I, H, N, A) :- lastRuleHead(I, H, N), ruleNeg(I, H, N, A).

#show lastRuleHead/3.
#show lastRulePos/4.
#show lastRuleNeg/4.
% #show rule_pos/4.
% #show rule_neg/4.
%#show last_rule/3.
% #show input/3.
% #show output/3.

