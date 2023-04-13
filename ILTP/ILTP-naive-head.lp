%%%% get the specific hypothesis%%%%%%%%%
rulePos(I, N, A) :- output(I, N), input(I, A, NA), NA >= N, N != 0.
ruleNeg(I, N, A) :- universe(A), output(I, N), not input(I, A, _), N != 0.


%%%%leave the unsubsumed rules  as the last rules only%%%%%%%%%%
%the head must exist, whereas the body does not.
lastRuleHead(J, N) :- output(J, N), not subsumed(J, N), N != 0.
lastRulePos(I, N, A) :- lastRuleHead(I, N), rulePos(I, N, A).
lastRuleNeg(I, N, A) :- lastRuleHead(I, N), ruleNeg(I, N, A).

#show lastRuleHead/2.
#show lastRulePos/3.
#show lastRuleNeg/3.
% #show rule_pos/4.
% #show rule_neg/4.
%#show last_rule/3.
% #show input/3.
% #show output/3.

