%%%% get the specific hypothesis%%%%%%%%%
rulePos(I, N, A) :- output(I, N), input(I, A, NA), NA >= N, N != 0.
ruleNeg(I, N, A) :- universe(A), output(I, N), not input(I, A, _), N != 0.


%%%subsumption between examples%%%%%%%
% get the relation I subsumes J
subsumeNec(I, J, N1, N2) :- output(I, N1), output(J, N2), I != J, N1 >= N2.
notSubsetPos(I, J, N1, N2) :- subsumeNec(I, J, N1, N2), rulePos(I, N1, A), not rulePos(J, N2, A).
notSubsetNeg(I, J, N1, N2) :- subsumeNec(I, J, N1, N2), ruleNeg(I, N1, A), not ruleNeg(J, N2, A).
subsumes(I, J, N1, N2) :- subsumeNec(I, J, N1, N2), not notSubsetPos(I, J, N1, N2), not notSubsetNeg(I, J, N1, N2).
% two kinds of sumsumptions
strictSubsumes(I, J, N1, N2) :- subsumes(I, J, N1, N2), not subsumes(J, I, N2, N1).
mutualSubsumes(I, J, N, N) :- subsumes(I, J, N, N), subsumes(J, I, N, N).%mutual subsuming relation is equivalence. So N1 == N2.
%the most genereal rule in the mutual subsumptions
winMutualSubsumes(C, N) :- output(I, N), C = #min{I; J: mutualSubsumes(I, J, N, _)}, #count{J: mutualSubsumes(I, J, N, _)} > 0.
%pick out the subsumed rules
subsumed(I, N) :- mutualSubsumes(I, J, N, N), not winMutualSubsumes(I, N).
subsumed(J, N2) :- strictSubsumes(I, J, N1, N2).


%%%subsumption between examples and backgrounds%%%%%%%%%%
% get the relation I subsumes J which corresponds to an example
bgSubsumeNec(I, J, N1, N2) :- bgRuleHead(I, N1), output(J, N2), N1 >= N2.
notBgSubsetPos(I, J, N1, N2) :- bgSubsumeNec(I, J, N1, N2), bgRulePos(I, A), not rulePos(J, N2, A).
notBgSubsetNeg(I, J, N1, N2) :- bgSubsumeNec(I, J, N1, N2), bgRuleNeg(I, A), not ruleNeg(J, N2, A).
subsumed(J, N2) :- bgSubsumeNec(I, J, N1, N2), not notBgSubsetPos(I, J, N1, N2), not notBgSubsetNeg(I, J, N1, N2).



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

