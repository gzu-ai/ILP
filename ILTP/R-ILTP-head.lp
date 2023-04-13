#defined bgNecessity/2.
#defined bgBodyPos/2.
#defined bgBodyNeg/2.

%%%% get the candidate solution space%%%%%%%%%
{rulePos(I, N, A)} :- output(I, N), input(I, A, NA), NA >= N, N != 0.
{ruleNeg(I, N, A)} :- universe(A), output(I, N), not input(I, A, _), N != 0. 



%%%% get the psuedo solution space%%%%%%%%%
%get the elements of psuedo solution spaces
limitation(I, N) :- output(I, N).
limitPos(I, N, A) :- limitation(I, N), input(I, A, NA), NA > N.
limitNeg(I, N, A) :- universe(A), limitation(I, N), not input(I, A, _). 

%%%% no rule occures in a  psuedo solution space%%%%%%%%%%%%%
limitStrictlySubsume(I, J, N1, N2) :- output(I, N1), limitation(J, N2), N1 > N2.
notLimitSubsetPos(I, J, N1, N2) :- limitStrictlySubsume(I, J, N1, N2), rulePos(I, N1, A), not limitPos(J, N2, A).
notLimitSubsetNeg(I, J, N1, N2) :- limitStrictlySubsume(I, J, N1, N2), ruleNeg(I, N1, A), not limitNeg(J, N2, A).
:- limitStrictlySubsume(I, J, N1, N2), not notLimitSubsetPos(I, J, N1, N2), not notLimitSubsetNeg(I, J, N1, N2).


%%%%leave the unsubsumed rules  as the last rules only%%%%%%%%%%
%the head must exist, whereas the body does not.
lastRuleHead(J, N) :- output(J, N), N != 0.
lastRulePos(I, N, A) :- lastRuleHead(I, N), rulePos(I, N, A).
lastRuleNeg(I, N, A) :- lastRuleHead(I, N), ruleNeg(I, N, A).

#show lastRuleHead/2.
#show lastRulePos/3.
#show lastRuleNeg/3.
%内容与ILTP一致
%#show subsumed/3.
%#show strictSubsumes/5.
%#show mutualSubsumes/5.
%#show winMutualSubsumes/3.
% #show limitStrictlySubsume/5.
% #show rule_pos/4.
% #show rule_neg/4.
%#show last_rule/3.
% #show input/3.
% #show output/3.


% get the three measures %
%the last program has C rules.
last_cardinality(C) :- C = #count{(I, N): lastRuleHead(I, N)}.
%the amount of all atoms appearing in the program.
%atoms_amount(Pos_cnt + Neg_cnt + Head_cnt) :- Pos_cnt = #count{(I, N, A): lastRulePos(I, N, A) },Neg_cnt = #count{(I, N, A): lastRuleNeg(I, N, A) },last_cardinality(Head_cnt).
%the count of atoms in the program
%atoms_cnt(Atoms_cnt) :- Atoms_cnt = #count{A: lastRulePos(I, N, A); A: lastRuleNeg(I, N, A); H:lastRuleHead(J, N2)}.

% point out the optimal directions %
%Change the two into weak constraints
%less rules as possible. This measure is more significant than atoms_amount whose priority is lower.
#minimize{C@3: last_cardinality(C)}.
%holistic less atoms including repetition (shorter rules) as possible
%#minimize{C@2: atoms_amount(C)}.
%holistic less atoms excluding repetition as possible
%#minimize{C@1: atoms_cnt(C)}.