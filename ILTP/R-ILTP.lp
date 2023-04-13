%%%% get the candidate solution space%%%%%%%%%
{rulePos(I, H, N, A)} :- output(I, H, N), input(I, A, NA), NA >= N.
{ruleNeg(I, H, N, A)} :- universe(A), output(I, H, N), not input(I, A, _). 



%%%% get the psuedo solution space%%%%%%%%%
%get the concomitant outputs
coOutput(I, A, 0) :- input(I, _, _), universe(A), not output(I, A, _).
%get the elements of psuedo solution spaces
limitation(I, H, N) :- coOutput(I, H, N).
limitation(I, H, N) :- output(I, H, N).
limitPos(I, H, N, A) :- limitation(I, H, N), input(I, A, NA), NA > N.
limitNeg(I, H, N, A) :- universe(A), limitation(I, H, N), not input(I, A, _). 

%%%% no rule occures in a  psuedo solution space%%%%%%%%%%%%%
limitStrictlySubsume(I, J, H, N1, N2) :- output(I, H, N1), limitation(J, H, N2), N1 > N2.
notLimitSubsetPos(I, J, H, N1, N2) :- limitStrictlySubsume(I, J, H, N1, N2), rulePos(I, H, N1, A), not limitPos(J, H, N2, A).
notLimitSubsetNeg(I, J, H, N1, N2) :- limitStrictlySubsume(I, J, H, N1, N2), ruleNeg(I, H, N1, A), not limitNeg(J, H, N2, A).
:- limitStrictlySubsume(I, J, H, N1, N2), not notLimitSubsetPos(I, J, H, N1, N2), not notLimitSubsetNeg(I, J, H, N1, N2).


%内容与ILTP一致
%%%subsumption between examples%%%%%%%
% get the relation I subsumes J
subsumeNec(I, J, H, N1, N2) :- output(I, H, N1), output(J, H, N2), I != J, N1 >= N2.
notSubsetPos(I, J, H, N1, N2) :- subsumeNec(I, J, H, N1, N2), rulePos(I, H, N1, A), not rulePos(J, H, N2, A).
notSubsetNeg(I, J, H, N1, N2) :- subsumeNec(I, J, H, N1, N2), ruleNeg(I, H, N1, A), not ruleNeg(J, H, N2, A).
subsumes(I, J, H, N1, N2) :- subsumeNec(I, J, H, N1, N2), not notSubsetPos(I, J, H, N1, N2), not notSubsetNeg(I, J, H, N1, N2).
% two kinds of sumsumptions
strictSubsumes(I, J, H, N1, N2) :- subsumes(I, J, H, N1, N2), not subsumes(J, I, H, N2, N1).
mutualSubsumes(I, J, H, N, N) :- subsumes(I, J, H, N, N), subsumes(J, I, H, N, N).%mutual subsuming relation is equivalence. So N1 == N2.
%the most genereal rule in the mutual subsumptions
winMutualSubsumes(C, H, N) :- output(I, H, N), C = #min{I; J: mutualSubsumes(I, J, H, N, _)}, #count{J: mutualSubsumes(I, J, H, N, _)} > 0.
%pick out the subsumed rules
subsumed(I, H, N) :- mutualSubsumes(I, J, H, N, N), not winMutualSubsumes(I, H, N).
subsumed(J, H, N2) :- strictSubsumes(I, J, H, N1, N2).


%%%subsumption between examples and backgrounds%%%%%%%%%%
% get the relation I subsumes J which corresponds to an example
bgSubsumeNec(I, J, H, N1, N2) :- bgRuleHead(I, H, N1), output(J, H, N2), N1 >= N2.
notBgSubsetPos(I, J, H, N1, N2) :- bgSubsumeNec(I, J, H, N1, N2), bgRulePos(I, A), not rulePos(J, H, N2, A).
notBgSubsetNeg(I, J, H, N1, N2) :- bgSubsumeNec(I, J, H, N1, N2), bgRuleNeg(I, A), not ruleNeg(J, H, N2, A).
subsumed(J, H, N2) :- bgSubsumeNec(I, J, H, N1, N2), not notBgSubsetPos(I, J, H, N1, N2), not notBgSubsetNeg(I, J, H, N1, N2).



%%%%leave the unsubsumed rules  as the last rules only%%%%%%%%%%
%the head must exist, whereas the body does not.
lastRuleHead(J, H, N) :- output(J, H, N), not subsumed(J, H, N).
lastRulePos(I, H, N, A) :- lastRuleHead(I, H, N), rulePos(I, H, N, A).
lastRuleNeg(I, H, N, A) :- lastRuleHead(I, H, N), ruleNeg(I, H, N, A).

#show lastRuleHead/3.
#show lastRulePos/4.
#show lastRuleNeg/4.
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
last_cardinality(C) :- C = #count{(I, H, N): lastRuleHead(I, H, N)}.
%the amount of all atoms appearing in the program.
%atoms_amount(Pos_cnt + Neg_cnt + Head_cnt) :- Pos_cnt = #count{(I, H, N, A): lastRulePos(I, H, N, A) },Neg_cnt = #count{(I, H, N, A): lastRuleNeg(I, H, N, A) },last_cardinality(Head_cnt).
%the count of atoms in the program
%atoms_cnt(Atoms_cnt) :- Atoms_cnt = #count{A: lastRulePos(I, H, N, A); A: lastRuleNeg(I, H, N, A); H:lastRuleHead(J, H, N2)}.

% point out the optimal directions %
%Change the two into weak constraints
%less rules as possible. This measure is more significant than atoms_amount whose priority is lower.
#minimize{C@3: last_cardinality(C)}.
%holistic less atoms including repetition (shorter rules) as possible
%#minimize{C@2: atoms_amount(C)}.
%holistic less atoms excluding repetition as possible
%#minimize{C@1: atoms_cnt(C)}.