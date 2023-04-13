%%%% get the candidate solution space%%%%%%%%%
1 {lastRuleHead(J, H, N); trimedRuleHead(J, H, N)} 1 :- output(J, H, N).
{rulePos(I, H, A)} :- output(I, H, N), input(I, A, NA), NA >= N.
{ruleNeg(I, H, A)} :- universe(A), output(I, H,_), not input(I, A, _). 



%%%% get the psuedo solution space%%%%%%%%%
%get the concomitant outputs
coOutput(I, A, 0) :- input(I, _, _), universe(A), not output(I, A, _).
%get the elements of psuedo solution spaces
limitation(I, H, N) :- coOutput(I, H, N).
limitation(I, H, N) :- output(I, H, N).
limitPos(I, H, A) :- limitation(I, H, N), input(I, A, NA), NA > N.
limitNeg(I, H, A) :- universe(A), limitation(I, H, N), not input(I, A, _). 

%%%% no rule occures in a  psuedo solution space%%%%%%%%%%%%%
limitStrictlySubsume(I, J, H) :- output(I, H, N1), limitation(J, H, N2), N1 > N2.
notLimitSubsetPos(I, J, H) :- limitStrictlySubsume(I, J, H), rulePos(I, H, A), not limitPos(J, H, A).
notLimitSubsetNeg(I, J, H) :- limitStrictlySubsume(I, J, H), ruleNeg(I, H, A), not limitNeg(J, H, A).
:- limitStrictlySubsume(I, J, H), not notLimitSubsetPos(I, J, H), not notLimitSubsetNeg(I, J, H).



hypoHead(I1,I2,H) :- output(I1,H,N1),output(I2,H,N2),I1 != I2,N1 >= N2.
notinPos(I1,I2,H) :- hypoHead(I1,I2,H),rulePos(I1,H,A),not rulePos(I2,H,A).
notinNeg(I1,I2,H) :- hypoHead(I1,I2,H),ruleNeg(I1,H,A),not ruleNeg(I2,H,A).
atomPosCoincide(I1,I2,H) :- hypoHead(I1,I2,H),not notinPos(I1,I2,H), not notinPos(I2,I1,H).
atomNegCoincide(I1,I2,H) :- hypoHead(I1,I2,H),not notinNeg(I1,I2,H), not notinNeg(I2,I1,H).
atomCombinable(I1,I2,H) :- lastRuleHead(I1, H, _), trimedRuleHead(I2, H, _),atomPosCoincide(I1,I2,H), atomNegCoincide(I1,I2,H).



% % get the relation that last rule I subsumes trimed rule J
% exSubsumeNec(I, J, H) :- lastRuleHead(I, H, N1), trimedRuleHead(J, H, N2), N1 >= N2.

% notAtomPosSubset(I, J, H) :- exSubsumeNec(I, J, H), rulePos(I, H, A), not rulePos(J, H, A).
% notAtomPosSupset(I, J, H) :- exSubsumeNec(I, J, H), rulePos(J, H, B), not rulePos(I, H, B).
% atomPosCoincide(I, J, H) :- exSubsumeNec(I, J, H), not notAtomPosSubset(I, J, H), not notAtomPosSupset(I, J, H).
% notAtomNegSubset(I, J, H) :- exSubsumeNec(I, J, H), ruleNeg(I, H, A), not ruleNeg(J, H, A).
% notAtomNegSupset(I, J, H) :- exSubsumeNec(I, J, H), ruleNeg(J, H, B), not ruleNeg(I, H, B).
% atomNegCoincide(I, J, H) :- exSubsumeNec(I, J, H), not notAtomNegSubset(I, J, H), not notAtomNegSupset(I, J, H).
% atomCombinable(I, J, H) :- exSubsumeNec(I, J, H), atomPosCoincide(I, J, H), atomNegCoincide(I, J, H).


% atomPosCoincide(I, J, H, N1, N2) :- exSubsumeNec(I, J, H, N1, N2), rulePos(I, H, N1, A), not rulePos(J, H, N2, A), rulePos(J, H, N2, B), not rulePos(I, H, N1, B).
% atomNegCoincide(I, J, H, N1, N2) :- exSubsumeNec(I, J, H, N1, N2), ruleNeg(I, H, N1, A), not ruleNeg(J, H, N2, A), ruleNeg(J, H, N2, B), not ruleNeg(I, H, N1, B).
% atomCombinable(I, J, H, N1, N2) :- atomPosCoincide(I, J, H, N1, N2), atomNegCoincide(I, J, H, N1, N2).

% strictSubsums(I1,H,N1) :- atomCombinable(I1,I2,H, N1, N2), N1 > N2.
% mutualSubsums(I1,H,N1) :- atomCombinable(I1,I2,H, N1, N2), N1 == N2.
% winMutualSubsum(C,H,N) :- mutualSubsums(I, H, N), C = #min{I; J: mutualSubsums(J, H, N)}.
% subsumed(I2, H, N2) :- atomCombinable(I1,I2,H, N1, N2), not strictSubsums(I2,H,N2), not winMutualSubsum(I2,H,N2).

subsumed(I2, H) :- atomCombinable(I1,I2,H).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % get the relation that last rule I subsumes trimed rule J
% exSubsumeNec(I, J, H, N1, N2) :- lastRuleHead(I, H, N1), trimedRuleHead(J, H, N2), N1 >= N2.
% notExSubsetPos(I, J, H, N1, N2) :- exSubsumeNec(I, J, H, N1, N2), rulePos(I, H, N1, A), not rulePos(J, H, N2, A).
% notExSubsetNeg(I, J, H, N1, N2) :- exSubsumeNec(I, J, H, N1, N2), ruleNeg(I, H, N1, A), not ruleNeg(J, H, N2, A).
% subsumed(J, H, N2) :- exSubsumeNec(I, J, H, N1, N2), not notExSubsetPos(I, J, H, N1, N2), not notExSubsetNeg(I, J, H, N1, N2).

% % get the relation that background rule I subsumes trimed rule J
% bgSubsumeNec(I, J, H, N1, N2) :- bgRuleHead(I, H, N1), trimedRuleHead(J, H, N2), N1 >= N2.
% notBgSubsetPos(I, J, H, N1, N2) :- bgSubsumeNec(I, J, H, N1, N2), bgRulePos(I, A), not rulePos(J, H, N2, A).
% notBgSubsetNeg(I, J, H, N1, N2) :- bgSubsumeNec(I, J, H, N1, N2), bgRuleNeg(I, A), not ruleNeg(J, H, N2, A).
% subsumed(J, H, N2) :- bgSubsumeNec(I, J, H, N1, N2), not notBgSubsetPos(I, J, H, N1, N2), not notBgSubsetNeg(I, J, H, N1, N2).


%a single transition must has corresponding rule at last if it is not covered by any existed rule.
:- trimedRuleHead(I, H, _), not subsumed(I, H).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%leave the unsubsumed rules  as the last rules only%%%%%%%%%%
%the head must exist, whereas the body does not.
%{lastRuleHead(J, H, N)} :- output(J, H, N).
lastRulePos(I, H, N, A) :- lastRuleHead(I, H, N), rulePos(I, H, N, A).
lastRuleNeg(I, H, N, A) :- lastRuleHead(I, H, N), ruleNeg(I, H, N, A).

#show lastRuleHead/3.
#show lastRulePos/4.
#show lastRuleNeg/4.
%内容与ILTP一致
%#show subsumed/3.
% #show input/3.
% #show output/3.
% #show subsumed/3.
% #show rulePos/4.
% #show ruleNeg/4.

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

%:- lastRuleHead(e2,6,56).

% input(e1,3,27).
% input(e1,4,86).
% input(e1,6,43).
% input(e1,9,84).
% output(e1,4,27).

% input(e2,2,77).
% input(e2,4,62).
% input(e2,5,28).
% input(e2,7,38).
% input(e2,9,90).
% input(e2,10,77).
% output(e2,4,28).

% lastRuleHead(e1,6,75).
% lastRuleHead(e2,6,56).
% lastRulePos(e1,6,75,8).
