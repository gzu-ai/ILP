hypotheses(I, N) :- output(I, N), N > 0.
{rulePos(I,A)} :- hypotheses(I,_), output(I,  N), input(I, A, NA), NA >= N.
{ruleNeg(I,A)} :- hypotheses(I,_), universe(A), not input(I, A, _). 

limitation(I, N) :- output(I, N).
limitPos(I, A) :- output(I,  N), input(I, A, NA), NA > N.
limitNeg(I, A) :- limitation(I,_), universe(A), not input(I, A, _).  
hypoLimitHead(I,J) :- hypotheses(I,_),limitation(J,_).
limitPosUnsat(I,J) :- rulePos(I,A),hypoLimitHead(I,J),not limitPos(J,A).
limitNegUnsat(I,J) :- ruleNeg(I,A),hypoLimitHead(I,J),not limitNeg(J,A).
necessityLimited(I,J) :- hypoLimitHead(I,J),not limitPosUnsat(I,J),not limitNegUnsat(I,J).
:- necessityLimited(I,J),limitation(J,N1),hypotheses(I,N2),N1 < N2.

hypoHead(I1,I2) :- hypotheses(I1,_),hypotheses(I2,_),I1 != I2.
notinPos(I1,I2) :- hypoHead(I1,I2),rulePos(I1,A),not rulePos(I2,A).
notinNeg(I1,I2) :- hypoHead(I1,I2),ruleNeg(I1,A),not ruleNeg(I2,A).
atomPosCoincide(I1,I2) :- hypoHead(I1,I2),not notinPos(I1,I2), not notinPos(I2,I1), I1 < I2.
atomNegCoincide(I1,I2) :- hypoHead(I1,I2),not notinNeg(I1,I2), not notinNeg(I2,I1), I1 < I2.
atomCombinable(I1,I2) :- atomPosCoincide(I1,I2), atomNegCoincide(I1,I2).
groupIdMember(I1,I2) :- atomCombinable(I1,I2), not atomCombinable(_,I1).
groupId(I) :- groupIdMember(I,_).
coveredRule(I2):- groupIdMember(I1,I2).
ruleCombinedHolistic(I1,L) :- groupId(I1), hypotheses(I1,N1),  L = #max{N1; N2: groupIdMember(I1, I2), hypotheses(I2,N2)}.
ruleCombinedHolistic(I,N) :- hypotheses(I,N), not groupId(I),not coveredRule(I).

hypoBgHead(I,J) :- hypotheses(I,_),bg_head(J).
combinablePB(I,J) :- hypoBgHead(I,J),ruleCombinedHolistic(I,L),bg_necessity(J,N), L <= N.
posPNotinB(I,J) :- combinablePB(I,J),rulePos(I,A),not bg_body_pos(J,A).
posBNotinP(I,J) :- combinablePB(I,J),not rulePos(I,A),bg_body_pos(J,A).
negPNotinB(I,J) :- combinablePB(I,J),ruleNeg(I,A),not bg_body_neg(J,A).
negBNotinP(I,J) :- combinablePB(I,J),not ruleNeg(I,A),bg_body_neg(J,A).
combinedPB(I,J) :- combinablePB(I,J),not posPNotinB(I,J),not posBNotinP(I,J),not negPNotinB(I,J),not negBNotinP(I,J).

lastRuleHead(I,L) :- ruleCombinedHolistic(I,L),hypotheses(I,_),not combinedPB(I,_).
lastRulePos(I, N, A) :- lastRuleHead(I, N), rulePos(I,A).
lastRuleNeg(I, N, A) :- lastRuleHead(I, N), ruleNeg(I,A).

% get the three measures %
%the last program has C rules.
lastCardinality(C) :- C = #count{(I, H): lastRuleHead(I, H)}.
%the amount of all atoms appearing in the program.
%atomsAmount(PosCnt + NegCnt + HeadCnt) :- PosCnt = #count{(I, N, A): lastRulePos(I, N, A) },NegCnt = #count{(I, N, A): lastRuleNeg(I, N, A) }, HeadCnt = #count{(I, H): lastRuleHead(I, H)}.
%the count of atoms in the program
%atomsCnt(AtomsCnt) :- AtomsCnt = #count{A: lastRulePos(I, N, A); B: lastRuleNeg(I, N, B); H:lastRuleHead(J, H)}.

% point out the optimal directions %
%Change the two into weak constraints
%less rules as possible. This measure is more significant than atomsAmount whose priority is lower.
#minimize{C@3: lastCardinality(C)}.

%holistic less atoms including repetition (shorter rules) as possible
%#minimize{C@2: atomsAmount(C)}.
%holistic less atoms excluding repetition as possible
%#minimize{C@1: atomsCnt(C)}.

#show lastRuleHead/2.
#show lastRulePos/3.
#show lastRuleNeg/3.



