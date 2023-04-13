%changed two places
% get the support %
hypotheses(I, N) :- output(I, N), N > 0.
{rulePos(I,A)} :- hypotheses(I,_), output(I,  N), input(I, A, NA), NA >= N.
{ruleNeg(I,A)} :- hypotheses(I,_), universe(A), not input(I, A, _). 

% get the limitation %
limitation(I, N) :- output(I, N).
limitPos(I, A) :- output(I,  N), input(I, A, NA), NA > N.
limitNeg(I, A) :- limitation(I,_), universe(A), not input(I, A, _).  

%analyze the four conditions of limitation which causes incompatibility, in order to set necessity constraints %
%hypothesis I and limitation J share the same head H
hypoLimitHead(I,J) :- hypotheses(I,_),limitation(J,_).
%positive body of I is not a subset of its corresponding limit_pos J. Namely, condition 1 of incompatibility unsatisfied
limitPosUnsat(I,J) :- rulePos(I,A),hypoLimitHead(I,J),not limitPos(J,A).
limitNegUnsat(I,J) :- ruleNeg(I,A),hypoLimitHead(I,J),not limitNeg(J,A).
%if the body part, viz. condition 1 and 2, of incompatibility is satisfied, the N of I must be constrainted by J. otherwise, free N.
necessityLimited(I,J) :- hypoLimitHead(I,J),not limitPosUnsat(I,J),not limitNegUnsat(I,J).
%it is not allowed that three conditions of incompatibility are inevitably satisfied simultaneously. 
%This constraint makes the solving faster, although the premised coherency makes this check meaningless.
:- necessityLimited(I,J),limitation(J,N1),hypotheses(I,N2),N1 < N2.


% cognize combinable rules  %
%hypothesis I1 and hypothesis J2 share the same head H. Rule out reflexivity via an order I1 != I2.
hypoHead(I1,I2) :- hypotheses(I1,_),hypotheses(I2,_),I1 != I2.
%I1 has some positive atom beyond the counterpart of I2.
notinPos(I1,I2) :- hypoHead(I1,I2),rulePos(I1,A),not rulePos(I2,A).
notinNeg(I1,I2) :- hypoHead(I1,I2),ruleNeg(I1,A),not ruleNeg(I2,A).
%the positive body of I1 is equivalent to the counterpart of I2. Rule out symmetry via an order I1 < I2.
atomPosCoincide(I1,I2) :- hypoHead(I1,I2),not notinPos(I1,I2), not notinPos(I2,I1), I1 < I2.
atomNegCoincide(I1,I2) :- hypoHead(I1,I2),not notinNeg(I1,I2), not notinNeg(I2,I1), I1 < I2.
%body I1 and body I2 are combinable.
atomCombinable(I1,I2) :- atomPosCoincide(I1,I2), atomNegCoincide(I1,I2).
% figure out the group identification and its boundaries %
%the body of I1 acts as an identification of group in which I2 is a member of it.
groupIdMember(I1,I2) :- atomCombinable(I1,I2), not atomCombinable(_,I1).
%I is a group id
groupId(I) :- groupIdMember(I,_).
%I2 has been covered by some group.
coveredRule(I2):- groupIdMember(I1,I2).
%In group I1, get the holistic lower bound L which meets all boundaries
ruleCombinedHolistic(I1,L) :- groupId(I1), hypotheses(I1,N1),  L = #max{N1; N2: groupIdMember(I1, I2), hypotheses(I2,N2)}.

%for those single hypothesis I, get the L and U
ruleCombinedHolistic(I,N) :- hypotheses(I,N), not groupId(I),not coveredRule(I).
%any holistic rule should satisfy U >= L
%:- ruleCombinedHolistic(_,L,U),L > U.


% take background into consideration %
%hypothesis I and background rule J share the same head H
hypoBgHead(I,J) :- hypotheses(I,_),bg_head(J).
%whether the hypothesis I in P  is removable because of the rule J in background B
combinablePB(I,J) :- hypoBgHead(I,J),ruleCombinedHolistic(I,L),bg_necessity(J,N), L <= N.
%at least one atom in the positive body of I belonging to P is outside the counterpart of J belonging to B
posPNotinB(I,J) :- combinablePB(I,J),rulePos(I,A),not bg_body_pos(J,A).
%at least one atom in the positive body of J belonging to B is outside the counterpart of I belonging to P
posBNotinP(I,J) :- combinablePB(I,J),not rulePos(I,A),bg_body_pos(J,A).
negPNotinB(I,J) :- combinablePB(I,J),ruleNeg(I,A),not bg_body_neg(J,A).
negBNotinP(I,J) :- combinablePB(I,J),not ruleNeg(I,A),bg_body_neg(J,A).
%the hypothesis I in P  is removed because of the rule J in background
combinedPB(I,J) :- combinablePB(I,J),not posPNotinB(I,J),not posBNotinP(I,J),not negPNotinB(I,J),not negBNotinP(I,J).


% settle the last rules and its cardinality, as well as atoms count %
%the last rule I whose head is H has the necessity L
lastRuleHead(I,L) :- ruleCombinedHolistic(I,L),hypotheses(I,_),not combinedPB(I,_).
lastRulePos(I, N, A) :- lastRuleHead(I, N), rulePos(I,A).
lastRuleNeg(I, N, A) :- lastRuleHead(I, N), ruleNeg(I,A).

% get the three measures %
%the last program has C rules.
lastCardinality(C) :- C = #count{(I, H, N): lastRuleHead(I, H, N)}.
%the amount of all atoms appearing in the program.
%atomsAmount(PosCnt + NegCnt + HeadCnt) :- PosCnt = #count{(I, H, N, A): lastRulePos(I, H, N, A) },NegCnt = #count{(I, H, N, A): lastRuleNeg(I, H, N, A) }, HeadCnt = #count{(I, H, N): lastRuleHead(I, H, N)}.
%the count of atoms in the program
%atomsCnt(AtomsCnt) :- AtomsCnt = #count{A: lastRulePos(I, H, N, A); A: lastRuleNeg(I, H, N, A); H:lastRuleHead(J, H, N2)}.

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



