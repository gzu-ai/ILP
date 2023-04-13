%changed two places
% get the support %
hypotheses(I, H) :- output(I, H, N), N > 0.
supportPos(I, A) :- hypotheses(I, _), output(I, _, N), input(I, A, NA), NA >= N.
supportNeg(I, A) :- hypotheses(I, _), universe(A), not input(I, A, _). 
supportNecessity(I,N) :- hypotheses(I, _), output(I, _, N).

% get the limitation %
limitation(I, H) :- output(I, H, _).
limitPos(I, A) :- output(I, _, N), input(I, A, NA), NA > N.
limitNeg(I, A) :- limitation(I,_), universe(A), not input(I, A, _).  
limitNecessity(I,N) :- output(I, _, N).

% set the body atoms of the rule candidate. generate only the subset of supporting set  %
{rule_pos(I,A)} :- supportPos(I,A).
{rule_neg(I,A)} :- supportNeg(I,A).

%analyze the four conditions of limitation which causes incompatibility, in order to set necessity constraints %
%hypothesis I and limitation J share the same head H
hypoLimitHead(I,J,H) :- hypotheses(I,H),limitation(J,H).
%positive body of I is not a subset of its corresponding limit_pos J. Namely, condition 1 of incompatibility unsatisfied
limitPosUnsat(I,J) :- rule_pos(I,A),hypoLimitHead(I,J,_),not limitPos(J,A).
limitNegUnsat(I,J) :- rule_neg(I,A),hypoLimitHead(I,J,_),not limitNeg(J,A).
%if the body part, viz. condition 1 and 2, of incompatibility is satisfied, the N of I must be constrainted by J. otherwise, free N.
necessityLimited(I,J) :- hypoLimitHead(I,J,_),not limitPosUnsat(I,J),not limitNegUnsat(I,J).
%it is not allowed that three conditions of incompatibility are inevitably satisfied simultaneously. 
%This constraint makes the solving faster, although the premised coherency makes this check meaningless.
:- necessityLimited(I,J),limitNecessity(J,N1),supportNecessity(I,N2),N1 < N2.


% cognize combinable rules  %
%hypothesis I1 and hypothesis J2 share the same head H. Rule out reflexivity via an order I1 != I2.
hypoHead(I1,I2) :- hypotheses(I1,H),hypotheses(I2,H),I1 != I2.
%I1 has some positive atom beyond the counterpart of I2.
notinPos(I1,I2) :- hypoHead(I1,I2),rule_pos(I1,A),not rule_pos(I2,A).
notinNeg(I1,I2) :- hypoHead(I1,I2),rule_neg(I1,A),not rule_neg(I2,A).
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
ruleCombinedHolistic(I1,L) :- groupId(I1), supportNecessity(I1,N1),  L = #max{N1; N2: groupIdMember(I1, I2), supportNecessity(I2,N2)}.

%for those single hypothesis I, get the L and U
ruleCombinedHolistic(I,N) :- hypotheses(I,_),supportNecessity(I,N), not groupId(I),not coveredRule(I).
%any holistic rule should satisfy U >= L
%:- ruleCombinedHolistic(_,L,U),L > U.


% take background into consideration %
%hypothesis I and background rule J share the same head H
hypoBgHead(I,J,H) :- hypotheses(I,H),bg_head(J,H).
%whether the hypothesis I in P  is removable because of the rule J in background B
combinablePB(I,J) :- hypoBgHead(I,J,_),ruleCombinedHolistic(I,L),bg_necessity(J,N), L <= N.
%at least one atom in the positive body of I belonging to P is outside the counterpart of J belonging to B
posPNotinB(I,J) :- combinablePB(I,J),rule_pos(I,A),not bg_body_pos(J,A).
%at least one atom in the positive body of J belonging to B is outside the counterpart of I belonging to P
posBNotinP(I,J) :- combinablePB(I,J),not rule_pos(I,A),bg_body_pos(J,A).
negPNotinB(I,J) :- combinablePB(I,J),rule_neg(I,A),not bg_body_neg(J,A).
negBNotinP(I,J) :- combinablePB(I,J),not rule_neg(I,A),bg_body_neg(J,A).
%the hypothesis I in P  is removed because of the rule J in background
combinedPB(I,J) :- combinablePB(I,J),not posPNotinB(I,J),not posBNotinP(I,J),not negPNotinB(I,J),not negBNotinP(I,J).


% settle the last rules and its cardinality, as well as atoms count %
%the last rule I whose head is H has the necessity L
last_rule(I,H,L) :- ruleCombinedHolistic(I,L),hypotheses(I,H),not combinedPB(I,_).


% get the three measures %
%the last program has C rules.count once for the entire group, while count once for every rule outside a group.
last_cardinality(C) :- C = #count{I: last_rule(I,_,_)}.
%the amount of all atoms appearing in the program.
%atoms_amount(Pos_cnt + Neg_cnt + Head_cnt) :- Pos_cnt = #count{(I,A): rule_pos(I,A),last_rule(I,_,_) },Neg_cnt = #count{(I,A): rule_neg(I,A),last_rule(I,_,_)},last_cardinality(Head_cnt).
%the count of atoms in the program
%atoms_cnt(Atoms_cnt) :- Atoms_cnt = #count{A: last_rule(I,_,_),rule_pos(I,A); A: last_rule(I,_,_),rule_neg(I,A); H:last_rule(I,_,_),hypotheses(I,H)}.


% point out the optimal directions %
%Change the two into weak constraints
%less rules as possible. This measure is more significant than atoms_amount whose priority is lower.
#minimize{C@3: last_cardinality(C)}.
%holistic less atoms including repetition (shorter rules) as possible
%#minimize{C@2: atoms_amount(C)}.
%holistic less atoms excluding repetition as possible
%#minimize{C@1: atoms_cnt(C)}.


#show rule_pos/2.
#show rule_neg/2.
#show last_rule/3.
