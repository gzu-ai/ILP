%changed two places
% get the support %
hypotheses(I, H) :- output(I, H, N), N > 0.
support_pos(I, A) :- hypotheses(I, _), output(I, _, N), input(I, A, NA), NA >= N.
support_neg(I, A) :- hypotheses(I, _), universe(A), not input(I, A, _). 
support_margin(I, A) :- hypotheses(I, _), output(I, _, N), input(I, A, NA), NA == N.
support_necessity(I,N) :- hypotheses(I, _), output(I, _, N).

% get the limitation %
limitation(I, H) :- output(I, H, _).
limit_pos(I, A) :- output(I, _, N), input(I, A, NA), NA > N.
limit_neg(I, A) :- limitation(I,_), universe(A), not input(I, A, _).  
limit_necessity(I,N) :- output(I, _, N).

% set the body atoms of the rule candidate. generate only the subset of supporting set  %
{rule_pos(I,A)} :- support_pos(I,A).
{rule_neg(I,A)} :- support_neg(I,A).

%analyze the four conditions of limitation which causes incompatibility, in order to set necessity constraints %
%hypothesis I and limitation J share the same head H
hypo_limit_head(I,J,H) :- hypotheses(I,H),limitation(J,H).
%positive body of I is not a subset of its corresponding limit_pos J. Namely, condition 1 of incompatibility unsatisfied
limit_pos_unsat(I,J) :- rule_pos(I,A),hypo_limit_head(I,J,_),not limit_pos(J,A).
limit_neg_unsat(I,J) :- rule_neg(I,A),hypo_limit_head(I,J,_),not limit_neg(J,A).
%if the body part, viz. condition 1 and 2, of incompatibility is satisfied, the N of I must be constrainted by J. otherwise, free N.
necessity_limited(I,J) :- hypo_limit_head(I,J,_),not limit_pos_unsat(I,J),not limit_neg_unsat(I,J).
%it is not allowed that three conditions of incompatibility are inevitably satisfied simultaneously. 
%This constraint makes the solving faster, although the premised coherency makes this check meaningless.
:- necessity_limited(I,J),limit_necessity(J,N1),support_necessity(I,N2),N1 < N2.
% set the two boundaries of necessity of the rule candidate %
%if lower_bound exists, it must be support_necessity.
lower(I,N) :- support_necessity(I,N).
% set the type of necessity of the rule candidate %
%if positive body separates from supporting margin, the necessity is fixed.
unfixed(I) :- hypotheses(I,_),rule_pos(I,A),support_margin(I,A). %a range
fixed(I) :- hypotheses(I,_),not unfixed(I).%a value
%:- fixed(I),unfixed(I).
%if the necessity is fixed, the upper bound is equal to the support necessity. So it is N
upper(I,N) :- fixed(I),support_necessity(I,N).%changed
%if the necessity is unfixed, and limited, then the only possible value of its upper bound is limit_necessity
%Here, a necessity could be constrainted by many limitation. So we choose the least.
upper(I,N_min) :- unfixed(I),N_min = #min{N: necessity_limited(I,J),limit_necessity(J,N)}.
%if the necessity is unfixed, and not limited, then set the default upper bound top.
upper(I,Y) :- unfixed(I), top(Y).
%This hard constraint is ensured by the way we get upper and the above hard constraint on support_necessity(I,N2)
%:- lower(I,N1),upper(I,N2),N2 < N1.

% cognize combinable rules  %
%hypothesis I1 and hypothesis J2 share the same head H. Rule out reflexivity via an order I1 != I2.
hypo_hypo_head(I1,I2,H) :- hypotheses(I1,H),hypotheses(I2,H),I1 != I2.
%I1 has some positive atom beyond the counterpart of I2.
notin_pos(I1,I2) :- hypo_hypo_head(I1,I2,_),rule_pos(I1,A),not rule_pos(I2,A).
notin_neg(I1,I2) :- hypo_hypo_head(I1,I2,_),rule_neg(I1,A),not rule_neg(I2,A).
%the positive body of I1 is equivalent to the counterpart of I2. Rule out symmetry via an order I1 < I2.
atom_pos_coincide(I1,I2) :- hypo_hypo_head(I1,I2,_),not notin_pos(I1,I2), not notin_pos(I2,I1), I1 < I2.
atom_neg_coincide(I1,I2) :- hypo_hypo_head(I1,I2,_),not notin_neg(I1,I2), not notin_neg(I2,I1), I1 < I2.
%body I1 and body I2 are combinable.
atom_combinable(I1,I2) :- atom_pos_coincide(I1,I2), atom_neg_coincide(I1,I2).
% figure out the group identification and its boundaries %
%the body of I1 acts as an identification of group in which I2 is a member of it.
group_id_member(I1,I2) :- atom_combinable(I1,I2), not atom_combinable(_,I1).
%I is a group id
group_id(I) :- group_id_member(I,_).
%I2 has been covered by some group.
covered_rule(I2):- group_id_member(I1,I2).
%In group I1, get the holistic lower bound L and upper bound U
rule_combined_holistic(I1,L,U) :- group_id(I1), 
L = #max{L2: group_id_member(I1,I2),lower(I2,L2); L1: lower(I1,L1)},
U = #min{U2: group_id_member(I1,I2),upper(I2,U2); U1: upper(I1,U1)}.
%for those single hypothesis I, get the L and U
rule_combined_holistic(I,L,U) :- hypotheses(I,_),not group_id(I),not covered_rule(I), lower(I,L), upper(I,U).
%any holistic rule should satisfy U >= L
:- rule_combined_holistic(_,L,U),L > U.

% take background into consideration %
%hypothesis I and background rule J share the same head H
hypo_bg_head(I,J,H) :- hypotheses(I,H),bg_head(J,H).
%whether the hypothesis I in P  is removable because of the rule J in background B
combinable_p_b(I,J) :- hypo_bg_head(I,J,_),rule_combined_holistic(I,L,U),bg_necessity(J,N), N >= L, N <= U.
%at least one atom in the positive body of I belonging to P is outside the counterpart of J belonging to B
pos_p_notin_b(I,J) :- combinable_p_b(I,J),rule_pos(I,A),not bg_body_pos(J,A).
%at least one atom in the positive body of J belonging to B is outside the counterpart of I belonging to P
pos_b_notin_p(I,J) :- combinable_p_b(I,J),not rule_pos(I,A),bg_body_pos(J,A).
neg_p_notin_b(I,J) :- combinable_p_b(I,J),rule_neg(I,A),not bg_body_neg(J,A).
neg_b_notin_p(I,J) :- combinable_p_b(I,J),not rule_neg(I,A),bg_body_neg(J,A).
%the hypothesis I in P  is removed because of the rule J in background
combined_p_b(I,J) :- combinable_p_b(I,J),not pos_p_notin_b(I,J),not pos_b_notin_p(I,J),not neg_p_notin_b(I,J),not neg_b_notin_p(I,J).
% settle the last rules and its cardinality, as well as atoms count %
%the last rule I whose head is H has the necessity L
last_rule(I,H,L) :- rule_combined_holistic(I,L,_),hypotheses(I,H),not combined_p_b(I,_).


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
