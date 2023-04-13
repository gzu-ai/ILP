
% universe(1 .. 3). 

% input(e1,1,60). 
% input(e1,2,20).
% output(e1,3,20).

% input(e2,2,30). 
% output(e2,3,30).

% input(e3,2,50). 
% output(e3,3,30).

% input(e4,2,60). 
% output(e4,3,80).

% universe(1 .. 2). 

% input(e1,1,60). 
% %output(e1,2,40).

% input(e2,1,30). 
% output(e2,2,30).

universe(1 .. 3). 


input(e1,1,60). 
input(e1,2,20).
output(e1,3,20).
output(e1,2,20).

input(e2,2,30). 
output(e2,3,30).
output(e2,2,20).

input(e3,2,50). 
output(e3,3,30).%inconsistent
output(e3,2,20).

%get the concomitant outputs
coOutput(I, A, 0) :- input(I, _, _), universe(A), not output(I, A, _).
anyOutput(I, A, N) :- coOutput(I, A, N).
anyOutput(I, A, N) :- output(I, A, N).

%the out necessity of J is less than the counterpart of K
less(J,K,A) :- anyOutput(J,A,N1),anyOutput(K,A,N2), N1 < N2.
%the input of J has A, but K doesn't. Namely, the input of J is not a subset of K. Therefore, confliction with the third requirement of incoherent
notNegSat(J,K,H) :- less(J,K,H),input(J,A,_),not input(K,A,_).
%the atoms of the inputs of J, whose necessity is greater than the settled output.
greatK(J,K,H,A) :- less(J,K,H),anyOutput(J,H,N1),input(J,A,N3),N3 > N1.
%the atoms of the inputs of K, whose necessity is greater than or equal to the settled output.
greatEqJ(J,K,H,A) :- less(J,K,H),anyOutput(K,H,N2),input(K,A,N4),N4 >= N2.
%conflict with the second requirement of incoherent.
notPosSat(J,K,H) :- greatEqJ(J,K,H,A) ,not greatK(J,K,H,A).
%three requirement of incoherency are satisfied, a contradiction
contradiction(J,K,H) :- less(J,K,H),not notNegSat(J,K,H), not notPosSat(J,K,H).
%this contradiction means incoherency
incoherent :- contradiction(_,_,_).
%if not incoherent, than coherent.
coherent :- not incoherent.


#show incoherent/0.
#show coherent/0.
%*
#show contradiction/2.


%trans(1 .. 15). input(1,a,6). output(1,b,0). input(2,a,6). output(2,a,0). input(3,a,6). output(3,s,6). input(4,a,3). input(4,s,4). output(4,b,4). input(5,a,3). input(5,s,4). output(5,a,0). input(6,a,3). input(6,s,4). output(6,s,1). input(7,s,7). output(7,b,5). input(8,s,7). output(8,a,3). input(9,s,7). output(9,s,1). input(10,b,2). output(10,b,0). input(11,b,2). output(11,a,0). input(12,b,2). output(12,s,1). input(13,b,8). output(13,b,0). input(14,b,8). output(14,a,0). input(15,b,8). output(15,s,1).

trans(16). input(16,s,6). output(16,a,4).

trans(8).input(8,s,7). output(8,a,3).

*%