:-use_module(library(lists)).
:-use_module(library(system)).
:-use_module(library(file_systems)).
:-op(900,fy,not).
:-dynamic fapt/3.
:-dynamic interogat/1.
:-dynamic scop/1.
:-dynamic interogabil/3.
:-dynamic regula/3.
:-dynamic intrebare_curenta/3.
:-dynamic statistica/2.

% ~~~~  inceput b)  ~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%predicat dinamic pentru parsarea fisierului info.txt
%acesta va avea urmatorul format: info_vehicul(Nume_vehicul,Descriere,Pasi,Imagine).
:-dynamic info_vehicul/4.

not(P):-P,!,fail.
not(_).

scrie_lista([]):-nl.
scrie_lista([H|T]) :-
write(H), tab(1),
scrie_lista(T).
             
afiseaza_fapte :-
write('Fapte existente în baza de cunostinte:'),
nl,nl, write(' (Atribut,valoare) '), nl,nl,
listeaza_fapte,nl.

listeaza_fapte:-  
fapt(av(Atr,Val),FC,_), 
write('('),write(Atr),write(','),
write(Val), write(')'),
write(','), write(' certitudine '),
FC1 is integer(FC),write(FC1),
nl,fail.
listeaza_fapte.

lista_float_int([],[]).
lista_float_int([Regula|Reguli],[Regula1|Reguli1]):-
(Regula \== utiliz,
Regula1 is integer(Regula);
Regula ==utiliz, Regula1=Regula),
lista_float_int(Reguli,Reguli1).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

un_pas(Rasp,OptiuniUrm,MesajUrm):-scop(Atr),(Rasp \== null,intreaba_acum(Rasp) ; true),
								determina1(Atr,OptiuniUrm,MesajUrm), afiseaza_scop(Atr).

intreaba_acum(Rasp):-intrebare_curenta(Atr,OptiuniV,MesajV),interogheaza1(Rasp,Atr,MesajV,OptiuniV,Istorie),nl,
asserta( interogat(av(Atr,_)) ).

interogheaza1(X,Atr,Mesaj,[da,nu],Istorie) :-
!,de_la_utiliz1(X,Istorie,[da,nu]),
det_val_fc(X,Val,FC),
asserta( fapt(av(Atr,Val),FC,[utiliz]) ).

interogheaza1(VLista,Atr,Mesaj,Optiuni,Istorie) :-
de_la_utiliz1(VLista,Optiuni,Istorie),
assert_fapt(Atr,VLista).


%de_la_utiliz1(+Rasp,?Istorie,+Lista_opt)
de_la_utiliz1(X,Istorie,Lista_opt) :-
proceseaza_raspuns([X],Istorie,Lista_opt).


determina1(Atr,OptiuniUrm,MesajUrm) :-
realizare_scop1(av(Atr,_),_,[scop(Atr)],OptiuniUrm,MesajUrm),!.
determina1(_,_,_).

realizare_scop1(not Scop,Not_FC,Istorie,OptiuniUrm,MesajUrm) :-
realizare_scop1(Scop,FC,Istorie,OptiuniUrm,MesajUrm),
Not_FC is - FC, !.
realizare_scop1(Scop,FC,_,_,_) :-
fapt(Scop,FC,_), !.
realizare_scop1(Scop,FC,Istorie,OptiuniUrm,MesajUrm) :-
pot_interoga1(Scop,Istorie,OptiuniUrm,MesajUrm),
!.

%realizare_scop1(Scop,FC,Istorie,OptiuniUrm,MesajUrm).

realizare_scop1(Scop,FC_curent,Istorie,OptiuniUrm,MesajUrm) :-
fg1(Scop,FC_curent,Istorie,OptiuniUrm,MesajUrm).


pot_interoga1(av(Atr,_),Istorie, Optiuni, Mesaj) :-
not interogat(av(Atr,_)),
interogabil(Atr,Optiuni,Mesaj),
retractall(intrebare_curenta(_,_,_)),
assert(intrebare_curenta(Atr, Optiuni,Mesaj)),!.


pornire1:-retractall(interogat(_)),
retractall(fapt(_,_,_)),
retractall(intrebare_curenta(_,_,_)),
retractall(scop(_)),
retractall(interogabil(_)),
retractall(regula(_,_,_)),
incarca('sist_expert.txt').


fg1(Scop,FC_curent,Istorie,OptiuniUrm,MesajUrm) :-
regula(N, premise(Lista), concluzie(Scop,FC)),
demonstreaza1(N,Lista,FC_premise,Istorie,OptiuniUrm,MesajUrm),
(nonvar(FC), nonvar(FC_premise),ajusteaza(FC,FC_premise,FC_nou),
actualizeaza(Scop,FC_nou,FC_curent,N),
FC_curent == 100; true),!.
fg1(Scop,FC,_,_,_) :- fapt(Scop,FC,_).



demonstreaza1(N,ListaPremise,Val_finala,Istorie,OptiuniUrm,MesajUrm) :-
dem1(ListaPremise,100,Val_finala,[N|Istorie],OptiuniUrm,MesajUrm),!.

dem1([],Val_finala,Val_finala,_,_,_).
dem1([H|T],Val_actuala,Val_finala,Istorie,OptiuniUrm,MesajUrm) :-
realizare_scop1(H,FC,Istorie,OptiuniUrm,MesajUrm),
(nonvar(FC),
Val_interm is min(Val_actuala,FC),
Val_interm >= 20,
dem1(T,Val_interm,Val_finala,Istorie,OptiuniUrm,MesajUrm) ;true).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pornire :-
retractall(interogat(_)),
retractall(fapt(_,_,_)),
retractall(intrebare_curenta(_,_,_)),
%pus la inceput ca sa existe baza de date cu statistica(X,Y) in momentul folosirii lui interogheaza
creare_director,
repeat,
write('Introduceti o optiune de mai jos: '),
nl,nl,
write(' (Incarca / Consulta / Reinitiaza / Afisare_fapte / Cum / Iesire) '),
nl,nl,write('|: '),citeste_linie([H|T]),
executa([H|T]), H == iesire.

executa([incarca]) :- 
incarca,!, nl,
write('Fisierul dorit a fost incarcat'),nl,
creare_fisier_statistica,
% citeste fisierul statistica_atribute
citeste_fis_r('output_vehicul/statistica_atribute.txt'). 

executa([consulta]) :- 
scopuri_princ, intreaba_descriere,!.
%%%%%%%%% scopuri_princ,!. %%%%%%%%%%%%

executa([reinitiaza]) :-
citeste_fis_r('output_vehicul/statistica_atribute.txt'), %%%%%
retractall(interogat(_)),
retractall(fapt(_,_,_)),!.

executa([afisare_fapte]) :-
afiseaza_fapte,!.

executa([cum|L]) :- cum(L),!.

executa([iesire]):-!.

executa([pasi]) :-
afiseaza_pasi,!.

executa([descriere]) :-
afiseaza_descriere,!.

executa([ambele]) :-
afiseaza_ambele,!.

executa([inapoi_meniu]) :- !.

executa([_|_]) :-
write('Comanda incorecta! '),nl.


% ~~~~  e)  ~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%creaza directorul in care se salveaza fisierele cu statistica_atribute.txt si demonstratiile
%directorul este format doar daca nu exista 
creare_director:- ((\+directory_exists('output_vehicul'))->(make_directory('output_vehicul'));true).

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%creaza fisierul cu statistica_atribute.txt
%fisierul este format doar daca nu exista
%o sa il apelam dupa incarca 
creare_fisier_statistica:- ((\+file_exists('output_vehicul/statistica_atribute.txt'))->(fis_statistica);true).

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%creaza fisierul 'statistica_atribute.txt' si initializeaza cu 0 numarul aparitiilor atributelor
fis_statistica:-
findall(st(Atr,0),( interogabil(Atr,_,_) ),L_st),
telling(Vechi),
tell('output_vehicul/statistica_atribute.txt'),
scrie_lista_st(L_st),
told,
tell(Vechi).

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%actualizeaza fisierul statistica_atribute.txt
actualizare_statistica:-
findall(st(Atr,Nr),retract(statistica(Atr,Nr)),L_st),
telling(Vechi),
tell('output_vehicul/statistica_atribute.txt'),
scrie_lista_st(L_st),
told,
tell(Vechi).

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%extrage datele din lista cu structurile st(Atr,Nr) si afiseaza Atr si Nr
scrie_lista_st([]):-nl.
scrie_lista_st([st(Atr,Nr)|T]) :-
write(Atr), write('.'), tab(3), write(Nr), write('.'), nl,
scrie_lista_st(T).

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%citeste termen cu termen din fisier (cu redirectare)
%citeste_fis_r(+Fisin)
citeste_fis_r(Fisin):- 
seeing(Input_curent),
see(Fisin),
repeat,
    read(X),
	((X==end_of_file)->(true);(read(Y),assert(statistica(X,Y)),fail)),
!,
seen,
see(Input_curent).

% ~~~~ sfarsit e) ~~~~


% ~~~~  f)  ~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%contruieste numele fisierului demonstratie
nume_fisier(NumeFisier,Val,FC) :-
atom_concat('output_vehicul/demonstratie###',Val,Rez),
atom_concat(Rez,'###',Rez1),
number_chars(FC,N),
atom_chars(NR,N),
atom_concat(Rez1,NR,Rez2),
atom_concat(Rez2,'.txt',NumeFisier).

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% sterge fisierele care incep cu 'demonstratie' din folderul output_vehicul
stergere_demonstratii:-
file_members_of_directory(output_vehicul,Lista_nume_fisiere),
decupare_nume_fisier(Lista_nume_fisiere,Lista_decupata),
parcurgere_lista(Lista_decupata).

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% decupeaza din Lista_nume_fisiere, care contine BaseName-FullName, doar BaseName-FullName
decupare_nume_fisier([A-_|T],[H_dec|T_dec]):- H_dec = A, decupare_nume_fisier(T,T_dec).
decupare_nume_fisier([],[]).

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%parcurge toata lista cu nume de fisiere
parcurgere_lista([H|RestL]):- ((prefix_atom(demonstratie,H))->(atom_concat('output_vehicul/',H,H1),delete_file(H1));true), parcurgere_lista(RestL).
parcurgere_lista([]).

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%verifica daca numele fierului incepe cu Atom1 ( demonstratie )

prefix_atom(Atom1,Atom2):-
% L1 o sa fie [d,e,m,o,n,s,t,r,a,t,i,e]
atom_chars(Atom1,L1),
atom_chars(Atom2,L2),
append(L1,_,L2).

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%crearea fisierelor cu demostratii

afiseaza_demonstratii :- 
scop(Atr),
telling(Vechi),
(fapt(av(Atr,Val),FC,_),
nume_fisier(NumeFisier,Val,FC),
tell(NumeFisier),
cum(av(Atr,Val)),
told,fail ; true),
tell(Vechi).

afiseaza_pasi :- 
scop(Atr),
(fapt(av(Atr,Val),FC,_),
atom_concat('Solutie: ', Val, R1),
    number_codes(FC,ListNr),
	atom_codes(Atom,ListNr),
atom_concat(R1, '  FC:', R2),
atom_concat(R2, Atom, R3),
nl, write(R3), nl,nl,nl,	
info_vehicul(Val,_,Pasi,_),
listeaza_pasi(Pasi,0),
nl, write('##########################'),nl,
fail ; true)
.

afiseaza_descriere :- 
scop(Atr),
(fapt(av(Atr,Val),FC,_),
	atom_concat('Solutie: ', Val, R1),
    number_codes(FC,ListNr),
	atom_codes(Atom,ListNr),
atom_concat(R1, '  FC:', R2),
atom_concat(R2, Atom, R3),
nl, write(R3), nl,nl,nl,
info_vehicul(Val,Desc,_,_),
write(Desc),nl,
nl, write('##########################'),nl,
fail ; true)
.

afiseaza_ambele :- 
scop(Atr),
(fapt(av(Atr,Val),FC,_),
	atom_concat('Solutie: ', Val, R1),
    number_codes(FC,ListNr),
	atom_codes(Atom,ListNr),
atom_concat(R1, '  FC:', R2),
atom_concat(R2, Atom, R3),
nl, write(R3), nl,nl,nl,
info_vehicul(Val,Desc,Pasi,_),
write(Desc),nl,
listeaza_pasi(Pasi,0),
nl, write('##########################'),nl,
fail ; true)
.

listeaza_pasi([],Num).
listeaza_pasi([H|T],Num):- Num1 is Num + 1, 
                           number_codes(Num1,ListNr),
                           atom_codes(Atom,ListNr),
                           atom_concat(Atom, '.<', R1),
                           atom_concat(R1,H,R2),
                           atom_concat(R2,'>', Pas),
                           write(Pas), nl, listeaza_pasi(T,Num1).

scopuri_princ :-
scop(Atr),determina(Atr), afiseaza_scop(Atr), info_vehicul(Atr,Desc,Pasi,Img),fail.
%apeleaza predicatul care scrie demonstatiile in fisiere
scopuri_princ:- stergere_demonstratii, afiseaza_demonstratii, actualizare_statistica. 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%intreaba_descriere 

intreaba_descriere:- 
repeat,
nl,nl,
write(' Afiseaza detalii? (pasi / descriere / ambele / inapoi_meniu) '),
nl,nl,write('|: '),citeste_linie([H|T]),
executa([H|T]), H == inapoi_meniu.

determina(Atr) :-
realizare_scop(av(Atr,_),_,[scop(Atr)]),!.
determina(_).

% ~~~~ sfarsit f) ~~~~

%am modificat modul in care functioneaza afisarea scopului: se creeaza o lista de fapte,
%ce se ordoneaza crescator, si apoi se apeleaza scrie_scop_lista pentru a o afisa in ordine
%inversa -> ordine descrescatoare

%afiseaza_scop(Atr) :-
%nl,fapt(av(Atr,Val),FC,_),
%FC >= 20,scrie_scop(av(Atr,Val),FC),
%nl,fail.

afiseaza_scop(Atr) :-
nl,creaza_lista_fapte(Atr,L),
ordonare_param_fapte(L,Lrez),
scrie_scop_lista(Lrez),
nl,fail.

afiseaza_scop(Atr) :-
nl,\+ fapt(av(Atr,Val),_,_),
write('Nu exista solutii pentru raspunsurile date!'),nl,fail.

afiseaza_scop(_):-nl,nl.

%creaza_lista_fapte(+Atr,-L) - creeaza o lista cu toate faptele

creaza_lista_fapte(Atr,L):- findall(fapt(av(Atr,Val),FC,_),fapt(av(Atr,Val),FC,_),L).

%sort_lista_fapte(+L,-L1) - sorteaza lista in functie de FC

sort_lista_fapte(L,L1):-setof(fapt(FC,av(Atr,Val),Hist), member(fapt(av(Atr,Val),FC,Hist),L),L1).    

%ordonare_param_fapte(+L,-Lrez) - reordoneaza parametrii predicatului (FC si av(...)
ordonare_param_fapte(L,Lrez):-sort_lista_fapte(L,L1),bagof(fapt(av(Atr,Val),FC,Hist), member(fapt(FC,av(Atr,Val),Hist),L1),Lrez).

%scrie_scop_lista(+L) - ia lista de fapte sortata si o afiseaza invers (pt afisarea descrescatoare in functie de FC)
scrie_scop_lista([H|T]):-
scrie_scop_lista(T), scrie_scop(H).
scrie_scop_lista([]).

%scrie_scop_lista_detalii(+L) - acelasi lucru ca si scrie_scop_lista, folosit insa pentru a scrie detaliile pentru fiecare solutie
scrie_scop_lista_detalii([H|T]):-
scrie_scop_lista_detalii(T), scrie_scop_detalii(H).
scrie_scop_lista_detalii([]).

%am adaugat verificarea FC >= 60 in scriere, deoarece am modificat 
%modul in care luam toate faptele
scrie_scop(fapt(av(Atr,Val),FC,_)) :-
FC >= 20,
transformare(av(Atr,Val), X),
scrie_lista(X),tab(2),
write(' '),
write('factorul de certitudine este '),
FC1 is integer(FC),write(FC1), nl.


scrie_scop(av(Atr,Val),FC) :-
transformare(av(Atr,Val), X),
scrie_lista(X),tab(2),
write(' '),
write('factorul de certitudine este '),
FC1 is integer(FC),write(FC1).

% verifica ca e deja in baza de cunostinte
realizare_scop(not Scop,Not_FC,Istorie) :-
realizare_scop(Scop,FC,Istorie),
Not_FC is - FC, !.
realizare_scop(Scop,FC,_) :-
fapt(Scop,FC,_), !.
% adaugat pt nu_conteaza
realizare_scop(av(Atr,_),FC,_) :-
fapt(av(Atr,nu_conteaza),FC,_), !.
realizare_scop(Scop,FC,Istorie) :-
pot_interoga(Scop,Istorie),
!,realizare_scop(Scop,FC,Istorie).
realizare_scop(Scop,FC_curent,Istorie) :-
fg(Scop,FC_curent,Istorie).
        
fg(Scop,FC_curent,Istorie) :-
regula(N, premise(Lista), concluzie(Scop,FC)),
demonstreaza(N,Lista,FC_premise,Istorie),
ajusteaza(FC,FC_premise,FC_nou),
actualizeaza(Scop,FC_nou,FC_curent,N),
FC_curent == 100,!.
fg(Scop,FC,_) :- fapt(Scop,FC,_).

pot_interoga(av(Atr,_),Istorie) :-
not interogat(av(Atr,_)),
interogabil(Atr,Optiuni,Mesaj),
interogheaza(Atr,Mesaj,Optiuni,Istorie),nl,
asserta( interogat(av(Atr,_)) ).

cum([]) :- write('Scop? '),nl,
write('|:'),citeste_linie(Linie),nl,
transformare(Scop,Linie), cum(Scop).
cum(L) :- 
transformare(Scop,L),nl, cum(Scop).
cum(not Scop) :- 
fapt(Scop,FC,Reguli),
lista_float_int(Reguli,Reguli1),
FC < -20,transformare(not Scop,PG),
append(PG,[a,fost,derivat,cu, ajutorul, 'regulilor: '|Reguli1],LL),
scrie_lista(LL),nl,afis_reguli(Reguli),fail.
cum(Scop) :-
fapt(Scop,FC,Reguli),
lista_float_int(Reguli,Reguli1),
FC > 20,transformare(Scop,PG),
append(PG,[a,fost,derivat,cu, ajutorul, 'regulilor: '|Reguli1],LL),
scrie_lista(LL),nl,afis_reguli(Reguli),
fail.
cum(_).

afis_reguli([]).
afis_reguli([N|X]) :-
afis_regula(N),
premisele(N),
afis_reguli(X).
afis_regula(N) :-
regula(N, premise(Lista_premise),
concluzie(Scop,FC)),NN is integer(N),
scrie_lista(['{\nNr_reg~',NN]),
scrie_lista_premise(Lista_premise),
write('concluzie('),
transformare3(Scop,Scop_tr),
FC1 is integer(FC),
append([' [fc '],[FC1],FC2),append(Scop_tr,FC2,LL),append(LL,[']'],LLL),
scrie_lista2(LL),write('])\n}.'),nl.

scrie_lista_premise([]).
scrie_lista_premise([H|T]) :-
write('daca('),
transformare2(H,H_tr),
scrie_lista2(H_tr), 
scrie_lista([')']),
scrie_lista_premise(T).

scrie_lista2([]).
scrie_lista2([H|T]) :-
write(H),
scrie_lista2(T).

transformare3(av(A,da),[A]) :- !.
transformare3(not av(A,da), ['\\+',A]) :- !.
transformare3(av(A,nu),['\\+',A]) :- !.
transformare3(av(A,V),[A,' = ',V]).

transformare2(av(A,da),[A]) :- !.
transformare2(not av(A,da), ['\\+',A]) :- !.
transformare2(av(A,nu),['\\+',A]) :- !.
transformare2(av(A,V),[A,' << ',V]).

transformare(av(A,da),[A]) :- !.
transformare(not av(A,da), [not,A]) :- !.
transformare(av(A,nu),[not,A]) :- !.
transformare(av(A,V),[A,este,V]).

%trans

premisele(N) :-
regula(N, premise(Lista_premise), _),
!, cum_premise(Lista_premise).
        
cum_premise([]).
cum_premise([Scop|X]) :-
cum(Scop),
cum_premise(X).
        
%atribut care pune intrebarea		
interogheaza(Atr,Mesaj,[da,nu],Istorie) :-
!,write(Mesaj),nl, write('da,nu,nu_stiu,nu_conteaza'),
% primea rasp utilizatorului si verifica lista de optiuni
de_la_utiliz(X,Istorie,[da,nu,nu_stiu,nu_conteaza]), 
det_val_fc(X,Val,FC),
asserta( fapt(av(Atr,Val),FC,[utiliz]) ),
% e) incrementarea numarului de aparitii ale unei intrebari ( Atr )
% stergem predicatul dinamic statistica care se unifica cu Atr si incementam cu 1 nr corespunzator
retract(statistica(Atr,Nr)),
Nr1 is Nr + 1,
assert(statistica(Atr,Nr1))
.
interogheaza(Atr,Mesaj,Optiuni,Istorie) :-
write(Mesaj),nl,
% ia optiunile, le afis pe ecran
append(Optiuni,[nu_stiu,nu_conteaza],Optiuni1),
citeste_opt(VLista,Optiuni1,Istorie),
assert_fapt(Atr,VLista),
% e) incrementarea numarului de aparitii ale unei intrebari ( Atr )
% stergem predicatul dinamic statistica care se unifica cu Atr si incementam cu 1 nr corespunzator
retract(statistica(Atr,Nr)),
Nr1 is Nr + 1,
assert(statistica(Atr,Nr1))
.


citeste_opt(X,Optiuni,Istorie) :-
append(['('],Optiuni,Opt1),
append(Opt1,[')'],Opt),
scrie_lista(Opt),
de_la_utiliz(X,Istorie,Optiuni).

de_la_utiliz(X,Istorie,Lista_opt) :-
repeat,write(': '),citeste_linie(X),
proceseaza_raspuns(X,Istorie,Lista_opt).

proceseaza_raspuns([de_ce],Istorie,_) :-nl,afis_istorie(Istorie),!,fail.

proceseaza_raspuns([X],_,Lista_opt):-
member(X,Lista_opt).
proceseaza_raspuns([X,fc,FC],_,Lista_opt):-
member(X,Lista_opt),float(FC).

assert_fapt(Atr,[Val,fc,FC]) :-
!,asserta( fapt(av(Atr,Val),FC,[utiliz]) ).
assert_fapt(Atr,[Val]) :-
asserta( fapt(av(Atr,Val),100,[utiliz])).

det_val_fc([nu],da,-100).
det_val_fc([nu,FC],da,NFC) :- NFC is -FC.
det_val_fc([nu,fc,FC],da,NFC) :- NFC is -FC.
det_val_fc([Val,FC],Val,FC).
det_val_fc([Val,fc,FC],Val,FC).
det_val_fc([Val],Val,100).
        
afis_istorie([]) :- nl.
afis_istorie([scop(X)|T]) :-
scrie_lista([scop,X]),!,
afis_istorie(T).
afis_istorie([N|T]) :-
afis_regula(N),!,afis_istorie(T).

demonstreaza(N,ListaPremise,Val_finala,Istorie) :-
dem(ListaPremise,100,Val_finala,[N|Istorie]),!.

dem([],Val_finala,Val_finala,_).
dem([H|T],Val_actuala,Val_finala,Istorie) :-
realizare_scop(H,FC,Istorie),
Val_interm is min(Val_actuala,FC),
Val_interm >= 20,
dem(T,Val_interm,Val_finala,Istorie).
 
actualizeaza(Scop,FC_nou,FC,RegulaN) :-
fapt(Scop,FC_vechi,_),
combina(FC_nou,FC_vechi,FC),
retract( fapt(Scop,FC_vechi,Reguli_vechi) ),
asserta( fapt(Scop,FC,[RegulaN | Reguli_vechi]) ),!.
actualizeaza(Scop,FC,FC,RegulaN) :-
asserta( fapt(Scop,FC,[RegulaN]) ).

ajusteaza(FC1,FC2,FC) :-
X is FC1 * FC2 / 100,
FC is round(X).
combina(FC1,FC2,FC) :-
FC1 >= 0,FC2 >= 0,
X is FC2*(100 - FC1)/100 + FC1,
FC is round(X).
combina(FC1,FC2,FC) :-
FC1 < 0,FC2 < 0,
X is - ( -FC1 -FC2 * (100 + FC1)/100),
FC is round(X).
combina(FC1,FC2,FC) :-
(FC1 < 0; FC2 < 0),
(FC1 > 0; FC2 > 0),
FCM1 is abs(FC1),FCM2 is abs(FC2),
MFC is min(FCM1,FCM2),
X is 100 * (FC1 + FC2) / (100 - MFC),
FC is round(X).

incarca :-
write('Introduceti numele fisierului cu reguli. '),nl, write('|:'),read(F),
file_exists(F),!,incarca(F).
incarca:-write('Nume incorect de fisier! '),nl,fail.

incarca_info_vehicule :-
write('Introduceti numele fisierului cu informatii. '),nl, write('|:'),read(F),
file_exists(F),!,incarca_info_vehicule(F).
incarca_info_vehicule:-write('Numele fisierului este incorect! '),nl,fail.

%file_exists('info.txt'),!,incarca_info_vehicule('info.txt').

incarca(F) :-
retractall(interogat(_)),retractall(fapt(_,_,_)),
retractall(scop(_)),retractall(interogabil(_,_,_)),
retractall(regula(_,_,_)),
see(F),incarca_reguli,seen, incarca_info_vehicule,!.

incarca_reguli :-
repeat,citeste_propozitie(L),
proceseaza(L),L == [end_of_file],nl.

% ~~~~ parsare b)  ~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%realizam parsarea pentru fisierul care contine informatiile despre vehicule

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% incarca_info_vehicule(+NumeFisierCuInformatii)
% stergem toate predicatele dinamice de forma info_vehicul, deschidem fisierul pt citire, apelam incarca_info, inchidem fisierul pt citire

incarca_info_vehicule(F) :-
retractall(info_vehicul(_,_,_,_)),
see(F),incarca_info,seen,!.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% incarca_info citeste propozitie cu propozitie fisierul
% info_vehicul(Nume_vehicul,Descriere,Pasi,Imagine).
% atom_concat(+Atom1,+Atom2,-AtomiConcatenati).

incarca_info :-
repeat,citeste_descriere(L), 
proceseaza(L),L == [end_of_file],nl.

citeste_descriere(L):- citeste_linie(Linie), 
                      % write(Linie), nl, nl,
					(
					   Linie = [end_of_file], ! , L = [end_of_file] ;
					   Linie =[H|_], H = '*', ! , L = [] ;
					   citeste_descriere(T1), append(Linie, T1, L)
					).

proceseaza([end_of_file]):-!.
proceseaza(L) :-
trad(R,L,[]),assertz(R), !.

trad(scop(X)) --> ['[',scop,']','{',X,'}'].

trad(interogabil(Atr,M,P)) --> 
['[','?',']','{','%',intrebare,'=','['], afiseaza(Atr,P), lista_optiuni(M),['%',atribut_intrebare,':',':','=', Atr, '}'].

trad(regula(N,premise(Daca),concluzie(Atunci,F))) --> identificator(N),daca(Daca),atunci(Atunci,F).

trad(info_vehicul(Nume_vehicul,Descriere,Pasi,Imagine)) --> 
[vehicul, ':','<',Nume_vehicul,'>',descriere,'<',Descriere,'>'],
lista_pasi(Pasi),
[Imagine,'>']
.

trad('Eroare la parsare'-L,L,_).

lista_optiuni(M) --> ['('],lista_de_optiuni(M).
lista_de_optiuni([Element]) -->  [Element,')',']'].
lista_de_optiuni([Element|T]) --> [Element,'^'],lista_de_optiuni(T).

afiseaza(_,P) -->  [P,'/'].
afiseaza(P,P) -->  [].
identificator(N) --> ['{',nr_reg,'~',N].  %regula

daca(Daca) --> [daca,'('],lista_premise(Daca).

lista_premise([Daca]) --> propoz(Daca),[')',concluzie,'('].
lista_premise([Prima|Celalalte]) --> propoz(Prima),[')',daca,'('],lista_premise(Celalalte). % modificat

atunci(av(Atr,Val),FC) --> [Atr,'=',Val,'[',fc, FC,']',')','}'].
atunci(av(Atr,Val),100) --> [Atr,'=',Val,')','}'].

propoz(not av(Atr,da)) --> ['\\','+',Atr].
propoz(av(Atr,Val)) --> [Atr,'<','<',Val]. % am inlocuit este cu <<
propoz(av(Atr,da)) --> [Atr].

lista_pasi(Pasi) --> [pasi,':'], lista_de_pasi(Pasi).

    lista_de_pasi([Element]) -->  
    ['(',Nr,')','<',Element,'>',imagine,'<']
    .
lista_de_pasi([Element|T]) -->
    ['(',Nr,')','<',Element,'>'],
    lista_de_pasi(T)
    .    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~--
% ~~~~ sfarsit parsare b) ~~~~


citeste_linie([Cuv|Lista_cuv]) :-
get_code(Car),
citeste_cuvant(Car, Cuv, Car1), 
rest_cuvinte_linie(Car1, Lista_cuv). 
      
% -1 este codul ASCII pt EOF

rest_cuvinte_linie(-1, []):-!.    
rest_cuvinte_linie(Car,[]) :-(Car==13;Car==10), !.
rest_cuvinte_linie(Car,[Cuv1|Lista_cuv]) :-
citeste_cuvant(Car,Cuv1,Car1),      
rest_cuvinte_linie(Car1,Lista_cuv).

citeste_propozitie([Cuv|Lista_cuv]) :-
get_code(Car),citeste_cuvant(Car, Cuv, Car1), 
rest_cuvinte_propozitie(Car1, Lista_cuv). 
     
rest_cuvinte_propozitie(-1, []):-!.    
rest_cuvinte_propozitie(Car,[]) :-Car==46, !.
rest_cuvinte_propozitie(Car,[Cuv1|Lista_cuv]) :-
citeste_cuvant(Car,Cuv1,Car1),      
rest_cuvinte_propozitie(Car1,Lista_cuv).

citeste_cuvant(-1,end_of_file,-1):-!.
citeste_cuvant(Caracter,Cuvant,Caracter1) :-   
caracter_cuvant(Caracter),!, 
name(Cuvant, [Caracter]),get_code(Caracter1).
citeste_cuvant(Caracter, Numar, Caracter1) :-
caracter_numar(Caracter),!,
citeste_tot_numarul(Caracter, Numar, Caracter1). 

citeste_tot_numarul(Caracter,Numar,Caracter1):-
determina_lista(Lista1,Caracter1),
append([Caracter],Lista1,Lista),
transforma_lista_numar(Lista,Numar).

determina_lista(Lista,Caracter1):-
get_code(Caracter), 
(caracter_numar(Caracter),
determina_lista(Lista1,Caracter1),
append([Caracter],Lista1,Lista); 
\+(caracter_numar(Caracter)),
Lista=[],Caracter1=Caracter). 

transforma_lista_numar([],0).
transforma_lista_numar([H|T],N):-
transforma_lista_numar(T,NN), 
lungime(T,L), Aux is exp(10,L),
HH is H-48,N is HH*Aux+NN.

lungime([],0).
lungime([_|T],L):-
lungime(T,L1),
L is L1+1.

tab(N):-N>0,write(' '), N1 is N-1, tab(N1).
tab(0).

% 39 este codul ASCII pt '


citeste_cuvant(Caracter,Cuvant,Caracter1) :-
Caracter==39,!,
pana_la_urmatorul_apostrof(Lista_caractere),
L=[Caracter|Lista_caractere],
name(Cuvant, L),get_code(Caracter1).        

pana_la_urmatorul_apostrof(Lista_caractere):-
get_code(Caracter),
(Caracter == 39,Lista_caractere=[Caracter];
Caracter\==39,
pana_la_urmatorul_apostrof(Lista_caractere1),
Lista_caractere=[Caracter|Lista_caractere1]).

citeste_cuvant(Caracter,Cuvant,Caracter1) :-          
caractere_in_interiorul_unui_cuvant(Caracter),!,              
((Caracter>64,Caracter<91),!,
Caracter_modificat is Caracter+32;
Caracter_modificat is Caracter),                              
citeste_intreg_cuvantul(Caractere,Caracter1),
name(Cuvant,[Caracter_modificat|Caractere]).        

citeste_intreg_cuvantul(Lista_Caractere,Caracter1) :-
get_code(Caracter),
(caractere_in_interiorul_unui_cuvant(Caracter),
((Caracter>64,Caracter<91),!, 
Caracter_modificat is Caracter+32;
Caracter_modificat is Caracter),
citeste_intreg_cuvantul(Lista_Caractere1, Caracter1),
Lista_Caractere=[Caracter_modificat|Lista_Caractere1]; \+(caractere_in_interiorul_unui_cuvant(Caracter)),
Lista_Caractere=[], Caracter1=Caracter).

citeste_cuvant(_,Cuvant,Caracter1) :-                
get_code(Caracter),       
citeste_cuvant(Caracter,Cuvant,Caracter1). 

caracter_cuvant(C):-member(C,[44,59,58,63,33,46,41,40,126,94,91,93,123,125,60,37,61,92,43,47,62,42]).

% am specificat codurile ASCII pentru , ; : ? ! . ~ ^ ) ( [ ] { } < % = \ + / > *

caractere_in_interiorul_unui_cuvant(C):-
C>64,C<91;C>47,C<58;
C==45;C==95;C>96,C<123.
caracter_numar(C):-C<58,C>=48.


