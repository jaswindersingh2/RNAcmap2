#Name: PKB202
#Length: ,34,
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
UGGUGCUUGUUAUUUCACCUAAAUCGAAAUAACG
[{......(]..(...............)..}).
EEEEEEEESIIISHHHHHHHHHHHHHHHSIIISE
KKNNNNNNNKNNNNNNNNNNNNNNNNNNNNNKNN
S1 9..9 "G" 33..33 "C"
S2 13..13 "U" 29..29 "A"
H1 14..28 "UUCACCUAAAUCGAA" (13,29) U:A 
I1.1 10..12 "UUA" (9,33) G:C PK{1}
I1.2 30..32 "UAA" (29,13) A:U PK{2}
E1 1..8 "UGGUGCUU" PK{1,2}
E2 34..34 "G" 
PK1 1bp 1..1 10..10 E1 1..8 I1.1 10..12
PK2 1bp 2..2 32..32 E1 1..8 I1.2 30..32
PK1.1 1 U 10 U
PK2.1 2 G 32 A
NCBP1 1 U 10 U PK1.1
NCBP2 2 G 32 A PK2.1
segment1 2bp 9..13 GUUAU 29..33 AUAAC
