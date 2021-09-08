#Name: PKB202
#Length: ,34,
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
UGGUGCUUGUUAUUUCACCUAAAUCGAAAUAACG
[{......(](.(.......(...)...).)}).
EEEEEEEESISISIIIIIIISHHHSIIISISISE
KKNNNNNNNKNNNNNNNNNNNNNNNNNNNNNKNN
S1 9..9 "G" 33..33 "C"
S2 11..11 "U" 31..31 "A"
S3 13..13 "U" 29..29 "A"
S4 21..21 "A" 25..25 "C"
H1 22..24 "AAU" (21,25) A:C 
I1.1 10..10 "U" (9,33) G:C PK{1}
I1.2 32..32 "A" (31,11) A:U PK{2}
I2.1 12..12 "A" (11,31) U:A 
I2.2 30..30 "U" (29,13) A:U 
I3.1 14..20 "UUCACCU" (13,29) U:A 
I3.2 26..28 "GAA" (25,21) C:A 
E1 1..8 "UGGUGCUU" PK{1,2}
E2 34..34 "G" 
PK1 1bp 1..1 10..10 E1 1..8 I1.1 10..10
PK2 1bp 2..2 32..32 E1 1..8 I1.2 32..32
PK1.1 1 U 10 U
PK2.1 2 G 32 A
NCBP1 21 A 25 C S4
NCBP2 1 U 10 U PK1.1
NCBP3 2 G 32 A PK2.1
segment1 4bp 9..21 GUUAUUUCACCUA 25..33 CGAAAUAAC
