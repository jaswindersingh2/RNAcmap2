#Name: PKB194
#Length: ,28,
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
UACGCUGUACAGUGCGUUAAACUGUACA
...[{(....](......}..).....)
EEEEESIIIIISHHHHHHHHHSIIIIIS
NNNKKNNNNNKNNNNNNNKNNNNNNNNN
S1 6..6 "U" 28..28 "A"
S2 12..12 "G" 22..22 "C"
H1 13..21 "UGCGUUAAA" (12,22) G:C PK{2}
I1.1 7..11 "GUACA" (6,28) U:A PK{1}
I1.2 23..27 "UGUAC" (22,12) C:G 
E1 1..5 "UACGC" PK{1,2}
PK1 1bp 4..4 11..11 E1 1..5 I1.1 7..11
PK2 1bp 5..5 19..19 E1 1..5 H1 13..21
PK1.1 4 G 11 A
PK2.1 5 C 19 A
NCBP1 4 G 11 A PK1.1
NCBP2 5 C 19 A PK2.1
segment1 2bp 6..12 UGUACAG 22..28 CUGUACA
