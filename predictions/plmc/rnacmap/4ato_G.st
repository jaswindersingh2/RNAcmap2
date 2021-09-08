#Name: 4ato_G
#Length: ,34,
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
AAAUUGGUGUAACCUUACCGUAGUAGGUGCUAAA
.......[...(.{(...].....}).)......
EEEEEEEEEEESIISHHHHHHHHHHSISEEEEEE
NNNNNNNKNNNNNKNNNNKNNNNNKNNNNNNNNN
S1 12..12 "A" 28..28 "U"
S2 15..15 "U" 26..26 "G"
H1 16..25 "UACCGUAGUA" (15,26) U:G PK{1,2}
I1.1 13..14 "CC" (12,28) A:U PK{2}
I1.2 27..27 "G" (26,15) G:U 
E1 1..11 "AAAUUGGUGUA" PK{1}
E2 29..34 "GCUAAA" 
PK1 1bp 8..8 19..19 E1 1..11 H1 16..25
PK2 1bp 14..14 25..25 I1.1 13..14 H1 16..25
PK1.1 8 U 19 C
PK2.1 14 C 25 A
NCBP1 8 U 19 C PK1.1
NCBP2 14 C 25 A PK2.1
segment1 2bp 12..15 ACCU 26..28 GGU
