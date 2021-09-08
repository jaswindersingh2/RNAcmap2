#Name: PKB64
#Length:  32 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
AAAUUGAAGAGUUUGAUCAUGGCUCAAUGAAU
[........(.............).(....])
EEEEEEEEESHHHHHHHHHHHHHSXSHHHHHS
KNNNNNNNNNNNNNNNNNNNNNNNNNNNNNKN
S1 10..10 "A" 24..24 "U"
S2 26..26 "A" 32..32 "U"
H1 11..23 "GUUUGAUCAUGGC" (10,24) A:U 
H2 27..31 "AUGAA" (26,32) A:U PK{1}
X1 25..25 "C" (24,10) U:A (26,32) A:U 
E1 1..9 "AAAUUGAAG" PK{1}
PK1 1bp 1..1 31..31 E1 1..9 H2 27..31
PK1.1 1 A 31 A
NCBP1 1 A 31 A PK1.1
segment1 1bp 10..10 A 24..24 U
segment2 1bp 26..26 A 32..32 U
