#Name: 5vof_A
#Length: ,36,
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GAGAGCCCCAGCGAGAUAAUACUUGGCCCCGCUCUU
[.(((...[...(....].....)....]..)))..
EESSSIIIIIIISHHHHHHHHHHSIIIIIIISSSEE
KNNNNNNNKNNNNNNNNKNNNNNNNNNNKNNNNNNN
S1 3..5 "GAG" 32..34 "CUC"
S2 13..13 "G" 24..24 "U"
H1 14..23 "AGAUAAUACU" (13,24) G:U PK{2}
I1.1 6..12 "CCCCAGC" (5,32) G:C PK{2}
I1.2 25..31 "GGCCCCG" (24,13) U:G PK{1}
E1 1..2 "GA" PK{1}
E2 35..36 "UU" 
PK1 1bp 1..1 29..29 E1 1..2 I1.2 25..31
PK2 1bp 9..9 18..18 I1.1 6..12 H1 14..23
PK1.1 1 G 29 C
PK2.1 9 C 18 A
NCBP1 9 C 18 A PK2.1
segment1 4bp 3..13 GAGCCCCAGCG 24..34 UGGCCCCGCUC
