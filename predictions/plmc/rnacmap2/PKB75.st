#Name: PKB75
#Length:  60 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
UGGCCGGCAUGGUCCCAGCCUCCUCGCUGGCGCCGGCUGGGCAACAUUCCGAGGGGACCG
..[..(....])...............(.........).......[.((.....)])...
EEEEESHHHHHSMMMMMMMMMMMMMMMSHHHHHHHHHSMMMMMMMMMSSHHHHHSBSEEE
NNKNNNNNNNKNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNNNNNNNKNNNN
S1 6..6 "G" 12..12 "G"
S2 28..28 "U" 38..38 "U"
S3 48..48 "U" 57..57 "A"
S4 49..49 "C" 55..55 "G"
H1 7..11 "GCAUG" (6,12) G:G PK{1}
H2 29..37 "GGCGCCGGC" (28,38) U:U 
H3 50..54 "CGAGG" (49,55) C:G 
B1 56..56 "G" (55,49) G:C (57,48) A:U PK{2}
M1.1 13..27 "UCCCAGCCUCCUCGC" (12,6) G:G (28,38) U:U 
M1.2 39..47 "GGGCAACAU" (38,28) U:U (48,57) U:A PK{2}
E1 1..5 "UGGCC" PK{1}
E2 58..60 "CCG" 
PK1 1bp 3..3 11..11 E1 1..5 H1 7..11
PK2 1bp 46..46 56..56 M1.2 39..47 B1 56..56
PK1.1 3 G 11 G
PK2.1 46 A 56 G
NCBP1 28 U 38 U S2
NCBP2 6 G 12 G S1
NCBP3 3 G 11 G PK1.1
NCBP4 46 A 56 G PK2.1
segment1 1bp 6..6 G 12..12 G
segment2 1bp 28..28 U 38..38 U
segment3 2bp 48..49 UC 55..57 GGA
