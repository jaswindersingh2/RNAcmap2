#Name: PKB68
#Length:  68 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GUUAAACCAUCACUGUGUACUGGCCAAUAAACCCAGUAUCCCGUUCGACCGAGCCCGCUUAUCGGUAU
..[{(.]....[.}....(.................)(.....]..........)(....))......
EEEESMMMMMMMMMMMMMSHHHHHHHHHHHHHHHHHSSHHHHHHHHHHHHHHHHSSHHHHSSEEEEEE
NNKKNNKNNNNKNKNNNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNNNNNN
S1 5..5 "A" 62..62 "U"
S2 19..19 "A" 37..37 "U"
S3 38..38 "A" 55..55 "C"
S4 56..56 "C" 61..61 "A"
H1 20..36 "CUGGCCAAUAAACCCAG" (19,37) A:U 
H2 39..54 "UCCCGUUCGACCGAGC" (38,55) A:C PK{3}
H3 57..60 "GCUU" (56,61) C:A 
M1.1 6..18 "ACCAUCACUGUGU" (5,62) A:U (19,37) A:U PK{1,2,3}
M1.2 38..37 "" (37,19) U:A (38,55) A:C 
M1.3 56..55 "" (55,38) C:A (56,61) C:A 
M1.4 62..61 "" (61,56) A:C (62,5) U:A 
E1 1..4 "GUUA" PK{1,2}
E2 63..68 "CGGUAU" 
PK1 1bp 3..3 7..7 E1 1..4 M1.1 6..18
PK2 1bp 4..4 14..14 E1 1..4 M1.1 6..18
PK3 1bp 12..12 44..44 M1.1 6..18 H2 39..54
PK1.1 3 U 7 C
PK2.1 4 A 14 U
PK3.1 12 A 44 U
NCBP1 38 A 55 C S3
NCBP2 56 C 61 A S4
NCBP3 3 U 7 C PK1.1
segment1 1bp 5..5 A 62..62 U
segment2 1bp 19..19 A 37..37 U
segment3 1bp 38..38 A 55..55 C
segment4 1bp 56..56 C 61..61 A
