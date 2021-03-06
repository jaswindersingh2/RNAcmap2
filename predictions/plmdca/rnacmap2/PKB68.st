#Name: PKB68
#Length:  68 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GUUAAACCAUCACUGUGUACUGGCCAAUAAACCCAGUAUCCCGUUCGACCGAGCCCGCUUAUCGGUAU
..(([.[....{.)....(...(.......].)...))[....}.......(..](....)).....]
EESSHHHHHHHHHSMMMMSIIISHHHHHHHHHSIIISSXXXXXXXXXXXXXSBBBSHHHHSSEEEEEE
NNNNKNKNNNNKNNNNNNNNNNNNNNNNNNKNNNNNNNKNNNNKNNNNNNNNNNKNNNNNNNNNNNNK
S1 3..3 "U" 38..38 "A"
S2 4..4 "A" 14..14 "U"
S3 19..19 "A" 37..37 "U"
S4 23..23 "G" 33..33 "C"
S5 52..52 "A" 62..62 "U"
S6 56..56 "C" 61..61 "A"
H1 5..13 "AACCAUCAC" (4,14) A:U PK{1,2,3}
H2 24..32 "CCAAUAAAC" (23,33) G:C PK{2}
H3 57..60 "GCUU" (56,61) C:A 
B1 53..55 "GCC" (52,62) A:U (56,61) C:A PK{4}
I1.1 20..22 "CUG" (19,37) A:U 
I1.2 34..36 "CAG" (33,23) C:G 
M1.1 4..3 "" (3,38) U:A (4,14) A:U 
M1.2 15..18 "GUGU" (14,4) U:A (19,37) A:U 
M1.3 38..37 "" (37,19) U:A (38,3) A:U 
X1 39..51 "UCCCGUUCGACCG" (38,3) A:U (52,62) A:U PK{3,4}
E1 1..2 "GU" 
E2 63..68 "CGGUAU" PK{1}
PK1 1bp 5..5 68..68 H1 5..13 E2 63..68
PK2 1bp 7..7 31..31 H1 5..13 H2 24..32
PK3 1bp 12..12 44..44 H1 5..13 X1 39..51
PK4 1bp 39..39 55..55 X1 39..51 B1 53..55
PK1.1 5 A 68 U
PK2.1 7 C 31 A
PK3.1 12 A 44 U
PK4.1 39 U 55 C
NCBP1 56 C 61 A S6
NCBP2 7 C 31 A PK2.1
NCBP3 39 U 55 C PK4.1
segment1 1bp 3..3 U 38..38 A
segment2 1bp 4..4 A 14..14 U
segment3 2bp 19..23 ACUGG 33..37 CCAGU
segment4 2bp 52..56 AGCCC 61..62 AU
