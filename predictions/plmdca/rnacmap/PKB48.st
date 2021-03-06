#Name: PKB48
#Length:  61 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
UAGGGCCGUCAGGGUUCUCCCGCCCUCCGUGAACCCAGGCUAAAAGUUAAGGUAGGGGGGC
.[........(....[...(]{......)(.........(...)...).]......)...}
EEEEEEEEEESMMMMMMMMSHHHHHHHHSSIIIIIIIIISHHHSIIISMMMMMMMMSEEEE
NKNNNNNNNNNNNNNKNNNNKKNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNNNNNNNNK
S1 11..11 "A" 57..57 "G"
S2 20..20 "C" 29..29 "G"
S3 30..30 "U" 48..48 "U"
S4 40..40 "C" 44..44 "A"
H1 21..28 "CGCCCUCC" (20,29) C:G PK{2,3}
H2 41..43 "UAA" (40,44) C:A 
I1.1 31..39 "GAACCCAGG" (30,48) U:U 
I1.2 45..47 "AGU" (44,40) A:C 
M1.1 12..19 "GGGUUCUC" (11,57) A:G (20,29) C:G PK{2}
M1.2 30..29 "" (29,20) G:C (30,48) U:U 
M1.3 49..56 "AAGGUAGG" (48,30) U:U (57,11) G:A PK{1}
E1 1..10 "UAGGGCCGUC" PK{1}
E2 58..61 "GGGC" PK{3}
PK1 1bp 2..2 50..50 E1 1..10 M1.3 49..56
PK2 1bp 16..16 21..21 M1.1 12..19 H1 21..28
PK3 1bp 22..22 61..61 H1 21..28 E2 58..61
PK1.1 2 A 50 A
PK2.1 16 U 21 C
PK3.1 22 G 61 C
NCBP1 30 U 48 U S3
NCBP2 40 C 44 A S4
NCBP3 11 A 57 G S1
NCBP4 2 A 50 A PK1.1
NCBP5 16 U 21 C PK2.1
segment1 1bp 11..11 A 57..57 G
segment2 1bp 20..20 C 29..29 G
segment3 2bp 30..40 UGAACCCAGGC 44..48 AAGUU
