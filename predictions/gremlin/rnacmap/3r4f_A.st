#Name: 3r4f_A
#Length:  66 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGCGUGUAUGUGGGGAAACCCUGAUUGAGUUCAGCCCACAUACUUUGUUGAUGGUUGUCAAUCGCC
[({..[{(................................)(]..}.)...............])}
ESMMMMMSHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHSSHHHHHSMMMMMMMMMMMMMMMMSE
KNKNNKKNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNKNNNNNNNNNNNNNNNNNKNK
S1 2..2 "G" 65..65 "C"
S2 8..8 "A" 41..41 "U"
S3 42..42 "A" 48..48 "U"
H1 9..40 "UGUGGGGAAACCCUGAUUGAGUUCAGCCCACA" (8,41) A:U 
H2 43..47 "CUUUG" (42,48) A:U PK{3,4}
M1.1 3..7 "CGUGU" (2,65) G:C (8,41) A:U PK{2,3,4}
M1.2 42..41 "" (41,8) U:A (42,48) A:U 
M1.3 49..64 "UGAUGGUUGUCAAUCG" (48,42) U:A (65,2) C:G PK{1}
E1 1..1 "G" PK{1}
E2 66..66 "C" PK{2}
PK1 1bp 1..1 64..64 E1 1..1 M1.3 49..64
PK2 1bp 3..3 66..66 M1.1 3..7 E2 66..66
PK3 1bp 6..6 43..43 M1.1 3..7 H2 43..47
PK4 1bp 7..7 46..46 M1.1 3..7 H2 43..47
PK1.1 1 G 64 G
PK2.1 3 C 66 C
PK3.1 6 G 43 C
PK4.1 7 U 46 U
NCBP1 1 G 64 G PK1.1
NCBP2 3 C 66 C PK2.1
NCBP3 7 U 46 U PK4.1
segment1 1bp 2..2 G 65..65 C
segment2 1bp 8..8 A 41..41 U
segment3 1bp 42..42 A 48..48 U
