#Name: PKB21
#Length:  38 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGGAUUCGAAUUCCCCCGGAACCCCCGGUUGGGGCCCA
.(((.......))).........[(......).(.].)
ESSSHHHHHHHSSSMMMMMMMMMMSHHHHHHSMSHHHS
NNNNNNNNNNNNNNNNNNNNNNNKNNNNNNNNNNNKNN
S1 2..4 "GGA" 12..14 "UCC"
S2 25..25 "C" 32..32 "G"
S3 34..34 "G" 38..38 "A"
H1 5..11 "UUCGAAU" (4,12) A:U 
H2 26..31 "CGGUUG" (25,32) C:G 
H3 35..37 "CCC" (34,38) G:A PK{1}
M1.1 15..24 "CCCGGAACCC" (14,2) C:G (25,32) C:G PK{1}
M1.2 33..33 "G" (32,25) G:C (34,38) G:A 
E1 1..1 "G" 
PK1 1bp 24..24 36..36 M1.1 15..24 H3 35..37
PK1.1 24 C 36 C
NCBP1 34 G 38 A S3
NCBP2 24 C 36 C PK1.1
segment1 3bp 2..4 GGA 12..14 UCC
segment2 1bp 25..25 C 32..32 G
segment3 1bp 34..34 G 38..38 A
