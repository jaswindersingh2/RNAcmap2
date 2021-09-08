# RNAcmap2
An improved fully automatic pipeline for predicting contact maps of RNAs by evolutionary coupling analysis


1. `conda activate venv_rnacmap2`

2. `./run.py --neff low --dca_method gremlin`


```
		       		  F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.007 	 0.008 		 0.007		    74		  0.0    
direct_infernal     0.168 	 0.216 		 0.139		    74		  2.0    
RNAcmap             0.187 	 0.242 		 0.155		    74		  3.0    
RNAcmap2            0.335 	 0.439 		 0.275		    74		  11.4    
RNAcmap_meta        0.222 	 0.287 		 0.186		    74		  5.0    
RNAcmap2_meta       0.416 	 0.542 		 0.344		    74		  15.8
```

3. `./run.py --neff median --dca_method gremlin`


```
			          F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.036 	 0.047 		 0.029		    31		  1.0    
direct_infernal     0.357 	 0.486 		 0.284		    31		  17.0    
RNAcmap             0.411 	 0.554 		 0.334		    31		  26.5    
RNAcmap2            0.515 	 0.678 		 0.443		    31		  105.0    
RNAcmap_meta        0.461 	 0.617 		 0.383		    31		  31.6    
RNAcmap2_meta       0.532 	 0.708 		 0.436		    31		  96.1
```

4. `./run.py --neff high --dca_method gremlin`


```
		       		  F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.074 	 0.104 		 0.058		    114		  2.2    
direct_infernal     0.589 	 0.805 		 0.467		    114		  321.4    
RNAcmap             0.633 	 0.862 		 0.502		    114		  636.5    
RNAcmap2            0.633 	 0.862 		 0.502		    114		  636.5    
RNAcmap_meta        0.630 	 0.857 		 0.501		    114		  605.1    
RNAcmap2_meta       0.630 	 0.857 		 0.501		    114		  605.1
```

5. `./run.py --neff all --dca_method gremlin`


```
		       		  F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.046 	 0.063 		 0.037		    219		  1.0    
direct_infernal     0.414 	 0.561 		 0.331		    219		  27.3    
RNAcmap             0.451 	 0.609 		 0.361		    219		  58.4    
RNAcmap2            0.515 	 0.693 		 0.417		    219		  169.0    
RNAcmap_meta        0.468 	 0.631 		 0.378		    219		  77.1    
RNAcmap2_meta       0.544 	 0.730 		 0.439		    219		  197.4
```

6. `./run.py --neff low --dca_method plmc`


```
		       		  F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.013 	 0.016 		 0.011		    74		  1.0    
direct_infernal     0.000 	 0.000 		 0.000		    74		  0.0    
RNAcmap             0.172 	 0.225 		 0.142		    74		  3.0    
RNAcmap2            0.351 	 0.459 		 0.289		    74		  11.4    
RNAcmap_meta        0.229 	 0.298 		 0.191		    74		  5.0    
RNAcmap2_meta       0.421 	 0.545 		 0.349		    74		  15.8
```

7. `./run.py --neff low --dca_method mfdca`


```
		       		  F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.017 	 0.022 		 0.014		    74		  1.0    
direct_infernal     0.000 	 0.000 		 0.000		    74		  0.0    
RNAcmap             0.203 	 0.266 		 0.168		    74		  3.0    
RNAcmap2            0.363 	 0.475 		 0.298		    74		  11.4    
RNAcmap_meta        0.247 	 0.324 		 0.204		    74		  5.0    
RNAcmap2_meta       0.451 	 0.588 		 0.373		    74		  15.8
```

8. `./run.py --neff low --dca_method plmdca`


```
		       		  F1	Precision	Sensitivity	No. of RNAs	Median Neff
blastn              0.009 	 0.012 		 0.008		    74		  1.0    
direct_infernal     0.000 	 0.000 		 0.000		    74		  0.0    
RNAcmap             0.191 	 0.249 		 0.158		    74		  3.0    
RNAcmap2            0.359 	 0.471 		 0.295		    74		  11.4    
RNAcmap_meta        0.232 	 0.303 		 0.192		    74		  5.0    
RNAcmap2_meta       0.439 	 0.570 		 0.363		    74		  15.8
```
